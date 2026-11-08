! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_tracer_transport
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: block_type,mpas_pool_type,MPAS_LOG_CRIT
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_dimension
 use mpas_log,only: mpas_log_write
 use mpas_atmphys_tracer_column,only: transport_budget_t,budget_reset,budget_accumulate,column_dry_mass
 use mpas_atmphys_pbl_gas,only: pbl_gas_init,pbl_gas_enabled,pbl_gas_active,pbl_gas_apply,pbl_gas_finalize
 use mpas_atmphys_cu_gas ,only: cu_gas_init ,cu_gas_enabled ,cu_gas_active ,cu_gas_apply ,cu_gas_finalize

 implicit none
 private
 public:: TRANSPORT_PBL,                     &
          TRANSPORT_CU,                      &
          transport_budget_t,                &
          tracer_transport_init,             &
          tracer_transport_assemble,         &
          tracer_transport_finalize,         &
          tracer_transport_setup,            &
          tracer_transport_teardown,         &
          tracer_transport_n,                &
          tracer_transport_enabled,          &
          tracer_transport_buffer_size,      &
          tracer_transport_assemble_columns, &
          tracer_transport_add_tend,         &
          tracer_transport_budget


!Vertical transport of independently transported MICM gases by PBL mixing and convection (shared foundation,
!P1 of docs/plans/pbl-transport.md).
!
!Once per model step, after chemistry has updated scalars(:,:,:,1), each owned column runs the deterministic
!sequence q1 = PBL(q0,dt), q2 = CU(q1,dt) on scratch copies, giving R_pbl = (q1-q0)/dt and R_cu = (q2-q1)/dt.
!physics_get_tend then adds (R_pbl+R_cu)*rho_zz once at each mapped scalar index. The prognostic state is
!never written here: the normal scalar advance is the only owner of gas state mutation.
!
!The map (gas n -> MPAS scalar index) is built by chemistry and reaches physics through block%dimensions, so
!libphys does not depend on libchem. One block per MPI task is assumed, as in chemistry.


 integer,parameter:: TRANSPORT_PBL = 1
 integer,parameter:: TRANSPORT_CU  = 2

 integer,save:: nTransport  = 0
 integer,save:: nVertLevels = 0
 integer,save:: nCellsSolve = 0
 integer,dimension(:),allocatable,save:: transport_to_scalar           ! (nTransport)

 real(kind=RKIND),dimension(:,:,:),allocatable,save:: r_pbl             ! (nTransport,nVertLevels,nCellsSolve) [s-1]
 real(kind=RKIND),dimension(:,:,:),allocatable,save:: r_cu              ! (nTransport,nVertLevels,nCellsSolve) [s-1]
 type(transport_budget_t),dimension(:),allocatable,save:: bud_pbl,bud_cu ! (nTransport), this rank, last assembly


 contains


!=================================================================================================================
 subroutine tracer_transport_init(block)
!=================================================================================================================
!read the map published by chemistry_init, initialize the transport processes, then allocate the tendency
!buffers if any process is enabled. Call once, after chemistry_init. Without chemistry (MUSICA=false or no MICM
!species) the map is absent or empty and every transport call is a no-op.
 type(block_type),intent(in):: block

 integer,pointer:: n_ptr,nVertLevels_ptr,nCellsSolve_ptr
 integer,dimension(:),pointer:: idx_ptr
 integer:: n

 nullify(n_ptr,idx_ptr,nVertLevels_ptr,nCellsSolve_ptr)
 call mpas_pool_get_dimension(block%dimensions,'nMicmTransport',n_ptr)
 call mpas_pool_get_dimension(block%dimensions,'nVertLevels'   ,nVertLevels_ptr)
 call mpas_pool_get_dimension(block%dimensions,'nCellsSolve'   ,nCellsSolve_ptr)

 n = 0
 if(associated(n_ptr)) n = n_ptr
 if(n > 0) then
    if(associated(block%next)) &
       call mpas_log_write('[transport] MICM gas transport requires one block per MPI task',messageType=MPAS_LOG_CRIT)
    call mpas_pool_get_dimension(block%dimensions,'micmTransportIndex',idx_ptr)
    if(.not.associated(idx_ptr)) &
       call mpas_log_write('[transport] nMicmTransport is set but micmTransportIndex is missing', &
                           messageType=MPAS_LOG_CRIT)
 endif

!processes first: setup sizes the tendency buffers by whether any process is enabled.
 call pbl_gas_init(block,n)
 call cu_gas_init (block,n)

 if(n > 0) then
    call tracer_transport_setup(idx_ptr(1:n),nVertLevels_ptr,nCellsSolve_ptr)
 else
    call tracer_transport_setup([integer::],nVertLevels_ptr,nCellsSolve_ptr)
 endif
 call mpas_log_write('[transport] nMicmTransport = $i; PBL gas transport enabled = $l; '// &
                     'convective gas transport enabled = $l; tendency buffer size = $i', &
                     intArgs=(/n,tracer_transport_buffer_size()/), &
                     logicArgs=(/pbl_gas_enabled(),cu_gas_enabled()/))

 end subroutine tracer_transport_init

!=================================================================================================================
 subroutine tracer_transport_assemble(mesh,state,dt)
!=================================================================================================================
!once per model step from physics_get_tend, after chemistry_step has updated scalars(:,:,:,1).
 type(mpas_pool_type),intent(in):: mesh,state
 real(kind=RKIND),intent(in):: dt

 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: zz,rho_zz,zgrid
 real(kind=RKIND),dimension(:),pointer:: areaCell

 if(.not.tracer_transport_enabled()) return

 call mpas_pool_get_array(state,'scalars' ,scalars,1)    ! post-chemistry, authoritative
 call mpas_pool_get_array(state,'rho_zz'  ,rho_zz ,2)    ! the tendency mass; equal to level 1 here, and the
                                                          ! level refreshed on the host under MPAS_OPENACC
 call mpas_pool_get_array(mesh ,'zz'      ,zz      )
 call mpas_pool_get_array(mesh ,'zgrid'   ,zgrid   )
 call mpas_pool_get_array(mesh ,'areaCell',areaCell)
 call tracer_transport_assemble_columns(scalars,zz,rho_zz,zgrid,areaCell,dt)

 end subroutine tracer_transport_assemble

!=================================================================================================================
 subroutine tracer_transport_finalize()
!=================================================================================================================

 call tracer_transport_teardown()
 call pbl_gas_finalize()
 call cu_gas_finalize()

 end subroutine tracer_transport_finalize

!=================================================================================================================
 subroutine tracer_transport_setup(map,nlevels,ncells)
!=================================================================================================================
!store the map and allocate the budget records, outside any threaded region. The (nTransport,nVertLevels,
!nCellsSolve) tendency buffers are allocated only when a process is enabled, so disabled transport costs no
!memory. Initialize the processes before calling this.
 integer,intent(in),dimension(:):: map
 integer,intent(in):: nlevels,ncells

 call tracer_transport_teardown()

 nTransport  = size(map)
 nVertLevels = nlevels
 nCellsSolve = ncells

 allocate(transport_to_scalar(nTransport))
 transport_to_scalar(:) = map(:)

 allocate(bud_pbl(nTransport),bud_cu(nTransport))
 if(tracer_transport_enabled()) then
    allocate(r_pbl(nTransport,nVertLevels,nCellsSolve))
    allocate(r_cu (nTransport,nVertLevels,nCellsSolve))
    r_pbl(:,:,:) = 0._RKIND
    r_cu (:,:,:) = 0._RKIND
 endif

 end subroutine tracer_transport_setup

!=================================================================================================================
 subroutine tracer_transport_teardown()
!=================================================================================================================
!release every buffer; a later setup starts clean.

 if(allocated(transport_to_scalar)) deallocate(transport_to_scalar)
 if(allocated(r_pbl)  ) deallocate(r_pbl  )
 if(allocated(r_cu)   ) deallocate(r_cu   )
 if(allocated(bud_pbl)) deallocate(bud_pbl)
 if(allocated(bud_cu) ) deallocate(bud_cu )
 nTransport  = 0
 nVertLevels = 0
 nCellsSolve = 0

 end subroutine tracer_transport_teardown

!=================================================================================================================
 integer function tracer_transport_n()
!=================================================================================================================

 tracer_transport_n = nTransport

 end function tracer_transport_n

!=================================================================================================================
 logical function tracer_transport_enabled()
!=================================================================================================================
!false when there is nothing to transport or no process can act; the caller then skips all transport work and
!leaves the scalar tendencies untouched.

 tracer_transport_enabled = nTransport > 0 .and. (pbl_gas_enabled() .or. cu_gas_enabled())

 end function tracer_transport_enabled

!=================================================================================================================
 integer function tracer_transport_buffer_size()
!=================================================================================================================
!number of values held in the PBL and convective tendency buffers (0 when transport is disabled).

 tracer_transport_buffer_size = 0
 if(allocated(r_pbl)) tracer_transport_buffer_size = tracer_transport_buffer_size + size(r_pbl)
 if(allocated(r_cu) ) tracer_transport_buffer_size = tracer_transport_buffer_size + size(r_cu)

 end function tracer_transport_buffer_size

!=================================================================================================================
 subroutine tracer_transport_assemble_columns(scalars,zz,rho_zz,zgrid,areaCell,dt)
!=================================================================================================================
!build R_pbl and R_cu for every owned column from the current (post-chemistry) gas state, and the per-rank
!budgets of this assembly. Columns where no process is active keep exactly zero tendencies.
 real(kind=RKIND),intent(in),dimension(:,:,:):: scalars       ! (num_scalars,nVertLevels,nCells), not modified
 real(kind=RKIND),intent(in),dimension(:,:):: zz,rho_zz       ! (nVertLevels,nCells)
 real(kind=RKIND),intent(in),dimension(:,:):: zgrid           ! (nVertLevels+1,nCells)
 real(kind=RKIND),intent(in),dimension(:):: areaCell          ! (nCells) [m2]
 real(kind=RKIND),intent(in):: dt                             ! model time step [s]

 real(kind=RKIND),dimension(:,:),allocatable:: q0,q1,q2
 real(kind=RKIND),dimension(nVertLevels):: m_d
 logical:: pbl_on,cu_on,pbl_any,cu_any
 integer:: iCell,k,n

 do n = 1, nTransport
    call budget_reset(bud_pbl(n))
    call budget_reset(bud_cu (n))
 enddo
 if(.not.tracer_transport_enabled()) return

 r_pbl(:,:,:) = 0._RKIND
 r_cu (:,:,:) = 0._RKIND

 pbl_any = pbl_gas_enabled()
 cu_any  = cu_gas_enabled()
 allocate(q0(nTransport,nVertLevels),q1(nTransport,nVertLevels),q2(nTransport,nVertLevels))

 do iCell = 1, nCellsSolve
    pbl_on = .false.
    cu_on  = .false.
    if(pbl_any) pbl_on = pbl_gas_active(iCell)
    if(cu_any ) cu_on  = cu_gas_active(iCell)
    if(.not.(pbl_on .or. cu_on)) cycle

    m_d(:) = column_dry_mass(zz(:,iCell),rho_zz(:,iCell),zgrid(:,iCell))
    do k = 1, nVertLevels
    do n = 1, nTransport
       q0(n,k) = scalars(transport_to_scalar(n),k,iCell)
    enddo
    enddo

    q1(:,:) = q0(:,:)
    if(pbl_on) call pbl_gas_apply(iCell,dt,m_d,q1)
    q2(:,:) = q1(:,:)
    if(cu_on ) call cu_gas_apply (iCell,dt,m_d,q2)

    r_pbl(:,:,iCell) = (q1(:,:)-q0(:,:))/dt
    r_cu (:,:,iCell) = (q2(:,:)-q1(:,:))/dt

    do n = 1, nTransport
       if(pbl_on) call budget_accumulate(bud_pbl(n),areaCell(iCell),m_d,q0(n,:),q1(n,:),dt)
       if(cu_on ) call budget_accumulate(bud_cu (n),areaCell(iCell),m_d,q1(n,:),q2(n,:),dt)
    enddo
 enddo

 deallocate(q0,q1,q2)

 end subroutine tracer_transport_assemble_columns

!=================================================================================================================
 subroutine tracer_transport_add_tend(mass,tend_scalars)
!=================================================================================================================
!add the PBL and convective gas tendencies exactly once per mapped scalar, after tend_scalars was cleared.
!Native moisture keeps its own path: host-bound species such as H2O are never in the map.
 real(kind=RKIND),intent(in),dimension(:,:):: mass               ! rho_zz at time level 2, as for native moisture
 real(kind=RKIND),intent(inout),dimension(:,:,:):: tend_scalars
 integer:: i,k,n

 if(.not.tracer_transport_enabled()) return

 do i = 1, nCellsSolve
 do k = 1, nVertLevels
 do n = 1, nTransport
    tend_scalars(transport_to_scalar(n),k,i) = tend_scalars(transport_to_scalar(n),k,i) &
                                             + (r_pbl(n,k,i)+r_cu(n,k,i))*mass(k,i)
 enddo
 enddo
 enddo

 end subroutine tracer_transport_add_tend

!=================================================================================================================
 function tracer_transport_budget(process,n) result(b)
!=================================================================================================================
!per-rank budget of gas n for TRANSPORT_PBL or TRANSPORT_CU from the last assembly.
 integer,intent(in):: process,n
 type(transport_budget_t):: b

 if(process == TRANSPORT_PBL) then
    b = bud_pbl(n)
 else
    b = bud_cu(n)
 endif

 end function tracer_transport_budget

!=================================================================================================================
 end module mpas_atmphys_tracer_transport
!=================================================================================================================
