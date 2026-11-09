! Test double for the PBL gas process: active in odd cells; moves half of the
! lowest layer's mass of every species into the layer above (conservative).
 module mpas_atmphys_pbl_gas
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: block_type

 implicit none
 private
 public:: pbl_gas_init,pbl_gas_enabled,pbl_gas_active,pbl_gas_apply,pbl_gas_finalize, &
         pbl_gas_begin_step,pbl_gas_refresh_time,pbl_double_active
 logical,save:: pbl_double_active = .true.  ! .false.: no column is active (an all-inactive rank)

 contains

 subroutine pbl_gas_init(block,nTransport)
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport
 end subroutine pbl_gas_init

 logical function pbl_gas_enabled()
 pbl_gas_enabled = .true.
 end function pbl_gas_enabled

 logical function pbl_gas_active(iCell)
 integer,intent(in):: iCell
 pbl_gas_active = pbl_double_active .and. mod(iCell,2) == 1
 end function pbl_gas_active

 subroutine pbl_gas_apply(iCell,dt,m_d,q)
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in):: dt,m_d(:)
 real(kind=RKIND),intent(inout):: q(:,:)
 real(kind=RKIND):: moved
 integer:: n
 do n = 1, size(q,1)
    moved  = 0.5_RKIND*q(n,1)*m_d(1)
    q(n,1) = q(n,1) - moved/m_d(1)
    q(n,2) = q(n,2) + moved/m_d(2)
 enddo
 end subroutine pbl_gas_apply

 subroutine pbl_gas_begin_step(dt)
 real(kind=RKIND),intent(in):: dt
 end subroutine pbl_gas_begin_step

 real(kind=RKIND) function pbl_gas_refresh_time()
 pbl_gas_refresh_time = -1._RKIND
 end function pbl_gas_refresh_time

 subroutine pbl_gas_finalize
 end subroutine pbl_gas_finalize

 end module mpas_atmphys_pbl_gas
