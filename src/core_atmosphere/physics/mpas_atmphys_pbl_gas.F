! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_pbl_gas
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_derived_types,only: block_type
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_dimension
 use mpas_log,only: mpas_log_write
 use mpas_atmphys_pbl_diffusion,only: pbl_diffuse_column
 use mpas_atmphys_tracer_column,only: tracer_transport_check_capability

 implicit none
 private
 public:: pbl_gas_init,               &
          pbl_gas_capability,         &
          pbl_gas_exchange,           &
          pbl_gas_cache_setup,        &
          pbl_gas_cache_store_column, &
          pbl_gas_cache_store_moist_column, &
          pbl_gas_enabled,            &
          pbl_gas_active,             &
          pbl_gas_apply,              &
          pbl_gas_finalize


!PBL vertical mixing of MICM gases (docs/plans/pbl-transport.md, P3). Interface used by
!mpas_atmphys_tracer_transport:
!   pbl_gas_init(block,nTransport)   resolve the scheme and check capability, once after chemistry_init
!   pbl_gas_enabled()                true when PBL gas transport can act in any column
!   pbl_gas_active(iCell)            true when the cached operator acts in column iCell this step
!   pbl_gas_apply(iCell,dt,m_d,q)    q(nTransport,nVertLevels) -> q after dt [kg/kg dry air], closed column
!   pbl_gas_finalize()
!
!Each model step, the PBL driver stores the scheme's scalar exchange as interface exchange G [kg m-2 s-1]
!(pbl_gas_cache_store_column) before chemistry runs; transport assembly in the same step applies the
!closed-column backward-Euler diffusion with the model dt to the post-chemistry gases. Only the operator is
!held, never a gas increment. Supported: YSU (with or without top-down mixing) at the default PBL interval,
!with a surface layer scheme; PBL off is a supported zero-transport case. MYNN (P4) and nondefault PBL
!intervals (P5) stop the model when gas transport is requested.


 logical,save:: enabled = .false.
 real(kind=RKIND),dimension(:,:),allocatable,save:: g_cache   ! (nVertLevels+1,nCellsSolve)
 logical,dimension(:),allocatable,save:: mask_cache            ! (nCellsSolve)


 contains


!=================================================================================================================
 subroutine pbl_gas_init(block,nTransport)
!=================================================================================================================
!block gives configs, structs and dimensions; nTransport is the number of independently transported gases.
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport

 logical,pointer:: requested
 character(len=StrKIND),pointer:: pbl_scheme,sfclayer_scheme,pbl_interval
 integer,pointer:: nVertLevels,nCellsSolve
 logical:: supported
 character(len=StrKIND):: reason

 call pbl_gas_finalize()

 nullify(requested,pbl_scheme,sfclayer_scheme,pbl_interval)
 call mpas_pool_get_config(block%configs,'config_micm_pbl_transport',requested)
 call mpas_pool_get_config(block%configs,'config_pbl_scheme'        ,pbl_scheme)
 call mpas_pool_get_config(block%configs,'config_sfclayer_scheme'   ,sfclayer_scheme)
 call mpas_pool_get_config(block%configs,'config_pbl_interval'      ,pbl_interval)
 if(.not.associated(requested)) return
 if(.not.requested) return

 if(nTransport == 0) then
    call mpas_log_write('[transport] PBL: config_micm_pbl_transport is set but no MICM gas is transported; no-op')
    return
 endif

 call pbl_gas_capability(pbl_scheme,sfclayer_scheme,pbl_interval,supported,reason)
 call tracer_transport_check_capability('PBL',.true.,pbl_scheme,supported,reason)

 enabled = trim(pbl_scheme) /= 'off'
 if(enabled) then
    call mpas_pool_get_dimension(block%dimensions,'nVertLevels',nVertLevels)
    call mpas_pool_get_dimension(block%dimensions,'nCellsSolve',nCellsSolve)
    call pbl_gas_cache_setup(nVertLevels,nCellsSolve)
 endif

 end subroutine pbl_gas_init

!=================================================================================================================
 pure subroutine pbl_gas_capability(pbl_scheme,sfclayer_scheme,pbl_interval,supported,reason)
!=================================================================================================================
!qualified PBL configurations for MICM gas transport.
 character(len=*),intent(in):: pbl_scheme,sfclayer_scheme,pbl_interval
 logical,intent(out):: supported
 character(len=*),intent(out):: reason

 supported = .false.
 reason    = ''
 select case(trim(pbl_scheme))
    case('off')
       supported = .true.
    case('bl_ysu')
       supported = .true.
    case('bl_mynn')
       reason = 'MYNN gas transport is not yet qualified (P4 of docs/plans/pbl-transport.md)'
    case default
       reason = 'no MICM gas adapter for this PBL scheme'
 end select

 if(supported .and. trim(pbl_scheme) /= 'off') then
    if(trim(sfclayer_scheme) == 'off') then
       supported = .false.
       reason    = 'an active PBL scheme requires a surface layer scheme'
    elseif(trim(pbl_interval) /= 'none') then
       supported = .false.
       reason    = 'config_pbl_interval must be ''none'' until held PBL operators are qualified (P5)'
    endif
 endif

 end subroutine pbl_gas_capability

!=================================================================================================================
 pure subroutine pbl_gas_exchange(kzq,m_d,dz,G,ok)
!=================================================================================================================
!interface exchange [kg m-2 s-1] for the closed-column diffusion kernel from a PBL scheme's scalar exchange
!coefficient kzq(k) [m2 s-1] at the interface between layers k and k+1 (kzq(nz) is ignored):
!   G(k+1) = 0.5*(m_d(k)+m_d(k+1))*kzq(k)/dzc**2 = rho_dry_interface*kzq(k)/dzc,  dzc = 0.5*(dz(k)+dz(k+1)),
!the dry-air form of YSU's dsig/(g*dza)*kzq/dza. G(1) = G(nz+1) = 0. ok is false, and G is zero, unless
!kzq >= 0, m_d > 0 and dz > 0 are finite at every level used.
 use,intrinsic:: ieee_arithmetic,only: ieee_is_finite
 real(kind=RKIND),intent(in),dimension(:):: kzq   ! (nz) [m2 s-1]
 real(kind=RKIND),intent(in),dimension(:):: m_d   ! (nz) dry layer mass [kg m-2]
 real(kind=RKIND),intent(in),dimension(:):: dz    ! (nz) layer thickness [m]
 real(kind=RKIND),intent(out),dimension(:):: G    ! (nz+1)
 logical,intent(out):: ok
 real(kind=RKIND):: dzc
 integer:: k,nz

 nz = size(m_d)
 G(:) = 0._RKIND
 ok = all(ieee_is_finite(m_d)) .and. all(ieee_is_finite(dz)) .and. all(ieee_is_finite(kzq(1:nz-1)))
 if(ok) ok = all(m_d > 0._RKIND) .and. all(dz > 0._RKIND) .and. all(kzq(1:nz-1) >= 0._RKIND)
 if(.not.ok) return

 do k = 1, nz-1
    dzc = 0.5_RKIND*(dz(k)+dz(k+1))
    G(k+1) = 0.5_RKIND*(m_d(k)+m_d(k+1))*kzq(k)/(dzc*dzc)
 enddo

 end subroutine pbl_gas_exchange

!=================================================================================================================
 subroutine pbl_gas_cache_setup(nlevels,ncells)
!=================================================================================================================
!allocate the per-column operator cache; every column starts inactive.
 integer,intent(in):: nlevels,ncells

 if(allocated(g_cache)   ) deallocate(g_cache)
 if(allocated(mask_cache)) deallocate(mask_cache)
 allocate(g_cache(nlevels+1,ncells),mask_cache(ncells))
 g_cache(:,:)  = 0._RKIND
 mask_cache(:) = .false.

 end subroutine pbl_gas_cache_setup

!=================================================================================================================
 subroutine pbl_gas_cache_store_column(iCell,kzq,m_d,dz,ok)
!=================================================================================================================
!store this step's operator for column iCell from the scheme's scalar exchange coefficient. An invalid
!column (ok = .false.) is left inactive; the caller decides whether to stop.
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in),dimension(:):: kzq,m_d,dz
 logical,intent(out):: ok

 call pbl_gas_exchange(kzq,m_d,dz,g_cache(:,iCell),ok)
 mask_cache(iCell) = ok .and. any(g_cache(:,iCell) > 0._RKIND)

 end subroutine pbl_gas_cache_store_column

!=================================================================================================================
 subroutine pbl_gas_cache_store_moist_column(iCell,kzq,rho,qv,dz,ok)
!=================================================================================================================
!as pbl_gas_cache_store_column, from the PBL driver's moist density rho = rho_dry*(1+qv) [kg m-3] (MPAS rho_p,
!with qv = max(0,qv)): the dry layer mass is rho/(1+qv)*dz = zz*rho_zz*dz.
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in),dimension(:):: kzq,rho,qv,dz
 logical,intent(out):: ok

 call pbl_gas_cache_store_column(iCell,kzq,rho/(1._RKIND+qv)*dz,dz,ok)

 end subroutine pbl_gas_cache_store_moist_column

!=================================================================================================================
 logical function pbl_gas_enabled()
!=================================================================================================================

 pbl_gas_enabled = enabled

 end function pbl_gas_enabled

!=================================================================================================================
 logical function pbl_gas_active(iCell)
!=================================================================================================================
 integer,intent(in):: iCell

 pbl_gas_active = .false.
 if(allocated(mask_cache)) pbl_gas_active = mask_cache(iCell)

 end function pbl_gas_active

!=================================================================================================================
 subroutine pbl_gas_apply(iCell,dt,m_d,q)
!=================================================================================================================
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in):: dt                       ! model time step [s] (not YSU's 2*dt)
 real(kind=RKIND),intent(in),dimension(:):: m_d         ! dry layer mass [kg m-2]
 real(kind=RKIND),intent(inout),dimension(:,:):: q      ! (nTransport,nVertLevels), post-chemistry

 call pbl_diffuse_column(dt,m_d,g_cache(:,iCell),q)

 end subroutine pbl_gas_apply

!=================================================================================================================
 subroutine pbl_gas_finalize
!=================================================================================================================

 if(allocated(g_cache)   ) deallocate(g_cache)
 if(allocated(mask_cache)) deallocate(mask_cache)
 enabled  = .false.

 end subroutine pbl_gas_finalize

!=================================================================================================================
 end module mpas_atmphys_pbl_gas
!=================================================================================================================
