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
 use mpas_atmphys_pbl_edmf,only: pbl_edmf_updraft_fluxes,pbl_edmf_downdraft_fluxes,pbl_edmf_column
 use mpas_atmphys_tracer_column,only: tracer_transport_check_capability
 use mpas_atmphys_utilities,only: physics_error_fatal
 use bl_mynn,only: mynn_edmf_nup,mynn_edmf_ndown

 implicit none
 private
 public:: pbl_gas_init,               &
          pbl_gas_capability,         &
          pbl_gas_exchange,           &
          pbl_gas_exchange_faces,     &
          pbl_gas_cache_setup,        &
          pbl_gas_cache_store_column, &
          pbl_gas_cache_store_moist_column, &
          pbl_gas_cache_store_mynn_column,  &
          pbl_gas_cache_store_mynn_plumes,  &
          pbl_gas_mynn_updrafts,      &
          pbl_gas_mynn_downdrafts,    &
          pbl_gas_enabled,            &
          pbl_gas_active,             &
          pbl_gas_apply,              &
          pbl_gas_cache_refreshed,    &
          pbl_gas_operator_valid,     &
          pbl_gas_begin_step,         &
          pbl_gas_refresh_time,       &
          pbl_gas_finalize


!PBL vertical mixing of MICM gases (docs/plans/pbl-transport.md, P3 and P4). Interface used by
!mpas_atmphys_tracer_transport:
!   pbl_gas_init(block,nTransport)   resolve the scheme and check capability, once after chemistry_init
!   pbl_gas_enabled()                true when PBL gas transport can act in any column
!   pbl_gas_active(iCell)            true when the cached operator acts in column iCell this step
!   pbl_gas_begin_step(dt)           once per assembly, before any pbl_gas_apply: stops the model unless the held
!                                    operator was refreshed within its PBL interval
!   pbl_gas_apply(iCell,dt,m_d,q)    q(nTransport,nVertLevels) -> q after dt [kg/kg dry air], closed column
!   pbl_gas_finalize()
!
!Each time the PBL scheme runs (every step at the default interval, else when the pbl alarm rings), the PBL
!driver stores the scheme's scalar exchange as interface exchange G [kg m-2 s-1] (pbl_gas_cache_store_column)
!before chemistry runs, then stamps the cache (pbl_gas_cache_refreshed). Every model step, transport assembly
!applies the held operator, closed-column backward-Euler diffusion with the model dt, to the post-chemistry
!gases, as the held meteorological PBL tendencies are applied between PBL calls. Only the operator is held,
!never a gas increment; it is not restarted, because the pbl alarm rings at a restart's first step and
!rebuilds it before use. Supported with a surface layer scheme: YSU (with
!or without top-down mixing) and MYNN, local diffusion with or without its mass-flux scheme (config_mynn_edmf,
!config_mynn_edmf_dd); PBL off is a supported zero-transport case. With MYNN's mass-flux scheme the driver also
!stores each updraft's and downdraft's dry mass flux and mass-consistent entrainment
!(pbl_gas_cache_store_mynn_plumes), and pbl_gas_apply follows the diffusion with the plume transport of
!mpas_atmphys_pbl_edmf, rebuilt from the post-chemistry gases. Any config_pbl_interval is supported (P5).


 logical,save:: enabled = .false.
 logical,save:: mynn_updrafts = .false.
 logical,save:: mynn_downdrafts = .false.
 real(kind=RKIND),dimension(:,:),allocatable,save:: g_cache   ! (nVertLevels+1,nCellsSolve)
 logical,dimension(:),allocatable,save:: mask_cache            ! (nCellsSolve)
 real(kind=RKIND),dimension(:,:,:),allocatable,save:: mu_cache ! (nVertLevels+1,nup,nCellsSolve) plume mass flux
 real(kind=RKIND),dimension(:,:,:),allocatable,save:: eu_cache ! (nVertLevels,nup,nCellsSolve) plume entrainment
 real(kind=RKIND),dimension(:,:,:),allocatable,save:: nd_cache ! (nVertLevels+1,ndown,nCellsSolve) downdraft flux
 real(kind=RKIND),dimension(:,:,:),allocatable,save:: ed_cache ! (nVertLevels,ndown,nCellsSolve) its entrainment

!refresh stamp of the held operator: model steps since the PBL scheme last refreshed it (-1 before the first
!refresh), and the model time [s] and PBL interval dt_pbl [s] of that refresh.
 integer,save:: refresh_age = -1
 real(kind=RKIND),save:: refresh_time = -1._RKIND
 real(kind=RKIND),save:: refresh_interval = 0._RKIND


 contains


!=================================================================================================================
 subroutine pbl_gas_init(block,nTransport)
!=================================================================================================================
!block gives configs, structs and dimensions; nTransport is the number of independently transported gases.
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport

 logical,pointer:: requested
 character(len=StrKIND),pointer:: pbl_scheme,sfclayer_scheme,pbl_interval
 integer,pointer:: mynn_edmf,mynn_edmf_dd
 integer,pointer:: nVertLevels,nCellsSolve
 integer:: edmf,edmf_dd
 logical:: supported
 character(len=StrKIND):: reason

 call pbl_gas_finalize()

 nullify(requested,pbl_scheme,sfclayer_scheme,pbl_interval,mynn_edmf,mynn_edmf_dd)
 call mpas_pool_get_config(block%configs,'config_micm_pbl_transport',requested)
 call mpas_pool_get_config(block%configs,'config_pbl_scheme'        ,pbl_scheme)
 call mpas_pool_get_config(block%configs,'config_sfclayer_scheme'   ,sfclayer_scheme)
 call mpas_pool_get_config(block%configs,'config_pbl_interval'      ,pbl_interval)
 call mpas_pool_get_config(block%configs,'config_mynn_edmf'         ,mynn_edmf)
 call mpas_pool_get_config(block%configs,'config_mynn_edmf_dd'      ,mynn_edmf_dd)
 if(.not.associated(requested)) return
 if(.not.requested) return

 if(nTransport == 0) then
    call mpas_log_write('[transport] PBL: config_micm_pbl_transport is set but no MICM gas is transported; no-op')
    return
 endif

 edmf    = 0
 edmf_dd = 0
 if(associated(mynn_edmf)   ) edmf    = mynn_edmf
 if(associated(mynn_edmf_dd)) edmf_dd = mynn_edmf_dd
 call pbl_gas_capability(pbl_scheme,sfclayer_scheme,pbl_interval,supported,reason)
 call tracer_transport_check_capability('PBL',.true.,pbl_scheme,supported,reason)
 enabled = trim(pbl_scheme) /= 'off'
 if(enabled) call mpas_log_write('[transport] PBL: gas operator refreshed at each PBL call (config_pbl_interval = '// &
                                 trim(pbl_interval)//') and applied every model step')
!MYNN runs its downdrafts only inside its updraft scheme.
 mynn_updrafts   = trim(pbl_scheme) == 'bl_mynn' .and. edmf == 1
 mynn_downdrafts = mynn_updrafts .and. edmf_dd == 1
 if(trim(pbl_scheme) == 'bl_mynn') then
    if(mynn_downdrafts) then
       call mpas_log_write('[transport] PBL: MYNN gas operator = local diffusion (exch_h) + updraft EDMF ($i plumes)'// &
                           ' + downdraft EDMF ($i downdrafts)',intArgs=(/mynn_edmf_nup,mynn_edmf_ndown/))
    elseif(mynn_updrafts) then
       call mpas_log_write('[transport] PBL: MYNN gas operator = local diffusion (exch_h) + updraft EDMF ($i plumes)', &
                           intArgs=(/mynn_edmf_nup/))
    else
       call mpas_log_write('[transport] PBL: MYNN gas operator = local diffusion (exch_h)')
    endif
 endif
 if(enabled) then
    call mpas_pool_get_dimension(block%dimensions,'nVertLevels',nVertLevels)
    call mpas_pool_get_dimension(block%dimensions,'nCellsSolve',nCellsSolve)
    call pbl_gas_cache_setup(nVertLevels,nCellsSolve,nup=merge(mynn_edmf_nup,0,mynn_updrafts), &
                             ndown=merge(mynn_edmf_ndown,0,mynn_downdrafts))
 endif

 end subroutine pbl_gas_init

!=================================================================================================================
 pure subroutine pbl_gas_capability(pbl_scheme,sfclayer_scheme,pbl_interval,supported,reason)
!=================================================================================================================
!qualified PBL configurations for MICM gas transport: YSU (P3) and MYNN in every mass-flux mode (P4), at any
!PBL interval (P5: the held operator is applied every step, pbl_gas_begin_step). pbl_interval is accepted
!so interval-specific restrictions can be stated here; none is needed today.
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
       supported = .true.
    case default
       reason = 'no MICM gas adapter for this PBL scheme'
 end select

 if(supported .and. trim(pbl_scheme) /= 'off') then
    if(trim(sfclayer_scheme) == 'off') then
       supported = .false.
       reason    = 'an active PBL scheme requires a surface layer scheme'
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
 pure subroutine pbl_gas_exchange_faces(kh,rho_d,dz,G,ok)
!=================================================================================================================
!interface exchange [kg m-2 s-1] for the closed-column diffusion kernel from an exchange coefficient on the
!bottom face of each layer, kh(k) [m2 s-1] between layers k-1 and k (MYNN's exch_h; kh(1), the surface, is
!ignored):
!   G(k) = rhoz_d(k)*kh(k)/dzk,  dzk = 0.5*(dz(k)+dz(k-1)),
!   rhoz_d(k) = (rho_d(k)*dz(k-1)+rho_d(k-1)*dz(k))/(dz(k-1)+dz(k)),   k = 2..nz,
!MYNN's khdz = rhoz*dfh with the dry density (exch_h = dfh*dzk). G(1) = G(nz+1) = 0. ok is false, and G is
!zero, unless kh(2:nz) >= 0, rho_d > 0 and dz > 0 are finite.
 use,intrinsic:: ieee_arithmetic,only: ieee_is_finite
 real(kind=RKIND),intent(in),dimension(:):: kh      ! (nz) [m2 s-1]
 real(kind=RKIND),intent(in),dimension(:):: rho_d   ! (nz) dry air density [kg m-3]
 real(kind=RKIND),intent(in),dimension(:):: dz      ! (nz) layer thickness [m]
 real(kind=RKIND),intent(out),dimension(:):: G      ! (nz+1)
 logical,intent(out):: ok
 real(kind=RKIND):: rhoz,dzk
 integer:: k,nz

 nz = size(rho_d)
 G(:) = 0._RKIND
 ok = all(ieee_is_finite(rho_d)) .and. all(ieee_is_finite(dz)) .and. all(ieee_is_finite(kh(2:nz)))
 if(ok) ok = all(rho_d > 0._RKIND) .and. all(dz > 0._RKIND) .and. all(kh(2:nz) >= 0._RKIND)
 if(.not.ok) return

 do k = 2, nz
    dzk  = 0.5_RKIND*(dz(k)+dz(k-1))
    rhoz = (rho_d(k)*dz(k-1)+rho_d(k-1)*dz(k))/(dz(k-1)+dz(k))
    G(k) = rhoz*kh(k)/dzk
 enddo

 end subroutine pbl_gas_exchange_faces

!=================================================================================================================
 subroutine pbl_gas_cache_setup(nlevels,ncells,nup,ndown)
!=================================================================================================================
!allocate the per-column operator cache, with nup updraft plumes and ndown downdrafts per column when present and
!positive; every column starts inactive.
 integer,intent(in):: nlevels,ncells
 integer,intent(in),optional:: nup,ndown

 if(allocated(g_cache)   ) deallocate(g_cache)
 if(allocated(mask_cache)) deallocate(mask_cache)
 if(allocated(mu_cache)  ) deallocate(mu_cache)
 if(allocated(eu_cache)  ) deallocate(eu_cache)
 if(allocated(nd_cache)  ) deallocate(nd_cache)
 if(allocated(ed_cache)  ) deallocate(ed_cache)
 allocate(g_cache(nlevels+1,ncells),mask_cache(ncells))
 g_cache(:,:)  = 0._RKIND
 mask_cache(:) = .false.
 if(present(nup)) then
    if(nup > 0) then
       allocate(mu_cache(nlevels+1,nup,ncells),eu_cache(nlevels,nup,ncells))
       mu_cache(:,:,:) = 0._RKIND
       eu_cache(:,:,:) = 0._RKIND
    endif
 endif
 if(present(ndown)) then
    if(ndown > 0) then
       allocate(nd_cache(nlevels+1,ndown,ncells),ed_cache(nlevels,ndown,ncells))
       nd_cache(:,:,:) = 0._RKIND
       ed_cache(:,:,:) = 0._RKIND
    endif
 endif

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
 subroutine pbl_gas_cache_store_mynn_column(iCell,kh,rho,qv,dz,ok)
!=================================================================================================================
!store this step's MYNN operator for column iCell from exch_h on the bottom face of each layer and the PBL
!driver's moist density rho = rho_dry*(1+qv) (pbl_gas_exchange_faces with rho_dry = rho/(1+qv)). An invalid
!column (ok = .false.) is left inactive; the caller decides whether to stop.
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in),dimension(:):: kh,rho,qv,dz
 logical,intent(out):: ok

 call pbl_gas_exchange_faces(kh,rho/(1._RKIND+qv),dz,g_cache(:,iCell),ok)
 mask_cache(iCell) = ok .and. any(g_cache(:,iCell) > 0._RKIND)

 end subroutine pbl_gas_cache_store_mynn_column

!=================================================================================================================
 subroutine pbl_gas_cache_store_mynn_plumes(iCell,aw_up,ent_up,rho,qv,dz,ok,aw_dn,ent_dn)
!=================================================================================================================
!store this step's MYNN updraft plumes for column iCell, after pbl_gas_cache_store_mynn_column: per plume, the
!exported area-weighted velocity aw_up(k,:) on the top face of layer k and mixing fraction ent_up(k,:), with the
!dry density rho/(1+qv) (pbl_edmf_updraft_fluxes); and, when the cache holds downdrafts, aw_dn(k,:) (downward,
!bottom face of layer k) and ent_dn(k,:) (pbl_edmf_downdraft_fluxes), which must then be present. An invalid
!column (ok = .false.) is left inactive; the caller decides whether to stop.
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in),dimension(:,:):: aw_up,ent_up   ! (nVertLevels,nup)
 real(kind=RKIND),intent(in),dimension(:):: rho,qv,dz
 logical,intent(out):: ok
 real(kind=RKIND),intent(in),dimension(:,:),optional:: aw_dn,ent_dn   ! (nVertLevels,ndown)
 logical:: ok_dn

 call pbl_edmf_updraft_fluxes(aw_up,ent_up,rho/(1._RKIND+qv),dz,mu_cache(:,:,iCell),eu_cache(:,:,iCell),ok)
 if(allocated(nd_cache)) then
    ok_dn = present(aw_dn) .and. present(ent_dn)
    if(ok_dn) call pbl_edmf_downdraft_fluxes(aw_dn,ent_dn,rho/(1._RKIND+qv),dz,nd_cache(:,:,iCell), &
                                             ed_cache(:,:,iCell),ok_dn)
    ok = ok .and. ok_dn
 endif
 if(ok) then
    mask_cache(iCell) = mask_cache(iCell) .or. any(mu_cache(:,:,iCell) > 0._RKIND)
    if(allocated(nd_cache)) mask_cache(iCell) = mask_cache(iCell) .or. any(nd_cache(:,:,iCell) > 0._RKIND)
 else
    mask_cache(iCell) = .false.
 endif

 end subroutine pbl_gas_cache_store_mynn_plumes

!=================================================================================================================
 logical function pbl_gas_mynn_updrafts()
!=================================================================================================================
!true when gas transport uses MYNN's updraft plumes, which the PBL driver must then export.

 pbl_gas_mynn_updrafts = enabled .and. mynn_updrafts

 end function pbl_gas_mynn_updrafts

!=================================================================================================================
 logical function pbl_gas_mynn_downdrafts()
!=================================================================================================================
!true when gas transport also uses MYNN's downdrafts, which the PBL driver must then export.

 pbl_gas_mynn_downdrafts = enabled .and. mynn_downdrafts

 end function pbl_gas_mynn_downdrafts

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
 real(kind=RKIND),dimension(size(m_d)+1,0):: no_flux     ! no downdrafts
 real(kind=RKIND),dimension(size(m_d),0):: no_entrainment
 logical:: ok,plumes
 integer:: nsub

 call pbl_diffuse_column(dt,m_d,g_cache(:,iCell),q)
 if(.not.allocated(mu_cache)) return
 plumes = any(mu_cache(:,:,iCell) > 0._RKIND)
 if(allocated(nd_cache)) plumes = plumes .or. any(nd_cache(:,:,iCell) > 0._RKIND)
 if(.not.plumes) return
 if(allocated(nd_cache)) then
    call pbl_edmf_column(dt,m_d,mu_cache(:,:,iCell),eu_cache(:,:,iCell),nd_cache(:,:,iCell),ed_cache(:,:,iCell), &
                         q,nsub,ok)
 else
    call pbl_edmf_column(dt,m_d,mu_cache(:,:,iCell),eu_cache(:,:,iCell),no_flux,no_entrainment,q,nsub,ok)
 endif
 if(.not.ok) call physics_error_fatal('[transport] PBL: MYNN mass-flux operator for MICM gases needs more '// &
                                      'than 1e4 substeps (non-finite or implausible plume mass flux)')

 end subroutine pbl_gas_apply

!=================================================================================================================
 subroutine pbl_gas_cache_refreshed(time_s,interval_s)
!=================================================================================================================
!stamp the cache after the PBL driver has stored every owned column in this PBL call: time_s is the model time
!of the call [s], interval_s the PBL interval dt_pbl [s] (config_dt at the default interval).
 real(kind=RKIND),intent(in):: time_s,interval_s

 refresh_age      = 0
 refresh_time     = time_s
 refresh_interval = interval_s

 end subroutine pbl_gas_cache_refreshed

!=================================================================================================================
 pure logical function pbl_gas_operator_valid(age,dt,interval)
!=================================================================================================================
!true when an operator refreshed age model steps ago (age 0: this step) may still be applied with time step dt:
!it covers [refresh, refresh + interval), so age*dt < interval. The pbl alarm rings at the first step at or
!after each interval, so the largest age reached is ceiling(interval/dt)-1. The relative margin keeps a step
!that lands on the next refresh time stale despite rounding.
 integer,intent(in):: age
 real(kind=RKIND),intent(in):: dt,interval

 pbl_gas_operator_valid = age >= 0
 if(pbl_gas_operator_valid) &
    pbl_gas_operator_valid = real(age,kind=RKIND)*dt < interval*(1._RKIND-1.e-10_RKIND)

 end function pbl_gas_operator_valid

!=================================================================================================================
 subroutine pbl_gas_begin_step(dt)
!=================================================================================================================
!once per transport assembly, before pbl_gas_apply: stop the model if the held operator has never been built
!or has outlived its PBL interval (a refresh was missed), then age it by one model step.
 real(kind=RKIND),intent(in):: dt

 if(refresh_age < 0) &
    call physics_error_fatal('[transport] PBL: the MICM gas operator is used before the PBL scheme has run')
 if(.not.pbl_gas_operator_valid(refresh_age,dt,refresh_interval)) &
    call physics_error_fatal('[transport] PBL: the held MICM gas operator is older than config_pbl_interval '// &
                             '(the PBL scheme did not refresh it)')
 refresh_age = refresh_age + 1

 end subroutine pbl_gas_begin_step

!=================================================================================================================
 real(kind=RKIND) function pbl_gas_refresh_time()
!=================================================================================================================
!model time [s] of the last refresh; -1 before the first.

 pbl_gas_refresh_time = refresh_time

 end function pbl_gas_refresh_time

!=================================================================================================================
 subroutine pbl_gas_finalize
!=================================================================================================================

 if(allocated(g_cache)   ) deallocate(g_cache)
 if(allocated(mask_cache)) deallocate(mask_cache)
 if(allocated(mu_cache)  ) deallocate(mu_cache)
 if(allocated(eu_cache)  ) deallocate(eu_cache)
 if(allocated(nd_cache)  ) deallocate(nd_cache)
 if(allocated(ed_cache)  ) deallocate(ed_cache)
 enabled          = .false.
 mynn_updrafts    = .false.
 mynn_downdrafts  = .false.
 refresh_age      = -1
 refresh_time     = -1._RKIND
 refresh_interval = 0._RKIND

 end subroutine pbl_gas_finalize

!=================================================================================================================
 end module mpas_atmphys_pbl_gas
!=================================================================================================================
