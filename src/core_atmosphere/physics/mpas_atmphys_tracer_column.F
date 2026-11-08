! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_tracer_column
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_log,only: mpas_log_write
 use mpas_atmphys_utilities,only: physics_error_fatal

 implicit none
 private
 public:: column_dry_mass,         &
          column_flux_divergence,  &
          column_species_mass,     &
          relative_mass_residual,  &
          transport_budget_t,      &
          budget_reset,            &
          budget_accumulate,       &
          capability_failure,      &
          transport_cache_meta_t,  &
          cache_validation,        &
          tracer_transport_check_capability, &
          transport_cache_validate


!Closed-column kernels shared by the MICM gas PBL and convective transport (P1 of docs/plans/pbl-transport.md).
!MPAS level k = 1 is at the ground; F(k) is the upward-positive tracer flux through the lower interface of
!layer k, so F(1) is the ground and F(nVertLevels+1) the model top. Units:
!   q   : kg species / kg dry air
!   R   : dq/dt [s-1]
!   m_d : dry layer mass per area = zz*rho_zz*dz [kg m-2], not the moist physics density rho_p*dz
!   F   : [kg species m-2 s-1], with m_d(k)*R(k) = F(k) - F(k+1)
!The kernels and checks are pure and can be unit tested without an initialized MPAS log. Only the two
!wrappers at the end, tracer_transport_check_capability and transport_cache_validate, log and stop the model;
!the scheme branches call them from their *_gas_init routines.


!schema of the per-species transport diagnostics, accumulated over owned columns on one rank:
 type transport_budget_t
    real(kind=RKIND):: mass_before   = 0._RKIND   ! sum of area*m_d*q before the process [kg]
    real(kind=RKIND):: mass_after    = 0._RKIND   ! sum of area*m_d*q after the process [kg]
    real(kind=RKIND):: increment     = 0._RKIND   ! mass_after - mass_before [kg]
    real(kind=RKIND):: boundary_flux = 0._RKIND   ! net inflow through column boundaries [kg s-1]; 0 when closed
    real(kind=RKIND):: residual      = 0._RKIND   ! largest per-column relative mass residual
    real(kind=RKIND):: qmin          = huge(1._RKIND)
    real(kind=RKIND):: qmax          = -huge(1._RKIND)
    integer:: active_columns         = 0
    real(kind=RKIND):: refresh_time  = -1._RKIND  ! model time of the operator refresh [s]; -1 if unknown
    real(kind=RKIND):: dt_applied    = 0._RKIND   ! time step the tendency was built with [s]
 end type transport_budget_t

!metadata of a held meteorological transport operator that must survive restart:
 type transport_cache_meta_t
    character(len=StrKIND):: process   = ''       ! 'PBL' or 'convection'
    character(len=StrKIND):: signature = ''       ! scheme, closure, interval and vertical-grid identity
    real(kind=RKIND):: refresh_time    = -1._RKIND ! [s] last meteorological refresh
    real(kind=RKIND):: valid_until     = -1._RKIND ! [s] next refresh or expiry
 end type transport_cache_meta_t


 contains


!=================================================================================================================
 pure function column_dry_mass(zz,rho_zz,zgrid) result(m_d)
!=================================================================================================================
!dry layer mass per area [kg m-2]. zgrid holds the size(zz)+1 interface heights.
 real(kind=RKIND),intent(in):: zz(:),rho_zz(:),zgrid(:)
 real(kind=RKIND):: m_d(size(zz))
 integer:: k

 do k = 1, size(zz)
    m_d(k) = zz(k)*rho_zz(k)*(zgrid(k+1)-zgrid(k))
 enddo

 end function column_dry_mass

!=================================================================================================================
 pure function column_flux_divergence(m_d,F) result(R)
!=================================================================================================================
!tendency [s-1] from interface fluxes F(1:size(m_d)+1); F(1) = F(nz+1) = 0 closes the column.
 real(kind=RKIND),intent(in):: m_d(:),F(:)
 real(kind=RKIND):: R(size(m_d))
 integer:: k

 do k = 1, size(m_d)
    R(k) = (F(k)-F(k+1))/m_d(k)
 enddo

 end function column_flux_divergence

!=================================================================================================================
 pure real(kind=RKIND) function column_species_mass(m_d,q)
!=================================================================================================================
!species mass per area [kg m-2].
 real(kind=RKIND),intent(in):: m_d(:),q(:)

 column_species_mass = sum(m_d*q)

 end function column_species_mass

!=================================================================================================================
 pure real(kind=RKIND) function relative_mass_residual(m_d,q_before,q_after)
!=================================================================================================================
!|change of column species mass| relative to max(initial mass, sum(m_d)*1e-30); a zero species stays finite.
 real(kind=RKIND),intent(in):: m_d(:),q_before(:),q_after(:)
 real(kind=RKIND):: mass_before,mass_after,norm

 mass_before = column_species_mass(m_d,q_before)
 mass_after  = column_species_mass(m_d,q_after)
 norm = max(mass_before,sum(m_d)*1.e-30_RKIND)
 relative_mass_residual = abs(mass_after-mass_before)/norm

 end function relative_mass_residual

!=================================================================================================================
 pure subroutine budget_reset(b)
!=================================================================================================================
 type(transport_budget_t),intent(out):: b

!intent(out) applies the default initialization of every component.
 b%active_columns = 0

 end subroutine budget_reset

!=================================================================================================================
 pure subroutine budget_accumulate(b,area,m_d,q_before,q_after,dt)
!=================================================================================================================
!add one closed column (area [m2]) to the budget of one species and one process.
 type(transport_budget_t),intent(inout):: b
 real(kind=RKIND),intent(in):: area,dt
 real(kind=RKIND),intent(in):: m_d(:),q_before(:),q_after(:)
 real(kind=RKIND):: mass_before,mass_after

 mass_before = area*column_species_mass(m_d,q_before)
 mass_after  = area*column_species_mass(m_d,q_after)

 b%mass_before    = b%mass_before + mass_before
 b%mass_after     = b%mass_after  + mass_after
 b%increment      = b%increment   + (mass_after-mass_before)
 b%residual       = max(b%residual,relative_mass_residual(m_d,q_before,q_after))
 b%qmin           = min(b%qmin,minval(q_after))
 b%qmax           = max(b%qmax,maxval(q_after))
 b%active_columns = b%active_columns + 1
 b%dt_applied     = dt

 end subroutine budget_accumulate

!=================================================================================================================
 pure logical function capability_failure(enabled,supported)
!=================================================================================================================
!true when gas transport is requested for a configuration that has no qualified adapter.
 logical,intent(in):: enabled,supported

 capability_failure = enabled .and. .not.supported

 end function capability_failure

!=================================================================================================================
 pure subroutine cache_validation(meta,expected_signature,restart_has_cache,enabled,ok,why)
!=================================================================================================================
!held operators are restarted; per-step gas tendencies are reconstructed, never replayed.
 type(transport_cache_meta_t),intent(in):: meta
 character(len=*),intent(in):: expected_signature
 logical,intent(in):: restart_has_cache,enabled
 logical,intent(out):: ok
 character(len=*),intent(out):: why

 ok  = .true.
 why = ''
 if(.not.enabled) return

 if(.not.restart_has_cache) then
    ok  = .false.
    why = trim(meta%process)//' gas transport is enabled but the restart file has no held-operator cache'
 elseif(trim(meta%signature) /= trim(expected_signature)) then
    ok  = .false.
    why = trim(meta%process)//' restart cache signature mismatch: '//trim(meta%signature)// &
          ' /= '//trim(expected_signature)
 endif

 end subroutine cache_validation

!=================================================================================================================
 subroutine tracer_transport_check_capability(process,enabled,scheme,supported,reason)
!=================================================================================================================
!log the resolved scheme and capability; stop with an actionable message when gas transport is requested for an
!unqualified configuration. Scheme-specific support tables live in mpas_atmphys_pbl_gas and mpas_atmphys_cu_gas.
 character(len=*),intent(in):: process,scheme,reason
 logical,intent(in):: enabled,supported

 call mpas_log_write('[transport] '//trim(process)//': resolved scheme = '//trim(scheme)// &
                     ', enabled = $l, supported = $l',logicArgs=(/enabled,supported/))
 if(capability_failure(enabled,supported)) &
    call physics_error_fatal('[transport] MICM gas '//trim(process)//' transport is not supported with '// &
                             trim(scheme)//': '//trim(reason))

 end subroutine tracer_transport_check_capability

!=================================================================================================================
 subroutine transport_cache_validate(meta,expected_signature,restart_has_cache,enabled)
!=================================================================================================================
!stop when an enabled process restarts without a compatible held-operator cache.
 type(transport_cache_meta_t),intent(in):: meta
 character(len=*),intent(in):: expected_signature
 logical,intent(in):: restart_has_cache,enabled
 logical:: ok
 character(len=StrKIND):: why

 call cache_validation(meta,expected_signature,restart_has_cache,enabled,ok,why)
 if(.not.ok) call physics_error_fatal('[transport] '//trim(why))

 end subroutine transport_cache_validate

!=================================================================================================================
 end module mpas_atmphys_tracer_column
!=================================================================================================================
