! Unit tests for the closed-column transport kernels, budget schema and pure
! capability/cache checks (P1).
program test_tracer_column

   use mpas_kind_types, only: RKIND, StrKIND
   use unit_assert, only: check, finish
   use mpas_atmphys_tracer_column, only: column_dry_mass, column_flux_divergence, column_species_mass, &
                                         relative_mass_residual, transport_budget_t, budget_reset, &
                                         budget_accumulate, capability_failure, transport_cache_meta_t, &
                                         cache_validation, tracer_transport_check_capability, &
                                         transport_cache_validate

   implicit none

   integer, parameter :: nz = 40
   real(kind=RKIND) :: m_d3(3), F4(4), R3(3)
   real(kind=RKIND) :: m_d(nz), F(nz+1), R(nz), q0(nz), q1(nz)
   real(kind=RKIND) :: scale, resid
   type(transport_budget_t) :: b
   type(transport_cache_meta_t) :: meta
   logical :: ok
   character(len=StrKIND) :: why
   integer :: k

   ! Dry layer mass is zz*rho_zz*dz on a nonuniform column.
   m_d3 = column_dry_mass([1.1_RKIND, 0.9_RKIND, 1.0_RKIND], [1.2_RKIND, 0.8_RKIND, 0.5_RKIND], &
                          [0.0_RKIND, 100.0_RKIND, 300.0_RKIND, 600.0_RKIND])
   call check(all(abs(m_d3 - [132.0_RKIND, 144.0_RKIND, 150.0_RKIND]) <= 1.0e-12_RKIND*150.0_RKIND), &
              'dry layer mass is zz*rho_zz*dz')

   ! Closed column: any interior fluxes conserve column mass to roundoff.
   do k = 1, nz
      m_d(k) = 50.0_RKIND + 40.0_RKIND*sin(0.7_RKIND*k)**2
   end do
   F(1) = 0.0_RKIND
   F(nz+1) = 0.0_RKIND
   do k = 2, nz
      F(k) = 1.0e-3_RKIND*cos(1.3_RKIND*k)
   end do
   R = column_flux_divergence(m_d, F)
   scale = sum(abs(F))
   call check(abs(sum(m_d*R)) <= 1.0e-12_RKIND*scale, 'closed-column flux divergence conserves mass')

   ! One upward flux through the interface between layers 2 and 3.
   m_d3 = [10.0_RKIND, 20.0_RKIND, 40.0_RKIND]
   F4 = [0.0_RKIND, 0.0_RKIND, 2.0e-4_RKIND, 0.0_RKIND]
   R3 = column_flux_divergence(m_d3, F4)
   call check(all(abs(m_d3*R3 - [0.0_RKIND, -2.0e-4_RKIND, 2.0e-4_RKIND]) <= 1.0e-18_RKIND), &
              'upward flux moves mass from the lower to the upper layer')

   q0 = 1.0e-9_RKIND
   call check(abs(column_species_mass(m_d, q0) - 1.0e-9_RKIND*sum(m_d)) <= 1.0e-12_RKIND*1.0e-9_RKIND*sum(m_d), &
              'column species mass is sum(m_d*q)')

   ! Residual normalization: zero species do not divide by zero.
   q0 = 0.0_RKIND
   q1 = 0.0_RKIND
   resid = relative_mass_residual(m_d, q0, q1)
   call check(resid == 0.0_RKIND, 'zero species has zero residual')
   q1 = 1.0e-40_RKIND
   resid = relative_mass_residual(m_d, q0, q1)
   call check(abs(resid - 1.0e-10_RKIND) <= 1.0e-12_RKIND*1.0e-10_RKIND, 'zero-species residual normalizes by sum(m_d)*1e-30')
   q0 = 2.0e-9_RKIND
   q1 = q0
   q1(1) = q0(1) + 1.0e-9_RKIND
   resid = relative_mass_residual(m_d, q0, q1)
   call check(abs(resid - m_d(1)*1.0e-9_RKIND/(2.0e-9_RKIND*sum(m_d))) <= 1.0e-12_RKIND*resid, &
              'nonzero species residual is relative to its initial mass')

   ! Budget schema: per-rank accumulation over columns.
   call budget_reset(b)
   call check(b%active_columns == 0 .and. b%mass_before == 0.0_RKIND .and. b%residual == 0.0_RKIND, &
              'budget reset clears the record')
   q0 = 1.0e-9_RKIND
   q1 = q0
   q1(1) = 0.0_RKIND
   q1(2) = q0(2) + q0(1)*m_d(1)/m_d(2)
   call budget_accumulate(b, 2.0_RKIND, m_d, q0, q1, 3.0_RKIND)
   call budget_accumulate(b, 1.0_RKIND, m_d, q0, q0, 3.0_RKIND)
   call check(b%active_columns == 2, 'budget counts active columns')
   call check(abs(b%mass_before - 3.0_RKIND*1.0e-9_RKIND*sum(m_d)) <= 1.0e-12_RKIND*b%mass_before, &
              'budget mass is area-weighted')
   call check(abs(b%increment) <= 1.0e-12_RKIND*b%mass_before .and. b%residual <= 1.0e-12_RKIND, &
              'conservative redistribution has no increment or residual')
   call check(b%qmin == 0.0_RKIND .and. abs(b%qmax - q1(2)) <= 1.0e-24_RKIND .and. b%dt_applied == 3.0_RKIND, &
              'budget records extremes and applied dt')
   call check(b%boundary_flux == 0.0_RKIND, 'closed columns have no boundary flux')

   ! Capability: only an enabled, unsupported configuration fails.
   call check(capability_failure(.true., .false.), 'enabled unsupported configuration fails')
   call check(.not. capability_failure(.false., .false.), 'disabled unsupported configuration is allowed')
   call check(.not. capability_failure(.true., .true.), 'enabled supported configuration is allowed')

   ! Held-operator cache metadata.
   meta%process = 'PBL'
   meta%signature = 'bl_ysu|none|60'
   meta%refresh_time = 0.0_RKIND
   meta%valid_until = 0.0_RKIND
   call cache_validation(meta, 'bl_ysu|none|60', .false., .false., ok, why)
   call check(ok, 'disabled transport needs no cache')
   call cache_validation(meta, 'bl_ysu|none|60', .false., .true., ok, why)
   call check(.not. ok .and. len_trim(why) > 0, 'enabled transport without a restart cache fails')
   call cache_validation(meta, 'bl_mynn|none|60', .true., .true., ok, why)
   call check(.not. ok .and. index(why, 'bl_ysu|none|60') > 0 .and. index(why, 'bl_mynn|none|60') > 0, &
              'cache signature mismatch names both signatures')
   call cache_validation(meta, 'bl_ysu|none|60', .true., .true., ok, why)
   call check(ok, 'matching cache signature is accepted')

   ! The logging/fatal wrappers shared with both scheme branches exist with the
   ! sketched signatures (compile-time check; they need an initialized MPAS log).
   if (.false.) then
      call tracer_transport_check_capability('PBL', .true., 'bl_ysu', .false., 'not qualified')
      call transport_cache_validate(meta, 'bl_ysu|none|60', .false., .true.)
   end if

   call finish('test_tracer_column')

end program test_tracer_column
