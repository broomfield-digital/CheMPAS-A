! Unit tests for the PBL gas process (P3): the capability gate and the per-column operator cache
! that the PBL driver fills and pbl_gas_apply uses.
program test_pbl_gas

   use mpas_kind_types, only: RKIND, StrKIND, I8KIND
   use unit_assert, only: check, finish
   use mpas_atmphys_pbl_gas
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   integer, parameter :: nz = 4
   logical :: supported, ok
   character(len=StrKIND) :: reason
   real(kind=RKIND) :: kzq(nz), m_d(nz), dz(nz), G(nz+1), q(2,nz), q_ref(2,nz), q0(2,nz), rho(nz), qv(nz)

   ! Capability gate (sketch §6.2; MYNN modes open stage by stage in P4, nondefault intervals wait for P5).
   call pbl_gas_capability('off', 'off', 'none', 0, 0, supported, reason)
   call check(supported, 'PBL off is a supported zero-transport case')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 1, 0, supported, reason)
   call check(supported, 'YSU with the Monin-Obukhov surface layer is supported')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov_rev', 'none', 1, 1, supported, reason)
   call check(supported, 'YSU with the revised Monin-Obukhov surface layer is supported (MYNN options ignored)')
   call pbl_gas_capability('bl_ysu', 'off', 'none', 0, 0, supported, reason)
   call check(.not. supported .and. index(reason, 'surface layer') > 0, 'active PBL without a surface layer is rejected')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', '00:30:00', 0, 0, supported, reason)
   call check(.not. supported .and. index(reason, 'config_pbl_interval') > 0, 'nondefault PBL interval is rejected')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 0, 0, supported, reason)
   call check(supported, 'MYNN local diffusion (config_mynn_edmf = 0) is supported')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 0, 1, supported, reason)
   call check(supported, 'MYNN with EDMF off is local diffusion whatever config_mynn_edmf_dd says')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 1, 0, supported, reason)
   call check(.not. supported .and. index(reason, 'MYNN updraft EDMF gas transport is not yet qualified') > 0, &
              'MYNN updraft EDMF is rejected until it is qualified')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', '00:30:00', 0, 0, supported, reason)
   call check(.not. supported .and. index(reason, 'config_pbl_interval') > 0, &
              'nondefault PBL interval is rejected with MYNN')
   call pbl_gas_capability('bl_other', 'sf_monin_obukhov', 'none', 0, 0, supported, reason)
   call check(.not. supported .and. len_trim(reason) > 0, 'unknown PBL scheme is rejected')

   ! Operator cache.
   call check(.not. pbl_gas_active(1), 'no column is active before the cache exists')
   call pbl_gas_cache_setup(nz, 3)
   kzq = 10.0_RKIND
   m_d = [120.0_RKIND, 100.0_RKIND, 80.0_RKIND, 60.0_RKIND]
   dz = [100.0_RKIND, 150.0_RKIND, 200.0_RKIND, 250.0_RKIND]
   call pbl_gas_cache_store_column(1, kzq, m_d, dz, ok)
   call check(ok .and. pbl_gas_active(1), 'a column with exchange is active')
   call pbl_gas_cache_store_column(2, 0.0_RKIND*kzq, m_d, dz, ok)
   call check(ok .and. .not. pbl_gas_active(2), 'a column without exchange is inactive')
   call check(.not. pbl_gas_active(3), 'a column never stored is inactive')

   q0(1,:) = [4.0e-9_RKIND, 1.0e-9_RKIND, 0.0_RKIND, 2.0e-9_RKIND]
   q0(2,:) = 0.0_RKIND
   q = q0
   call pbl_gas_apply(1, 30.0_RKIND, m_d, q)
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   q_ref = q0
   call pbl_diffuse_column(30.0_RKIND, m_d, G, q_ref)
   call check(all(transfer(q, [0_I8KIND]) == transfer(q_ref, [0_I8KIND])), &
              'apply is the diffusion kernel with the cached exchange (bitwise)')
   call check(relative_mass_residual(m_d, q0(1,:), q(1,:)) <= 1.0e-12_RKIND .and. all(q(2,:) == 0.0_RKIND), &
              'apply conserves mass and keeps a zero gas zero')

   ! The PBL driver hands over moist density rho = rho_dry*(1+qv) (MPAS rho_p, with qv >= 0); the
   ! stored operator uses the dry layer mass rho/(1+qv)*dz.
   rho = [1.20_RKIND, 1.10_RKIND, 1.00_RKIND, 0.90_RKIND]
   qv = [1.2e-2_RKIND, 8.0e-3_RKIND, 3.0e-3_RKIND, 0.0_RKIND]
   call pbl_gas_cache_store_moist_column(3, kzq, rho, qv, dz, ok)
   call check(ok .and. pbl_gas_active(3), 'a column stored from moist density is active')
   q = q0
   call pbl_gas_apply(3, 30.0_RKIND, rho/(1.0_RKIND + qv)*dz, q)
   call pbl_gas_exchange(kzq, rho/(1.0_RKIND + qv)*dz, dz, G, ok)
   q_ref = q0
   call pbl_diffuse_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, G, q_ref)
   call check(all(transfer(q, [0_I8KIND]) == transfer(q_ref, [0_I8KIND])), &
              'moist-density store uses the dry layer mass rho/(1+qv)*dz (bitwise)')
   call pbl_gas_cache_store_moist_column(3, kzq, rho, -1.0_RKIND + 0.0_RKIND*qv, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(3), 'a column with qv = -1 (no dry mass) is rejected')

   ! MYNN hands over exch_h on the bottom face of each layer (P4a); the store uses the dry density.
   call pbl_gas_cache_store_mynn_column(2, kzq, rho, qv, dz, ok)
   call check(ok .and. pbl_gas_active(2), 'a MYNN column with exchange is active')
   q = q0
   call pbl_gas_apply(2, 30.0_RKIND, rho/(1.0_RKIND + qv)*dz, q)
   call pbl_gas_exchange_faces(kzq, rho/(1.0_RKIND + qv), dz, G, ok)
   q_ref = q0
   call pbl_diffuse_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, G, q_ref)
   call check(all(transfer(q, [0_I8KIND]) == transfer(q_ref, [0_I8KIND])), &
              'MYNN store: bottom-face exchange with the dry density rho/(1+qv) (bitwise)')
   call pbl_gas_cache_store_mynn_column(2, -kzq, rho, qv, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(2), 'a MYNN column with negative exchange is rejected')

   kzq(2) = -1.0_RKIND
   call pbl_gas_cache_store_column(1, kzq, m_d, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(1), 'an invalid column is reported and left inactive')

   call pbl_gas_finalize()
   call check(.not. pbl_gas_active(1), 'finalize releases the cache')
   call check(.not. pbl_gas_enabled(), 'transport is not enabled without initialization')

   call finish('test_pbl_gas')

end program test_pbl_gas
