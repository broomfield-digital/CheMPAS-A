! Unit tests for the PBL gas process (P3): the capability gate and the per-column operator cache
! that the PBL driver fills and pbl_gas_apply uses.
program test_pbl_gas

   use mpas_kind_types, only: RKIND, StrKIND, I8KIND
   use unit_assert, only: check, finish
   use mpas_atmphys_pbl_gas
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column
   use mpas_atmphys_pbl_edmf, only: pbl_edmf_updraft_fluxes, pbl_edmf_downdraft_fluxes, pbl_edmf_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   integer, parameter :: nz = 4
   logical :: supported, ok
   character(len=StrKIND) :: reason
   real(kind=RKIND) :: kzq(nz), m_d(nz), dz(nz), G(nz+1), q(2,nz), q_ref(2,nz), q0(2,nz), rho(nz), qv(nz)
   real(kind=RKIND) :: aw(nz,2), ent(nz,2), Mu(nz+1,2), Eu(nz,2), awd(nz,1), entd(nz,1), Nd(nz+1,1), Ed(nz,1)
   integer :: nsub

   ! Capability gate (sketch §6.2; any PBL interval since P5).
   call pbl_gas_capability('off', 'off', 'none', 'off', 1, supported, reason)
   call check(supported, 'PBL off is a supported zero-transport case')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'off', 1, supported, reason)
   call check(supported, 'YSU with the Monin-Obukhov surface layer is supported')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov_rev', 'none', 'off', 1, supported, reason)
   call check(supported, 'YSU with the revised Monin-Obukhov surface layer is supported')
   call pbl_gas_capability('bl_ysu', 'off', 'none', 'off', 1, supported, reason)
   call check(.not. supported .and. index(reason, 'surface layer') > 0, 'active PBL without a surface layer is rejected')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', '00:00:09', 'off', 1, supported, reason)
   call check(supported, 'YSU at a nondefault PBL interval is supported (held operator, P5)')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 'off', 1, supported, reason)
   call check(supported, 'MYNN is supported (local diffusion, updrafts and downdrafts)')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', '00:30:00', 'off', 1, supported, reason)
   call check(supported, 'MYNN at a nondefault PBL interval is supported (held operator, P5)')
   call pbl_gas_capability('bl_other', 'sf_monin_obukhov', 'none', 'off', 1, supported, reason)
   call check(.not. supported .and. len_trim(reason) > 0, 'unknown PBL scheme is rejected')

   ! Parameterized convection carries no MICM gas yet (#41), and OpenMP threading is not qualified (P5).
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'cu_ntiedtke', 1, supported, reason)
   call check(.not. supported .and. index(reason, '#41') > 0, 'YSU with new Tiedtke convection is rejected (#41)')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'cu_tiedtke', 1, supported, reason)
   call check(.not. supported .and. index(reason, '#41') > 0, 'YSU with Tiedtke convection is rejected')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'cu_kain_fritsch', 1, supported, reason)
   call check(.not. supported .and. index(reason, '#41') > 0, 'YSU with Kain-Fritsch convection is rejected')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 'cu_grell_freitas', 1, supported, reason)
   call check(.not. supported .and. index(reason, '#41') > 0, &
              'MYNN with Grell-Freitas (the convection_permitting suite) is rejected')
   call pbl_gas_capability('bl_mynn', 'sf_mynn', 'none', 'off', 1, supported, reason)
   call check(supported, 'convection off is supported (an explicit off under any suite)')
   call pbl_gas_capability('off', 'off', 'none', 'cu_ntiedtke', 1, supported, reason)
   call check(supported, 'PBL off with convection is a supported zero-transport case')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'off', 2, supported, reason)
   call check(.not. supported .and. index(reason, 'OpenMP') > 0, 'two OpenMP threads are rejected')
   call pbl_gas_capability('bl_ysu', 'sf_monin_obukhov', 'none', 'off', 1, supported, reason)
   call check(supported, 'one thread is supported')

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

   ! MYNN updraft plumes (P4b): a cache with two plumes per column.
   call pbl_gas_cache_setup(nz, 3, nup=2)
   aw(:,1) = [0.010_RKIND, 0.012_RKIND, 0.006_RKIND, 0.0_RKIND]
   aw(:,2) = [0.008_RKIND, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND]
   ent(:,1) = [0.0_RKIND, 0.1_RKIND, 0.2_RKIND, 0.0_RKIND]
   ent(:,2) = 0.0_RKIND
   call pbl_gas_cache_store_mynn_column(1, 0.0_RKIND*kzq, rho, qv, dz, ok)
   call check(ok .and. .not. pbl_gas_active(1), 'a MYNN column without exchange or plumes is inactive')
   call pbl_gas_cache_store_mynn_plumes(1, aw, ent, rho, qv, dz, ok)
   call check(ok .and. pbl_gas_active(1), 'a MYNN column with plumes but no exchange is active')
   call pbl_gas_cache_store_mynn_column(2, kzq, rho, qv, dz, ok)
   call pbl_gas_cache_store_mynn_plumes(2, aw, ent, rho, qv, dz, ok)
   q = q0
   call pbl_gas_apply(2, 30.0_RKIND, rho/(1.0_RKIND + qv)*dz, q)
   call pbl_gas_exchange_faces(kzq, rho/(1.0_RKIND + qv), dz, G, ok)
   q_ref = q0
   call pbl_diffuse_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, G, q_ref)
   call pbl_edmf_updraft_fluxes(aw, ent, rho/(1.0_RKIND + qv), dz, Mu, Eu, ok)
   call pbl_edmf_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, Mu, Eu, Nd(:,1:0), Ed(:,1:0), q_ref, nsub, ok)
   call check(all(transfer(q, [0_I8KIND]) == transfer(q_ref, [0_I8KIND])), &
              'MYNN apply is local diffusion, then the updraft plumes, with the dry density (bitwise)')
   call check(any(q(1,:) /= q0(1,:)) .and. all(q(2,:) == 0.0_RKIND), 'the plumes act; a zero gas stays zero')
   aw(2,1) = -1.0_RKIND
   call pbl_gas_cache_store_mynn_plumes(2, aw, ent, rho, qv, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(2), 'a column with invalid plumes is reported and left inactive')
   aw(2,1) = 0.012_RKIND

   ! MYNN downdrafts (P4c): one downdraft from layer 3 to the face of layer 2.
   call pbl_gas_cache_setup(nz, 3, nup=2, ndown=1)
   awd(:,1) = [0.0_RKIND, 0.005_RKIND, 0.004_RKIND, 0.0_RKIND]
   entd(:,1) = [0.0_RKIND, 0.1_RKIND, 0.0_RKIND, 0.0_RKIND]
   call pbl_gas_cache_store_mynn_column(1, 0.0_RKIND*kzq, rho, qv, dz, ok)
   call pbl_gas_cache_store_mynn_plumes(1, 0.0_RKIND*aw, ent, rho, qv, dz, ok, awd, entd)
   call check(ok .and. pbl_gas_active(1), 'a MYNN column with only a downdraft is active')
   call pbl_gas_cache_store_mynn_column(2, kzq, rho, qv, dz, ok)
   call pbl_gas_cache_store_mynn_plumes(2, aw, ent, rho, qv, dz, ok, awd, entd)
   q = q0
   call pbl_gas_apply(2, 30.0_RKIND, rho/(1.0_RKIND + qv)*dz, q)
   call pbl_gas_exchange_faces(kzq, rho/(1.0_RKIND + qv), dz, G, ok)
   q_ref = q0
   call pbl_diffuse_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, G, q_ref)
   call pbl_edmf_updraft_fluxes(aw, ent, rho/(1.0_RKIND + qv), dz, Mu, Eu, ok)
   call pbl_edmf_downdraft_fluxes(awd, entd, rho/(1.0_RKIND + qv), dz, Nd, Ed, ok)
   call pbl_edmf_column(30.0_RKIND, rho/(1.0_RKIND + qv)*dz, Mu, Eu, Nd, Ed, q_ref, nsub, ok)
   call check(all(transfer(q, [0_I8KIND]) == transfer(q_ref, [0_I8KIND])), &
              'MYNN apply is local diffusion, then updrafts and downdrafts together (bitwise)')
   call pbl_gas_cache_store_mynn_plumes(2, aw, ent, rho, qv, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(2), 'a downdraft cache without downdraft data is reported')
   awd(1,1) = 0.001_RKIND
   call pbl_gas_cache_store_mynn_plumes(2, aw, ent, rho, qv, dz, ok, awd, entd)
   call check(.not. ok .and. .not. pbl_gas_active(2), 'a downdraft into the ground is reported and left inactive')
   call pbl_gas_cache_setup(nz, 3)

   kzq(2) = -1.0_RKIND
   call pbl_gas_cache_store_column(1, kzq, m_d, dz, ok)
   call check(.not. ok .and. .not. pbl_gas_active(1), 'an invalid column is reported and left inactive')

   ! Held operator (P5): valid for the steps before the next refresh, judged by its age in model steps.
   call check(pbl_gas_operator_valid(0, 3.0_RKIND, 3.0_RKIND), 'default interval: valid in the refresh step')
   call check(.not. pbl_gas_operator_valid(1, 3.0_RKIND, 3.0_RKIND), 'default interval: stale one step later')
   call check(pbl_gas_operator_valid(2, 3.0_RKIND, 9.0_RKIND), '9 s interval: valid two steps after refresh')
   call check(.not. pbl_gas_operator_valid(3, 3.0_RKIND, 9.0_RKIND), '9 s interval: stale three steps after')
   call check(pbl_gas_operator_valid(3, 3.0_RKIND, 10.0_RKIND), '10 s interval: valid three steps after refresh')
   call check(.not. pbl_gas_operator_valid(4, 3.0_RKIND, 10.0_RKIND), '10 s interval: stale four steps after')
   call check(.not. pbl_gas_operator_valid(-1, 3.0_RKIND, 9.0_RKIND), 'never refreshed: not valid')

   call pbl_gas_finalize()
   call check(pbl_gas_refresh_time() == -1.0_RKIND, 'no refresh time before the first refresh')
   call pbl_gas_cache_refreshed(27.0_RKIND, 9.0_RKIND)
   call check(pbl_gas_refresh_time() == 27.0_RKIND, 'the refresh records the model time')
   call pbl_gas_begin_step(3.0_RKIND)
   call pbl_gas_begin_step(3.0_RKIND)
   call pbl_gas_begin_step(3.0_RKIND)
   call check(.true., 'three steps of a 9 s operator after a refresh are accepted')

   call pbl_gas_finalize()
   call check(pbl_gas_refresh_time() == -1.0_RKIND, 'finalize clears the refresh stamp')
   call check(.not. pbl_gas_active(1), 'finalize releases the cache')
   call check(.not. pbl_gas_enabled(), 'transport is not enabled without initialization')

   call finish('test_pbl_gas')

end program test_pbl_gas
