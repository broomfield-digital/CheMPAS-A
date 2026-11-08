! Unit tests for the MICM gas PBL operator construction (P3): interface exchange
! G(k+1) = 0.5*(m_d(k)+m_d(k+1))*K(k)/dzc**2 = rho_dry_interface*K/dzc, with K(k)
! the scheme's scalar exchange coefficient between layers k and k+1 (YSU, P3); and the MYNN form (P4a)
! G(k) = rhoz_d(k)*K(k)/dzk with K(k) on the bottom face of layer k and MYNN's own face density
! rhoz_d(k) = (rho_d(k)*dz(k-1)+rho_d(k-1)*dz(k))/(dz(k-1)+dz(k)).
program test_pbl_gas_exchange

   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mpas_kind_types, only: RKIND
   use unit_assert, only: check, finish
   use mpas_atmphys_pbl_gas, only: pbl_gas_exchange, pbl_gas_exchange_faces

   implicit none

   integer, parameter :: nz = 5
   real(kind=RKIND) :: kzq(nz), m_d(nz), dz(nz), G(nz+1), expect, kh(nz), rho_d(nz), rhoz, dzk
   logical :: ok, all_ok
   integer :: k

   kzq = 10.0_RKIND
   m_d = 100.0_RKIND
   dz = 200.0_RKIND
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   call check(ok, 'valid uniform column is accepted')
   call check(G(1) == 0.0_RKIND .and. G(nz+1) == 0.0_RKIND, 'ground and top interfaces are closed')
   call check(all(abs(G(2:nz) - 0.025_RKIND) <= 1.0e-15_RKIND), 'uniform column: G = rho*K/dz')

   kzq = [1.0_RKIND, 2.0_RKIND, 3.0_RKIND, 4.0_RKIND, 99.0_RKIND]
   m_d = [120.0_RKIND, 100.0_RKIND, 80.0_RKIND, 60.0_RKIND, 40.0_RKIND]
   dz  = [100.0_RKIND, 150.0_RKIND, 200.0_RKIND, 250.0_RKIND, 300.0_RKIND]
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   all_ok = ok
   do k = 1, nz-1
      expect = 0.5_RKIND*(m_d(k) + m_d(k+1))*kzq(k)/(0.5_RKIND*(dz(k) + dz(k+1)))**2
      all_ok = all_ok .and. abs(G(k+1) - expect) <= 1.0e-14_RKIND*expect
   end do
   call check(all_ok, 'nonuniform column: interface k+1 uses K(k), the mean dry mass and the center spacing')
   call check(G(nz+1) == 0.0_RKIND, 'K above the top interface is ignored')

   kzq = 10.0_RKIND
   m_d = 100.0_RKIND
   dz = 200.0_RKIND
   kzq(2) = -1.0_RKIND
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   call check(.not. ok, 'negative exchange coefficient is rejected')
   kzq(2) = ieee_value(1.0_RKIND, ieee_quiet_nan)
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   call check(.not. ok, 'non-finite exchange coefficient is rejected')
   kzq(2) = 10.0_RKIND
   m_d(3) = 0.0_RKIND
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   call check(.not. ok, 'nonpositive layer mass is rejected')
   m_d(3) = 100.0_RKIND
   dz(4) = -5.0_RKIND
   call pbl_gas_exchange(kzq, m_d, dz, G, ok)
   call check(.not. ok, 'nonpositive layer thickness is rejected')

   ! MYNN face convention (exch_h(k) on the bottom face of layer k; exch_h(1) is the surface).
   kh = [5.0_RKIND, 1.0_RKIND, 2.0_RKIND, 3.0_RKIND, 4.0_RKIND]
   rho_d = [1.2_RKIND, 1.1_RKIND, 1.0_RKIND, 0.9_RKIND, 0.8_RKIND]
   dz  = [100.0_RKIND, 150.0_RKIND, 200.0_RKIND, 250.0_RKIND, 300.0_RKIND]
   call pbl_gas_exchange_faces(kh, rho_d, dz, G, ok)
   all_ok = ok .and. G(1) == 0.0_RKIND .and. G(nz+1) == 0.0_RKIND
   do k = 2, nz
      dzk = 0.5_RKIND*(dz(k) + dz(k-1))
      rhoz = (rho_d(k)*dz(k-1) + rho_d(k-1)*dz(k))/(dz(k-1) + dz(k))
      expect = rhoz*kh(k)/dzk
      all_ok = all_ok .and. abs(G(k) - expect) <= 1.0e-14_RKIND*expect
   end do
   call check(all_ok, 'MYNN: face k uses K(k), MYNN''s face density and the center spacing; surface K ignored')
   kh(3) = -1.0_RKIND
   call pbl_gas_exchange_faces(kh, rho_d, dz, G, ok)
   call check(.not. ok .and. all(G == 0.0_RKIND), 'MYNN: negative exchange coefficient is rejected')
   kh(3) = ieee_value(1.0_RKIND, ieee_quiet_nan)
   call pbl_gas_exchange_faces(kh, rho_d, dz, G, ok)
   call check(.not. ok, 'MYNN: non-finite exchange coefficient is rejected')
   kh(3) = 2.0_RKIND
   rho_d(2) = 0.0_RKIND
   call pbl_gas_exchange_faces(kh, rho_d, dz, G, ok)
   call check(.not. ok, 'MYNN: nonpositive density is rejected')
   rho_d(2) = 1.1_RKIND
   kh(1) = -7.0_RKIND
   call pbl_gas_exchange_faces(kh, rho_d, dz, G, ok)
   call check(ok .and. G(1) == 0.0_RKIND, 'MYNN: the surface value is ignored')

   call finish('test_pbl_gas_exchange')

end program test_pbl_gas_exchange
