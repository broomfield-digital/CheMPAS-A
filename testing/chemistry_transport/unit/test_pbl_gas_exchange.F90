! Unit tests for the MICM gas PBL operator construction (P3): interface exchange
! G(k+1) = 0.5*(m_d(k)+m_d(k+1))*K(k)/dzc**2 = rho_dry_interface*K/dzc, with K(k)
! the scheme's scalar exchange coefficient between layers k and k+1.
program test_pbl_gas_exchange

   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mpas_kind_types, only: RKIND
   use unit_assert, only: check, finish
   use mpas_atmphys_pbl_gas, only: pbl_gas_exchange

   implicit none

   integer, parameter :: nz = 5
   real(kind=RKIND) :: kzq(nz), m_d(nz), dz(nz), G(nz+1), expect
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

   call finish('test_pbl_gas_exchange')

end program test_pbl_gas_exchange
