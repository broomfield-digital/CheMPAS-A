! Unit tests for the MYNN EDMF gas kernels (P4b of docs/plans/pbl-transport.md), on synthetic columns:
!   - pbl_edmf_updraft_fluxes: dry plume mass flux on faces from MYNN's per-plume UPA*UPW*Psig_w, and the
!     entrainment that makes each plume mass-consistent (D >= 0), matching MYNN's mixing fraction where it can;
!   - pbl_edmf_updraft_plume: the rebuilt plume concentration, equal to MYNN's linear mixing (1-eps)*c+eps*q
!     wherever MYNN's own fraction is mass-consistent;
!   - pbl_edmf_column: plume flux plus upwind compensating subsidence, substepped so that each layer update
!     is a convex combination: conservation, bounds, zero and uniform profiles, an independent layer-budget
!     form, substep equivalence and rejection of a bogus operator.
program test_pbl_edmf

   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mpas_kind_types, only: RKIND, I8KIND
   use unit_assert, only: check, finish
   use mpas_atmphys_pbl_edmf, only: pbl_edmf_updraft_fluxes, pbl_edmf_updraft_plume, pbl_edmf_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   integer, parameter :: nz = 30, np = 8, nt = 4
   real(RKIND) :: dz(nz), zc(nz), rho_d(nz), m_d(nz), aw(nz,np), eps(nz,np)
   real(RKIND) :: M(nz+1,np), E(nz,np), q0(nt,nz), q(nt,nz), qa(nt,nz), qb(nt,nz), c(nt,nz)
   real(RKIND) :: rhoz, resid, tol, worst, mtot(nz+1), etot(nz), cmax, dt, D
   logical :: ok, all_ok
   integer :: k, i, n, nsub, it, n_mynn

   call column(dz, zc, rho_d)
   m_d = rho_d*dz
   call plumes(aw, eps)

   ! ---- face mass flux and mass-consistent entrainment
   call pbl_edmf_updraft_fluxes(aw, eps, rho_d, dz, M, E, ok)
   call check(ok, 'valid plumes are accepted')
   all_ok = all(M(1,:) == 0.0_RKIND) .and. all(M(nz+1,:) == 0.0_RKIND)
   do k = 1, nz-1
      rhoz = (rho_d(k+1)*dz(k) + rho_d(k)*dz(k+1))/(dz(k) + dz(k+1))
      all_ok = all_ok .and. all(abs(M(k+1,:) - rhoz*aw(k,:)) <= 1.0e-15_RKIND*rhoz*aw(k,:))
   end do
   call check(all_ok, 'plume mass flux on face k+1 is the dry face density times aw(k); surface and top closed')
   all_ok = .true.
   do i = 1, np
      do k = 1, nz
         if (M(k+1,i) == 0.0_RKIND) then
            all_ok = all_ok .and. E(k,i) == 0.0_RKIND
         else
            D = M(k,i) + E(k,i) - M(k+1,i)
            all_ok = all_ok .and. D >= -1.0e-15_RKIND*M(k+1,i)
            all_ok = all_ok .and. abs(E(k,i) - max(eps(k,i)*M(k,i)/(1.0_RKIND - eps(k,i)), M(k+1,i) - M(k,i))) &
                     <= 1.0e-15_RKIND*(E(k,i) + M(k+1,i))
         end if
      end do
   end do
   call check(all_ok, 'entrainment is max(eps*M_in/(1-eps), M_out-M_in): detrainment >= 0; none where the plume ends')
   call check(all(E(1,:) == M(2,:)), 'every plume is rooted in layer 1 (all of its flux entrained there)')
   call check(any(E(2:nz,:) > 0.0_RKIND .and. abs(E(2:nz,:) - (M(3:nz+1,:) - M(2:nz,:))) <= 1.0e-15_RKIND*E(2:nz,:) &
                  .and. M(3:nz+1,:) > 0.0_RKIND .and. M(2:nz,:) > 0.0_RKIND), &
              'the test column has an accelerating plume whose entrainment is raised')

   ! ---- rebuilt plume concentration
   call profiles(zc, q0)
   all_ok = .true.
   n_mynn = 0
   do i = 1, np
      call pbl_edmf_updraft_plume(M(:,i), E(:,i), q0, c)
      all_ok = all_ok .and. all(c(:,1) == q0(:,1) .or. M(2,i) == 0.0_RKIND)
      do k = 2, nz
         if (M(k+1,i) > 0.0_RKIND .and. M(k,i) + E(k,i) - M(k+1,i) > 1.0e-12_RKIND*M(k+1,i) .and. &
             E(k,i) > max(0.0_RKIND, M(k+1,i) - M(k,i))) then
            ! MYNN's fraction is mass-consistent here: the plume is MYNN's linear mixing
            n_mynn = n_mynn + 1
            all_ok = all_ok .and. all(abs(c(:,k) - ((1.0_RKIND - eps(k,i))*c(:,k-1) + eps(k,i)*q0(:,k))) &
                                      <= 1.0e-14_RKIND*maxval(abs(q0)))
         end if
         if (M(k+1,i) == 0.0_RKIND) all_ok = all_ok .and. all(c(:,k) == 0.0_RKIND)
      end do
   end do
   write(*,'(a,i0)') 'plume layers where MYNN''s fraction is mass-consistent: ', n_mynn
   call check(all_ok .and. n_mynn > 0, 'rebuilt plume is MYNN''s linear mixing where MYNN''s fraction is mass-consistent')

   ! ---- column operator: one substep, against the layer-budget form
   mtot = sum(M, dim=2)
   etot = sum(E, dim=2)
   cmax = maxval((etot + mtot(1:nz))/m_d)
   dt = 0.9_RKIND/cmax
   q = q0
   call pbl_edmf_column(dt, m_d, M, E, q, nsub, ok)
   call check(ok .and. nsub == 1, 'a step within the outgoing-mass bound takes one substep')
   qb = q0
   do i = 1, np
      call pbl_edmf_updraft_plume(M(:,i), E(:,i), q0, c)
      do k = 1, nz
         D = M(k,i) + E(k,i) - M(k+1,i)
         if (M(k+1,i) > 0.0_RKIND) then
            qb(:,k) = qb(:,k) + dt/m_d(k)*(D*c(:,k) - E(k,i)*q0(:,k))
         else if (M(k,i) > 0.0_RKIND) then
            qb(:,k) = qb(:,k) + dt/m_d(k)*M(k,i)*c(:,k-1)    ! plume ends in layer k: all of it detrains
         end if
      end do
   end do
   do k = 1, nz
      qb(:,k) = qb(:,k) - dt/m_d(k)*mtot(k)*q0(:,k)
      if (k < nz) qb(:,k) = qb(:,k) + dt/m_d(k)*mtot(k+1)*q0(:,k+1)
   end do
   worst = 0.0_RKIND
   do n = 1, nt-1
      worst = max(worst, maxval(abs(q(n,:) - qb(n,:)))/maxval(q0(n,:)))
   end do
   write(*,'(a,es10.3)') 'worst difference from the layer-budget form, relative to the column maximum = ', worst
   call check(worst <= 1.0e-13_RKIND, 'flux form equals the independent layer-budget form')

   ! ---- column operator: properties with several substeps
   dt = 2.5_RKIND/cmax
   q = q0
   call pbl_edmf_column(dt, m_d, M, E, q, nsub, ok)
   call check(ok .and. nsub == 3, 'a step of 2.5 times the outgoing-mass bound takes three substeps')
   resid = 0.0_RKIND
   do n = 1, nt-1
      resid = max(resid, relative_mass_residual(m_d, q0(n,:), q(n,:)))
   end do
   call check(resid <= 1.0e-12_RKIND, 'every gas is conserved to roundoff')
   all_ok = .true.
   do n = 1, nt
      tol = 1.0e-20_RKIND + 1.0e-12_RKIND*maxval(q0(n,:))
      all_ok = all_ok .and. all(q(n,:) >= minval(q0(n,:)) - tol) .and. all(q(n,:) <= maxval(q0(n,:)) + tol)
   end do
   call check(all_ok, 'every gas stays within its initial range (nonnegative, no new extrema)')
   call check(all(q(nt,:) == 0.0_RKIND), 'a zero gas stays exactly zero')
   call check(q(1,1) < q0(1,1) .and. sum(m_d*zc*q(1,:)) > sum(m_d*zc*q0(1,:)), &
              'the plumes carry a surface-confined gas upward')
   qa = q0
   do it = 1, 3
      call pbl_edmf_column(dt/3.0_RKIND, m_d, M, E, qa, nsub, ok)
   end do
   call check(all(transfer(qa, [0_I8KIND]) == transfer(q, [0_I8KIND])), &
              'three substeps equal three steps of a third (bitwise)')
   q(1,:) = 4.0e-9_RKIND
   call pbl_edmf_column(dt, m_d, M, E, q(1:1,:), nsub, ok)
   call check(all(abs(q(1,:) - 4.0e-9_RKIND) <= 1.0e-20_RKIND + 1.0e-12_RKIND*4.0e-9_RKIND), &
              'a uniform gas is unchanged')
   q = q0
   do it = 1, 1000
      call pbl_edmf_column(dt, m_d, M, E, q, nsub, ok)
   end do
   resid = 0.0_RKIND
   do n = 1, nt-1
      resid = max(resid, relative_mass_residual(m_d, q0(n,:), q(n,:)))
   end do
   call check(resid <= 1.0e-10_RKIND .and. all(q(nt,:) == 0.0_RKIND), 'conserved over 1,000 calls')
   q = q0
   call pbl_edmf_column(1.0e9_RKIND, m_d, M, E, q, nsub, ok)
   call check(.not. ok .and. all(transfer(q, [0_I8KIND]) == transfer(q0, [0_I8KIND])), &
              'a step beyond 1e4 times the outgoing-mass bound is rejected and leaves the gases unchanged')

   ! ---- invalid plume data
   call expect_invalid(-1.0_RKIND, 0.1_RKIND, 'negative aw is rejected')
   call expect_invalid(0.01_RKIND, 1.0_RKIND, 'an entrainment fraction of 1 is rejected')
   call expect_invalid(0.01_RKIND, -0.1_RKIND, 'a negative entrainment fraction is rejected')
   call expect_invalid(ieee_value(1.0_RKIND, ieee_quiet_nan), 0.1_RKIND, 'non-finite aw is rejected')
   aw(nz,1) = 0.1_RKIND
   call pbl_edmf_updraft_fluxes(aw, eps, rho_d, dz, M, E, ok)
   call check(.not. ok .and. all(M == 0.0_RKIND) .and. all(E == 0.0_RKIND), 'flux through the model top is rejected')

   call finish('test_pbl_edmf')

contains

   subroutine column(dz, zc, rho_d)
      real(RKIND), intent(out) :: dz(nz), zc(nz), rho_d(nz)
      real(RKIND) :: z
      integer :: k
      z = 0.0_RKIND
      do k = 1, nz
         dz(k) = 40.0_RKIND + 15.0_RKIND*(k-1)
         zc(k) = z + 0.5_RKIND*dz(k)
         z = z + dz(k)
         rho_d(k) = 1.2_RKIND*exp(-zc(k)/8000.0_RKIND)
      end do
   end subroutine column

   subroutine plumes(aw, eps)
      ! Eight plumes ending at different heights. Velocity (and so aw) first increases with height, the
      ! accelerating case where MYNN's fraction is not mass-consistent, and then decreases.
      real(RKIND), intent(out) :: aw(nz,np), eps(nz,np)
      integer :: i, k, ktop
      aw = 0.0_RKIND
      eps = 0.0_RKIND
      do i = 1, np
         ktop = 6 + 2*i
         do k = 1, ktop
            aw(k,i) = 0.004_RKIND*(1.0_RKIND + 0.5_RKIND*sin(3.0_RKIND*k/ktop))*(1.0_RKIND + 0.1_RKIND*i)
            if (k > 1) eps(k,i) = min(0.85_RKIND, 0.02_RKIND*(1.0_RKIND + mod(3*k + i, 5)))
         end do
      end do
   end subroutine plumes

   subroutine profiles(zc, q0)
      real(RKIND), intent(in) :: zc(nz)
      real(RKIND), intent(out) :: q0(nt,nz)
      integer :: k
      q0(1,:) = 2.0e-9_RKIND*exp(-zc/300.0_RKIND) + 1.0e-12_RKIND
      q0(2,:) = merge(1.0e-9_RKIND, 1.0e-11_RKIND, zc < 1000.0_RKIND)
      do k = 1, nz
         q0(3,k) = 5.0e-8_RKIND*(1.0_RKIND + 0.6_RKIND*sin(1.7_RKIND*k))
      end do
      q0(4,:) = 0.0_RKIND
   end subroutine profiles

   subroutine expect_invalid(bad_aw, bad_eps, label)
      real(RKIND), intent(in) :: bad_aw, bad_eps
      character(len=*), intent(in) :: label
      real(RKIND) :: aw2(nz,np), eps2(nz,np)
      logical :: ok2
      call plumes(aw2, eps2)
      aw2(3,2) = bad_aw
      eps2(4,5) = bad_eps
      if (bad_eps >= 0.0_RKIND .and. bad_eps < 1.0_RKIND) eps2(4,5) = 0.1_RKIND
      call pbl_edmf_updraft_fluxes(aw2, eps2, rho_d, dz, M, E, ok2)
      call check(.not. ok2 .and. all(M == 0.0_RKIND) .and. all(E == 0.0_RKIND), label)
   end subroutine expect_invalid

end program test_pbl_edmf
