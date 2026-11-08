! Unit tests for the closed-column implicit diffusion kernel (P2): agreement
! with exact reference solutions (columns/diffusion_*.txt, from
! make_diffusion_fixtures.py), conservation, monotonicity, and the analytic
! constant-coefficient decay of discrete cosine modes. Agreement with the
! reference is required to max(1e-12, 32*eps*kappa) relative, where
! kappa = max(dt*(G(k)+G(k+1))/m_d(k)): roundoff in any backward-stable solve
! grows with this conditioning. The frozen thresholds (mass, uniform profile,
! nonnegativity) are checked at their fixed values.
program test_pbl_diffusion

   use mpas_kind_types, only: RKIND, I8KIND, StrKIND
   use unit_assert, only: check, finish
   use mpas_atmphys_tracer_column, only: relative_mass_residual
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column

   implicit none

   character(len=*), parameter :: cases(6) = [character(len=24) :: 'nonuniform', 'pulses', 'sharp_gradient', &
                                              'zero_coefficients', 'large_timestep', 'uniform']
   character(len=StrKIND) :: dir
   real(kind=RKIND), allocatable :: m_d(:), G(:), q0(:,:), qref(:,:), q(:,:), qp(:,:), G2(:)
   real(kind=RKIND) :: dt, tol, worst, floor_value, kappa
   integer :: ic, n, nz, nt, j, k, step, status
   logical :: ok

   call get_environment_variable('CHEMPAS_UNIT_COLUMNS', dir, status=status)
   call check(status == 0, 'fixture directory is set (CHEMPAS_UNIT_COLUMNS)')
   if (status /= 0) call finish('test_pbl_diffusion')

   do ic = 1, size(cases)
      call read_case(trim(dir)//'/diffusion_'//trim(cases(ic))//'.txt', ok)
      call check(ok, trim(cases(ic))//': reference fixture is readable')
      if (.not. ok) cycle
      q = q0
      call pbl_diffuse_column(dt, m_d, G, q)

      kappa = 0.0_RKIND
      do k = 1, nz
         kappa = max(kappa, dt*(merge(G(k), 0.0_RKIND, k > 1) + merge(G(k+1), 0.0_RKIND, k < nz))/m_d(k))
      end do
      ok = .true.
      do n = 1, nt
         tol = max(1.0e-12_RKIND, 32.0_RKIND*epsilon(1.0_RKIND)*kappa)*maxval(abs(qref(n,:))) + 1.0e-40_RKIND
         ok = ok .and. all(abs(q(n,:) - qref(n,:)) <= tol)
      end do
      call check(ok, trim(cases(ic))//': matches the exact dense solve')

      worst = 0.0_RKIND
      do n = 1, nt
         worst = max(worst, relative_mass_residual(m_d, q0(n,:), q(n,:)))
      end do
      call check(worst <= 1.0e-12_RKIND, trim(cases(ic))//': mass residual per call <= 1e-12')

      ok = .true.
      do n = 1, nt
         floor_value = -(1.0e-20_RKIND + 1.0e-12_RKIND*maxval(q0(n,:)))
         ok = ok .and. all(q(n,:) >= floor_value)
         if (all(q0(n,:) == 0.0_RKIND)) ok = ok .and. all(q(n,:) == 0.0_RKIND)
      end do
      call check(ok, trim(cases(ic))//': nonnegative, zero species stay exactly zero')

      ok = .true.
      do n = 1, nt
         ok = ok .and. minval(q(n,:)) >= minval(q0(n,:)) - 1.0e-12_RKIND*maxval(abs(q0(n,:))) &
                 .and. maxval(q(n,:)) <= maxval(q0(n,:)) + 1.0e-12_RKIND*maxval(abs(q0(n,:)))
      end do
      call check(ok, trim(cases(ic))//': maximum principle')

      ! Species order does not change any result bit.
      qp = q0(nt:1:-1,:)
      call pbl_diffuse_column(dt, m_d, G, qp)
      call check(bitwise_equal(qp(nt:1:-1,:), q), trim(cases(ic))//': species permutation is bitwise invariant')

      ! Values supplied at the ground and top interfaces are ignored (closed column).
      G2 = G
      G2(1) = 5.0_RKIND
      G2(nz+1) = 5.0_RKIND
      qp = q0
      call pbl_diffuse_column(dt, m_d, G2, qp)
      call check(bitwise_equal(qp, q), trim(cases(ic))//': ground and top interfaces are closed')

      if (trim(cases(ic)) == 'uniform') then
         call check(all(abs(q - q0) <= 1.0e-20_RKIND + 1.0e-12_RKIND*q0), 'uniform profile is unchanged')
      end if

      if (trim(cases(ic)) == 'nonuniform') then
         q = q0
         do step = 1, 1000
            call pbl_diffuse_column(dt, m_d, G, q)
         end do
         worst = 0.0_RKIND
         do n = 1, nt
            worst = max(worst, relative_mass_residual(m_d, q0(n,:), q(n,:)))
         end do
         call check(worst <= 1.0e-10_RKIND, 'mass residual after 1000 calls <= 1e-10')
      end if
   end do

   ! Extreme stiffness (dt*G/m_d = 1e6, beyond any reference solve): conservation must not degrade with
   ! stiffness, per call or over 1000 calls, and the result stays nonnegative.
   nz = 60
   deallocate(m_d, G, q, q0)
   allocate(m_d(nz), G(nz+1), q(2,nz), q0(2,nz))
   do k = 1, nz
      m_d(k) = 50.0_RKIND + 150.0_RKIND*sin(0.37_RKIND*k)**2
      q0(1,k) = 1.0e-9_RKIND*(1.0_RKIND + sin(real(k, RKIND)))
      q0(2,k) = merge(5.0e-9_RKIND, 0.0_RKIND, k <= 5)
   end do
   dt = 600.0_RKIND
   G = 1.0e6_RKIND*sum(m_d)/nz/dt
   q = q0
   call pbl_diffuse_column(dt, m_d, G, q)
   worst = max(relative_mass_residual(m_d, q0(1,:), q(1,:)), relative_mass_residual(m_d, q0(2,:), q(2,:)))
   call check(worst <= 1.0e-12_RKIND, 'extreme stiffness: mass residual per call <= 1e-12')
   call check(all(q >= -(1.0e-20_RKIND + 1.0e-12_RKIND*5.0e-9_RKIND)), 'extreme stiffness: nonnegative')
   do step = 2, 1000
      call pbl_diffuse_column(dt, m_d, G, q)
   end do
   worst = max(relative_mass_residual(m_d, q0(1,:), q(1,:)), relative_mass_residual(m_d, q0(2,:), q(2,:)))
   call check(worst <= 1.0e-10_RKIND, 'extreme stiffness: mass residual after 1000 calls <= 1e-10')

   ! Stiff columns keep the shape properties too: a uniform profile is unchanged within the frozen uniform
   ! tolerance, and no value leaves the initial range (maximum principle), at dt*G/m_d = 1e4 and 1e6.
   do j = 1, 2
      G = merge(1.0e4_RKIND, 1.0e6_RKIND, j == 1)*sum(m_d)/nz/dt
      q(1,:) = 3.0e-9_RKIND
      q(2,:) = q0(2,:)
      call pbl_diffuse_column(dt, m_d, G, q)
      call check(all(abs(q(1,:) - 3.0e-9_RKIND) <= 1.0e-20_RKIND + 1.0e-12_RKIND*3.0e-9_RKIND), &
                 trim(merge('dt*G/m_d = 1e4', 'dt*G/m_d = 1e6', j == 1))//': uniform profile is unchanged')
      call check(minval(q(2,:)) >= -1.0e-12_RKIND*5.0e-9_RKIND .and. &
                 maxval(q(2,:)) <= 5.0e-9_RKIND*(1.0_RKIND + 1.0e-12_RKIND), &
                 trim(merge('dt*G/m_d = 1e4', 'dt*G/m_d = 1e6', j == 1))//': maximum principle')
   end do
   deallocate(q0)

   ! Constant coefficients: the closed-column cosine modes are eigenvectors, and backward Euler
   ! multiplies mode j by 1/(1 + dt*G/m*4*sin(pi*j/(2*nz))**2). The mean (j = 0) is unchanged.
   nz = 16
   deallocate(m_d, G, q)
   allocate(m_d(nz), G(nz+1), q(3,nz))
   m_d = 100.0_RKIND
   G = 0.7_RKIND
   G(1) = 0.0_RKIND
   G(nz+1) = 0.0_RKIND
   dt = 50.0_RKIND
   do j = 1, 3
      do k = 1, nz
         q(j,k) = 2.0_RKIND + cos(acos(-1.0_RKIND)*j*(k - 0.5_RKIND)/nz)
      end do
   end do
   call pbl_diffuse_column(dt, m_d, G, q)
   ok = .true.
   do j = 1, 3
      do k = 1, nz
         ok = ok .and. abs(q(j,k) - (2.0_RKIND + cos(acos(-1.0_RKIND)*j*(k - 0.5_RKIND)/nz) / &
              (1.0_RKIND + dt*0.7_RKIND/100.0_RKIND*4.0_RKIND*sin(acos(-1.0_RKIND)*j/(2.0_RKIND*nz))**2))) &
              <= 1.0e-12_RKIND*3.0_RKIND
      end do
   end do
   call check(ok, 'cosine modes decay at the analytic backward-Euler rate')

   call finish('test_pbl_diffusion')

contains

   subroutine read_case(path, found)
      character(len=*), intent(in) :: path
      logical, intent(out) :: found
      character(len=256) :: comment
      integer :: u, row, ios

      if (allocated(m_d)) deallocate(m_d, G, q0, qref)
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      found = ios == 0
      if (.not. found) return
      read(u, '(a)') comment
      read(u, *) nz, nt
      read(u, *) dt
      allocate(m_d(nz), G(nz+1), q0(nt,nz), qref(nt,nz))
      read(u, *) m_d
      read(u, *) G
      do row = 1, nt
         read(u, *) q0(row,:)
      end do
      do row = 1, nt
         read(u, *) qref(row,:)
      end do
      close(u)
   end subroutine read_case

   logical function bitwise_equal(a, b)
      real(kind=RKIND), intent(in) :: a(:,:), b(:,:)

      bitwise_equal = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
   end function bitwise_equal

end program test_pbl_diffusion
