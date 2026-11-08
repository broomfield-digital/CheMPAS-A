! Run the production pbl_diffuse_column on the reference fixtures and a stiffness sweep, and write the results
! as text for plotting. Usage: dump_pbl_diffusion COLUMNS_DIR OUT_DIR
program dump_pbl_diffusion

   use mpas_kind_types, only: RKIND
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column

   implicit none

   character(len=*), parameter :: cases(6) = [character(len=24) :: 'nonuniform', 'pulses', 'sharp_gradient', &
                                              'zero_coefficients', 'large_timestep', 'uniform']
   character(len=512) :: cols, outdir
   real(kind=RKIND), allocatable :: m_d(:), G(:), q(:,:), q0(:,:)
   real(kind=RKIND) :: dt, kappa, resid, uniform_change
   character(len=256) :: comment
   integer :: ic, n, nz, nt, u, v, k, step, ip

   call get_command_argument(1, cols)
   call get_command_argument(2, outdir)

   do ic = 1, size(cases)
      open(newunit=u, file=trim(cols)//'/diffusion_'//trim(cases(ic))//'.txt', status='old', action='read')
      read(u, '(a)') comment
      read(u, *) nz, nt
      read(u, *) dt
      allocate(m_d(nz), G(nz+1), q(nt,nz))
      read(u, *) m_d
      read(u, *) G
      do n = 1, nt
         read(u, *) q(n,:)
      end do
      close(u)
      call pbl_diffuse_column(dt, m_d, G, q)
      open(newunit=v, file=trim(outdir)//'/kernel_'//trim(cases(ic))//'.txt', action='write')
      do n = 1, nt
         write(v, '(*(es25.17))') q(n,:)
      end do
      close(v)
      deallocate(m_d, G, q)
   end do

   ! Stiffness sweep on the 60-layer column of test_pbl_diffusion.
   nz = 60
   allocate(m_d(nz), G(nz+1), q(2,nz), q0(2,nz))
   do k = 1, nz
      m_d(k) = 50.0_RKIND + 150.0_RKIND*sin(0.37_RKIND*k)**2
      q0(1,k) = 1.0e-9_RKIND*(1.0_RKIND + sin(real(k, RKIND)))
      q0(2,k) = merge(5.0e-9_RKIND, 0.0_RKIND, k <= 5)
   end do
   dt = 600.0_RKIND
   open(newunit=v, file=trim(outdir)//'/kernel_stiffness.txt', action='write')
   write(v, '(a)') '# kappa  mass_residual_after_1000_calls  uniform_profile_max_relative_change'
   do ip = 0, 24
      kappa = 10.0_RKIND**(0.25_RKIND*ip)
      G = kappa*sum(m_d)/nz/dt
      q = q0
      do step = 1, 1000
         call pbl_diffuse_column(dt, m_d, G, q)
      end do
      resid = 0.0_RKIND
      do n = 1, 2
         resid = max(resid, abs(sum(m_d*q(n,:)) - sum(m_d*q0(n,:)))/sum(m_d*q0(n,:)))
      end do
      q(1,:) = 3.0e-9_RKIND
      q(2,:) = 3.0e-9_RKIND
      call pbl_diffuse_column(dt, m_d, G, q)
      uniform_change = maxval(abs(q - 3.0e-9_RKIND))/3.0e-9_RKIND
      write(v, '(3es25.17)') kappa, resid, uniform_change
   end do
   close(v)

end program dump_pbl_diffusion
