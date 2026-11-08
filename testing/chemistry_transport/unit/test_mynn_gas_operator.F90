! Column reference test for the MYNN gas operator (P4 of docs/plans/pbl-transport.md).
!
! Synthetic stable and convective columns run through the MPAS MYNN wrapper (mynn_column.F90).
! P4a, local diffusion (config_mynn_edmf = 0):
!   A. In a dry column (qv = 0), the CheMPAS-A operator built from MYNN's exported exch_h (bottom face
!      of layer k) with MYNN's face density (pbl_gas_exchange_faces) and the P2 backward-Euler kernel at
!      MYNN's dt reproduces MYNN's own mixing of a passive scalar (water-friendly aerosol, uniform above
!      3 km) to roundoff: exch_h is the operator MYNN applies to scalars. (MYNN holds its top layer
!      fixed; the flat top makes that boundary identical to the closed one. With EDMF off MYNN's aerosol
!      mixing is local diffusion only.)
!   B. The operator conserves dry mass of every gas to roundoff, keeps a zero gas zero, leaves a
!      uniform gas unchanged and mixes a surface-confined gas upward.
program test_mynn_gas_operator

   use mpas_kind_types, only: RKIND
   use unit_assert, only: check, finish
   use mynn_column, only: column_state, column_output, make_column, run_mynn, nz
   use mpas_atmphys_pbl_gas, only: pbl_gas_exchange_faces
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   real(RKIND), parameter :: delt = 30.0_RKIND
   real(RKIND) :: gsum_stable, gsum_convective

   call local_case('stable', gsum_stable)
   call local_case('convective', gsum_convective)
   call check(gsum_convective > gsum_stable, 'convective MYNN exchange exceeds the stable exchange')

   call finish('test_mynn_gas_operator')

contains

   subroutine local_case(regime, gsum)
      character(len=*), intent(in) :: regime
      real(RKIND), intent(out) :: gsum

      type(column_state) :: col
      type(column_output) :: out
      real(RKIND) :: ozone(nz), aerosol(nz), rnwfa(nz), rho_d(nz), m_d(nz), G(nz+1), x_mynn(nz), oz(1,nz)
      real(RKIND) :: q0(4,nz), q(4,nz), u(1,nz), worst, resid
      logical :: ok
      integer :: n

      call make_column(regime, .false., col)
      ozone = 5.0e-8_RKIND*(0.5_RKIND + 0.5_RKIND*min(col%zc/3000.0_RKIND, 1.0_RKIND))
      aerosol = 1.0e9_RKIND*(1.0_RKIND + max(0.0_RKIND, 1.0_RKIND - col%zc/3000.0_RKIND))
      call run_mynn(col, 0, 0, 0, delt, out, nwfa=aerosol, rnwfa=rnwfa)
      call check(out%errflg == 0, regime//': MYNN succeeds')
      call check(out%exch_h(1) == 0.0_RKIND .and. all(out%exch_h >= 0.0_RKIND) .and. &
                 any(out%exch_h > 0.0_RKIND), regime//': exch_h is zero at the surface face and nonnegative')

      rho_d = col%rho/(1.0_RKIND + col%qv)
      m_d = rho_d*col%dz
      call pbl_gas_exchange_faces(out%exch_h, rho_d, col%dz, G, ok)
      call check(ok, regime//': MYNN exchange builds a valid operator')
      gsum = sum(G)

      ! A. MYNN's own scalar operator, applied to the aerosol
      oz(1,:) = aerosol
      call pbl_diffuse_column(delt, m_d, G, oz)
      x_mynn = aerosol + delt*rnwfa
      worst = maxval(abs(oz(1,:) - x_mynn)/aerosol)
      write(*,'(a,a,es10.3)') regime, ': worst relative difference from MYNN''s scalar mixing = ', worst
      call check(worst <= 1.0e-12_RKIND, regime//': exch_h reproduces MYNN''s own scalar mixing')
      call check(any(abs(rnwfa) > 0.0_RKIND), regime//': MYNN mixes the aerosol profile')

      ! B. gases with structure, and a zero gas
      q0(1,:) = 2.0e-9_RKIND*exp(-col%zc/300.0_RKIND) + 1.0e-12_RKIND
      q0(2,:) = merge(1.0e-9_RKIND, 1.0e-11_RKIND, col%zc < 1000.0_RKIND)
      q0(3,:) = ozone
      q0(4,:) = 0.0_RKIND
      q = q0
      call pbl_diffuse_column(delt, m_d, G, q)
      resid = 0.0_RKIND
      do n = 1, 3
         resid = max(resid, relative_mass_residual(m_d, q0(n,:), q(n,:)))
      end do
      call check(resid <= 1.0e-12_RKIND, regime//': dry mass of every gas is conserved')
      call check(all(q(4,:) == 0.0_RKIND), regime//': a zero gas stays zero')
      call check(q(1,1) < q0(1,1) .and. q(1,3) > q0(1,3), regime//': the surface-confined gas is mixed upward')
      u(1,:) = 3.0e-9_RKIND
      call pbl_diffuse_column(delt, m_d, G, u)
      call check(all(abs(u(1,:) - 3.0e-9_RKIND) <= 1.0e-20_RKIND + 1.0e-12_RKIND*3.0e-9_RKIND), &
                 regime//': a uniform gas is unchanged')
   end subroutine local_case

end program test_mynn_gas_operator
