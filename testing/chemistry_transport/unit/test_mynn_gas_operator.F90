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
! P4b, updraft EDMF (config_mynn_edmf = 1), moist convective column:
!   C. Local diffusion followed by the plumes rebuilt from MYNN's export conserves every gas, keeps every
!      gas within its initial range, a zero gas zero and a uniform gas unchanged.
!   D. The plumes lift a surface-rich gas beyond what local diffusion does, as MYNN's own mixing of
!      water-friendly aerosol (local diffusion and updraft mass flux) does; the ratio of the two lifts is
!      reported (the discretizations differ: MYNN uses centered implicit subsidence with a raised
!      diffusivity, the gas operator upwind explicit substeps).
! P4c, downdraft EDMF (config_mynn_edmf_dd = 1), stratocumulus-topped column with strong cloud-top cooling:
!   E. With MYNN's updrafts and downdrafts the operator conserves and bounds every gas, keeps a zero gas zero
!      and a uniform gas unchanged, and the downdrafts carry a gas confined to the cloud-top layer (where
!      they start) downward.
program test_mynn_gas_operator

   use mpas_kind_types, only: RKIND
   use unit_assert, only: check, finish
   use mynn_column, only: column_state, column_output, make_column, run_mynn, cloud_top_layer, nz, kme
   use bl_mynn, only: mynn_edmf_nup, mynn_edmf_ndown
   use mpas_atmphys_pbl_edmf, only: pbl_edmf_updraft_fluxes, pbl_edmf_downdraft_fluxes, pbl_edmf_column
   use mpas_atmphys_pbl_gas, only: pbl_gas_exchange_faces
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   real(RKIND), parameter :: delt = 30.0_RKIND
   real(RKIND) :: gsum_stable, gsum_convective

   call local_case('stable', gsum_stable)
   call local_case('convective', gsum_convective)
   call check(gsum_convective > gsum_stable, 'convective MYNN exchange exceeds the stable exchange')
   call edmf_case()
   call downdraft_case()

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

   subroutine edmf_case()
      type(column_state) :: col
      type(column_output) :: out
      real(RKIND), dimension(:,:,:,:), allocatable :: aw, ent
      real(RKIND) :: aerosol(nz), rnwfa(nz), rho_d(nz), m_d(nz), G(nz+1), Mu(nz+1,mynn_edmf_nup)
      real(RKIND) :: Eu(nz,mynn_edmf_nup), q0(5,nz), q(5,nz), ql(5,nz), resid, tol, lift_ours, lift_mynn
      real(RKIND) :: x_mynn(nz), xl(1,nz), none_n(nz+1,0), none_e(nz,0)
      logical :: ok, bounded
      integer :: n, nsub

      call make_column('convective', .true., col)
      aerosol = 1.0e9_RKIND*(1.0_RKIND + max(0.0_RKIND, 1.0_RKIND - col%zc/3000.0_RKIND))
      allocate(aw(1,kme,1,mynn_edmf_nup), ent(1,kme,1,mynn_edmf_nup))
      call run_mynn(col, 1, 0, 0, delt, out, nwfa=aerosol, rnwfa=rnwfa, aw_up=aw, ent_up=ent)
      call check(out%errflg == 0 .and. any(aw(1,1:nz,1,:) > 0.0_RKIND), 'EDMF: MYNN succeeds and its plumes transport')
      rho_d = col%rho/(1.0_RKIND + col%qv)
      m_d = rho_d*col%dz
      call pbl_gas_exchange_faces(out%exch_h, rho_d, col%dz, G, ok)
      call pbl_edmf_updraft_fluxes(aw(1,1:nz,1,:), ent(1,1:nz,1,:), rho_d, col%dz, Mu, Eu, ok)
      call check(ok, 'EDMF: the exported plumes build a valid operator')

      ! C. gases with structure, a uniform and a zero gas
      q0(1,:) = 2.0e-9_RKIND*exp(-col%zc/300.0_RKIND) + 1.0e-12_RKIND
      q0(2,:) = merge(1.0e-9_RKIND, 1.0e-11_RKIND, col%zc < 1000.0_RKIND)
      q0(3,:) = 5.0e-8_RKIND*(0.5_RKIND + 0.5_RKIND*min(col%zc/3000.0_RKIND, 1.0_RKIND))
      q0(4,:) = 3.0e-9_RKIND
      q0(5,:) = 0.0_RKIND
      ql = q0
      call pbl_diffuse_column(delt, m_d, G, ql)
      q = ql
      call pbl_edmf_column(delt, m_d, Mu, Eu, none_n, none_e, q, nsub, ok)
      write(*,'(a,i0)') 'EDMF: substeps at the 30 s step = ', nsub
      call check(ok, 'EDMF: the column step is accepted')
      resid = 0.0_RKIND
      bounded = .true.
      do n = 1, 3
         resid = max(resid, relative_mass_residual(m_d, q0(n,:), q(n,:)))
         tol = 1.0e-20_RKIND + 1.0e-12_RKIND*maxval(q0(n,:))
         bounded = bounded .and. all(q(n,:) >= minval(q0(n,:)) - tol) .and. all(q(n,:) <= maxval(q0(n,:)) + tol)
      end do
      call check(resid <= 1.0e-12_RKIND, 'EDMF: dry mass of every gas is conserved')
      call check(bounded, 'EDMF: every gas stays within its initial range')
      call check(all(abs(q(4,:) - 3.0e-9_RKIND) <= 1.0e-20_RKIND + 1.0e-12_RKIND*3.0e-9_RKIND), &
                 'EDMF: a uniform gas is unchanged')
      call check(all(q(5,:) == 0.0_RKIND), 'EDMF: a zero gas stays zero')

      ! D. lift of a surface-rich profile: plumes on top of local diffusion, ours and MYNN's
      lift_ours = sum(m_d*col%zc*q(1,:)) - sum(m_d*col%zc*ql(1,:))
      xl(1,:) = aerosol
      call pbl_diffuse_column(delt, m_d, G, xl)
      x_mynn = aerosol + delt*rnwfa
      lift_mynn = sum(m_d*col%zc*x_mynn) - sum(m_d*col%zc*xl(1,:))
      call pbl_edmf_column(delt, m_d, Mu, Eu, none_n, none_e, xl, nsub, ok)
      lift_ours = sum(m_d*col%zc*xl(1,:)) - (sum(m_d*col%zc*x_mynn) - lift_mynn)
      write(*,'(a,es10.3)') 'EDMF: lift of the aerosol beyond local diffusion, gas operator / MYNN = ', &
                            lift_ours/lift_mynn
      call check(lift_ours > 0.0_RKIND .and. lift_mynn > 0.0_RKIND, &
                 'EDMF: the plumes lift a surface-rich profile beyond local diffusion, as MYNN''s own mixing does')
   end subroutine edmf_case

   subroutine downdraft_case()
      type(column_state) :: col
      type(column_output) :: out
      real(RKIND), dimension(:,:,:,:), allocatable :: aw, ent, awd, entd
      real(RKIND) :: rho_d(nz), m_d(nz), G(nz+1), Mu(nz+1,mynn_edmf_nup), Eu(nz,mynn_edmf_nup)
      real(RKIND) :: Nd(nz+1,mynn_edmf_ndown), Ed(nz,mynn_edmf_ndown), q0(4,nz), q(4,nz), qu(1,nz), qd(1,nz)
      real(RKIND) :: resid, tol
      logical :: ok, okd, bounded
      integer :: n, nsub, ktop

      call make_column('cloudtop', .true., col)
      allocate(aw(1,kme,1,mynn_edmf_nup), ent(1,kme,1,mynn_edmf_nup))
      allocate(awd(1,kme,1,mynn_edmf_ndown), entd(1,kme,1,mynn_edmf_ndown))
      call run_mynn(col, 1, 1, 0, delt, out, aw_up=aw, ent_up=ent, aw_dn=awd, ent_dn=entd)
      call check(out%errflg == 0 .and. any(awd(1,1:nz,1,:) > 0.0_RKIND), &
                 'downdrafts: MYNN succeeds and its downdrafts transport')
      rho_d = col%rho/(1.0_RKIND + col%qv)
      m_d = rho_d*col%dz
      call pbl_gas_exchange_faces(out%exch_h, rho_d, col%dz, G, ok)
      call pbl_edmf_updraft_fluxes(aw(1,1:nz,1,:), ent(1,1:nz,1,:), rho_d, col%dz, Mu, Eu, ok)
      call pbl_edmf_downdraft_fluxes(awd(1,1:nz,1,:), entd(1,1:nz,1,:), rho_d, col%dz, Nd, Ed, okd)
      call check(ok .and. okd, 'downdrafts: the exported plumes build a valid operator')

      q0(1,:) = 2.0e-9_RKIND*exp(-col%zc/300.0_RKIND) + 1.0e-12_RKIND
      ktop = cloud_top_layer(col)
      q0(2,:) = 1.0e-11_RKIND
      q0(2,ktop) = 1.0e-9_RKIND                                          ! confined to the cloud-top layer
      q0(3,:) = 3.0e-9_RKIND
      q0(4,:) = 0.0_RKIND
      q = q0
      call pbl_diffuse_column(delt, m_d, G, q)
      call pbl_edmf_column(delt, m_d, Mu, Eu, Nd, Ed, q, nsub, ok)
      write(*,'(a,i0)') 'downdrafts: substeps at the 30 s step = ', nsub
      resid = 0.0_RKIND
      bounded = ok
      do n = 1, 2
         resid = max(resid, relative_mass_residual(m_d, q0(n,:), q(n,:)))
         tol = 1.0e-20_RKIND + 1.0e-12_RKIND*maxval(q0(n,:))
         bounded = bounded .and. all(q(n,:) >= minval(q0(n,:)) - tol) .and. all(q(n,:) <= maxval(q0(n,:)) + tol)
      end do
      call check(resid <= 1.0e-12_RKIND .and. bounded, 'downdrafts: every gas is conserved and within its range')
      call check(all(abs(q(3,:) - 3.0e-9_RKIND) <= 1.0e-20_RKIND + 1.0e-12_RKIND*3.0e-9_RKIND) .and. &
                 all(q(4,:) == 0.0_RKIND), 'downdrafts: a uniform gas is unchanged and a zero gas stays zero')
      qu(1,:) = q0(2,:)
      qd(1,:) = q0(2,:)
      call pbl_edmf_column(delt, m_d, Mu, Eu, Nd(:,1:0), Ed(:,1:0), qu, nsub, ok)
      call pbl_edmf_column(delt, m_d, Mu, Eu, Nd, Ed, qd, nsub, ok)
      call check(sum(m_d*col%zc*qd(1,:)) < sum(m_d*col%zc*qu(1,:)), &
                 'downdrafts carry a gas confined to the cloud-top layer downward')
   end subroutine downdraft_case

end program test_mynn_gas_operator
