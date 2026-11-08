! The MPAS MYNN wrapper (physics_wrf/module_bl_mynn.F) hands MMM-physics' per-plume updraft export out as
! edmf_aw_up and edmf_ent_up (P4b of docs/plans/pbl-transport.md). One moist convective column through the
! wrapper with the updraft scheme and its diagnostics on (mynn_column.F90), kme = kte+1:
!   - the export is MYNN's own mass flux: summed over plumes, rhoz(k)*edmf_aw_up(k) equals the product of
!     the wrapper's diagnostics edmf_a(k)*edmf_w(k) (= s_aw(k+1)), with rhoz MYNN's interpolation of the
!     moist density to the top face of layer k; the mixing fractions lie in [0,1); level kme is not written;
!   - with the export absent, and with unallocated allocatable arrays passed for it (the PBL driver when gas
!     transport is off), every other output is bitwise unchanged.
program test_mynn_wrapper_export

   use mpas_kind_types, only: RKIND, I8KIND
   use unit_assert, only: check, finish
   use mynn_column, only: column_state, column_output, make_column, run_mynn, nz, kme
   use bl_mynn, only: mynn_edmf_nup

   implicit none

   real(RKIND), parameter :: delt = 30.0_RKIND, sentinel = -7.0_RKIND
   type(column_state) :: col
   type(column_output) :: out(3)
   real(RKIND), dimension(:,:,:,:), allocatable :: aw, ent, aw_none, ent_none
   real(RKIND) :: rhoz, flux, worst
   logical :: all_ok
   integer :: k, c

   call make_column('convective', .true., col)
   allocate(aw(1,kme,1,mynn_edmf_nup), ent(1,kme,1,mynn_edmf_nup))
   aw = sentinel
   ent = sentinel
   call run_mynn(col, 1, 0, 1, delt, out(1), aw_up=aw, ent_up=ent)
   call run_mynn(col, 1, 0, 1, delt, out(2))
   call run_mynn(col, 1, 0, 1, delt, out(3), aw_up=aw_none, ent_up=ent_none)
   call check(all(out(:)%errflg == 0), 'wrapper calls succeed')
   ! MYNN reports the largest mass flux as maxmf, negated when the plumes stay dry (no condensation)
   call check(abs(out(1)%maxmf) > 0.0_RKIND .and. any(aw(1,1:nz,1,:) > 0.0_RKIND), 'the updraft scheme transports')

   worst = 0.0_RKIND
   all_ok = .true.
   do k = 1, nz-1
      rhoz = (col%rho(k)*col%dz(k+1) + col%rho(k+1)*col%dz(k))/(col%dz(k+1) + col%dz(k))
      flux = rhoz*sum(aw(1,k,1,:))
      if (out(1)%edmf_a(k) > 0.0_RKIND) then
         worst = max(worst, abs(flux - out(1)%edmf_a(k)*out(1)%edmf_w(k))/(out(1)%edmf_a(k)*out(1)%edmf_w(k)))
      else
         all_ok = all_ok .and. flux == 0.0_RKIND
      end if
   end do
   write(*,'(a,es10.3)') 'worst relative difference of sum(rhoz*aw) from edmf_a*edmf_w = ', worst
   call check(worst <= 1.0e-12_RKIND .and. all_ok, &
              'summed over plumes, rhoz*edmf_aw_up(k) is MYNN''s updraft mass flux on the top face of layer k')
   call check(all(aw(1,1:nz,1,:) >= 0.0_RKIND) .and. all(aw(1,nz,1,:) == 0.0_RKIND), &
              'the export is nonnegative and zero at the model top')
   call check(all(ent(1,1:nz,1,:) >= 0.0_RKIND .and. ent(1,1:nz,1,:) < 1.0_RKIND) .and. &
              any(ent(1,1:nz,1,:) > 0.0_RKIND), 'mixing fractions lie in [0,1)')
   call check(all(ent(1,1:nz,1,:) == 0.0_RKIND .or. aw(1,1:nz,1,:) > 0.0_RKIND), &
              'a mixing fraction is set only where the plume leaves the layer')
   call check(all(aw(1,kme,1,:) == sentinel) .and. all(ent(1,kme,1,:) == sentinel), 'level kme is not written')
   do c = 2, 3
      call check(same(out(c), out(1)), merge('export absent: every other output is unchanged (bitwise)       ', &
                                             'export unallocated: every other output is unchanged (bitwise)  ', c == 2))
   end do

   call finish('test_mynn_wrapper_export')

contains

   logical function same(a, b)
      type(column_output), intent(in) :: a, b
      same = bits(a%exch_h, b%exch_h) .and. bits(a%exch_m, b%exch_m) .and. bits(a%ru, b%ru) .and. &
             bits(a%rv, b%rv) .and. bits(a%rth, b%rth) .and. bits(a%rqv, b%rqv) .and. bits(a%rqc, b%rqc) .and. &
             bits(a%qke, b%qke) .and. bits(a%edmf_a, b%edmf_a) .and. bits(a%edmf_w, b%edmf_w) .and. &
             bits([a%pblh, a%maxmf], [b%pblh, b%maxmf]) .and. a%kpbl == b%kpbl
   end function same

   logical function bits(a, b)
      real(RKIND), intent(in) :: a(:), b(:)
      bits = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
   end function bits

end program test_mynn_wrapper_export
