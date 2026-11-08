! The MPAS YSU wrapper (physics_wrf/module_bl_ysu.F) hands MMM-physics' xkzqx export out as exch_q
! (P3 of docs/plans/pbl-transport.md). One convective column, run through the wrapper with the
! driver's (ims:ime, kms:kme, jms:jme) arrays, kme = kte+1:
!   - exch_q(i,k,j) equals bl_ysu_run's xkzqx(i,k), k = kts..kte, bitwise: the interface between
!     layers k and k+1, with zero at the top; level kme is not written;
!   - with exch_q absent, and with an unallocated allocatable array passed for it (the PBL driver
!     when gas transport is off), every other output is bitwise unchanged. (The ice tendency is
!     excluded: with flag_qi = .false. YSU does not set it, and the wrapper copies out undefined
!     values, which the MPAS driver does not use.)
program test_ysu_wrapper_export

   use ccpp_kind_types, only: kind_phys
   use mpas_kind_types, only: I8KIND
   use unit_assert, only: check, finish
   use bl_ysu, only: bl_ysu_run
   use module_bl_ysu, only: ysu

   implicit none

   integer, parameter :: nz = 30, kme = nz + 1
   real(kind_phys), parameter :: g = 9.81, rd = 287.0, rv = 461.6, cp = 3.5*287.0, xlv = 2.5e6
   real(kind_phys), parameter :: karman = 0.4, p0 = 1.0e5, dt = 30.0, sentinel = -7.0

   real(kind_phys), dimension(1,kme,1) :: u, v, t, qv, qc, qi, p, pi, dz, rthraten, pint
   real(kind_phys), dimension(1,kme,1,3) :: ru, rv_, rth, rqv, rqc, rqi, exh, exm
   real(kind_phys), dimension(:,:,:), allocatable :: exq, unallocated
   real(kind_phys), dimension(1,1) :: psfc, znt, ust, wspd, psim, psih, xland, hfx, qfx, br, uoce, voce, ctopo
   real(kind_phys), dimension(1,1,3) :: hpbl, wstar, delta, u10, v10
   integer, dimension(1,1,3) :: kpbl
   real(kind_phys) :: zc, theta, thv, pii(kme), kap
   ! direct bl_ysu_run reference
   real(kind_phys), dimension(1,nz) :: utnp, vtnp, ttnp, qvtnp, qctnp, qitnp, exch_h, exch_m, xkzq
   real(kind_phys), dimension(1,nz,0) :: qmix0, qmixtnp0
   real(kind_phys), dimension(1) :: hpbl1, wstar1, delta1, u101, v101
   integer, dimension(1) :: kpbl1
   character(len=256) :: errmsg
   character(len=80) :: label
   integer :: errflg, k, c

   kap = rd/cp
   u = 0.0
   v = 0.0
   t = 0.0
   qv = 0.0
   p = 0.0
   pi = 0.0
   dz = 0.0
   pii(1) = 1.0
   do k = 1, nz
      dz(1,k,1) = 50.0 + 10.0*(k-1)
      zc = sum(dz(1,1:k-1,1)) + 0.5*dz(1,k,1)
      theta = merge(300.0_kind_phys, 300.0 + 4.0e-3*(zc - 1000.0), zc < 1000.0)
      qv(1,k,1) = 8.0e-3*exp(-zc/2500.0)
      u(1,k,1) = 5.0 + 5.0e-3*zc
      v(1,k,1) = 2.0
      thv = theta*(1.0 + 0.608*qv(1,k,1))
      pii(k+1) = pii(k) - g*dz(1,k,1)/(cp*thv)
      pi(1,k,1) = 0.5*(pii(k) + pii(k+1))
      p(1,k,1) = p0*pi(1,k,1)**(1.0/kap)
      t(1,k,1) = theta*pi(1,k,1)
   end do
   do k = 1, kme
      pint(1,k,1) = p0*pii(k)**(1.0/kap)
   end do
   qc = 0.0
   qi = 0.0
   rthraten = 0.0
   psfc = p0
   znt = 0.1
   ust = 0.4
   wspd = sqrt(u(1,1,1)**2 + v(1,1,1)**2)
   psim = -1.0
   psih = -1.5
   xland = 1.0
   hfx = 250.0
   qfx = 1.0e-4
   br = -1.0
   uoce = 0.0
   voce = 0.0
   ctopo = 1.0
   u10 = 0.8*u(1,1,1)
   v10 = 0.8*v(1,1,1)

   ! three wrapper calls: exch_q present, absent, and an unallocated allocatable (the PBL driver's
   ! kzq_gas_p when gas transport is off)
   allocate(exq(1,kme,1))
   exq = sentinel
   do c = 1, 3
      select case (c)
      case (1)
         call run_wrapper(c, exq)
      case (2)
         call run_wrapper(c)
      case (3)
         call run_wrapper(c, unallocated)
      end select
   end do

   xkzq = -1.0
   u101 = 0.8*u(1,1,1)
   v101 = 0.8*v(1,1,1)
   call bl_ysu_run(ux=u(:,1:nz,1), vx=v(:,1:nz,1), tx=t(:,1:nz,1), qvx=qv(:,1:nz,1), qcx=qc(:,1:nz,1), &
                   qix=qi(:,1:nz,1), nmix=0, qmix=qmix0, p2d=p(:,1:nz,1), p2di=pint(:,:,1), pi2d=pi(:,1:nz,1), &
                   f_qc=.true., f_qi=.false., utnp=utnp, vtnp=vtnp, ttnp=ttnp, qvtnp=qvtnp, qctnp=qctnp, &
                   qitnp=qitnp, qmixtnp=qmixtnp0, cp=cp, g=g, rovcp=rd/cp, rd=rd, rovg=rd/g, ep1=rv/rd-1.0, &
                   ep2=rd/rv, karman=karman, xlv=xlv, rv=rv, dz8w2d=dz(:,1:nz,1), psfcpa=psfc(:,1), &
                   znt=znt(:,1), ust=ust(:,1), hpbl=hpbl1, psim=psim(:,1), psih=psih(:,1), xland=xland(:,1), &
                   hfx=hfx(:,1), qfx=qfx(:,1), wspd=wspd(:,1), br=br(:,1), dt=dt, kpbl1d=kpbl1, &
                   exch_hx=exch_h, exch_mx=exch_m, wstar=wstar1, delta=delta1, u10=u101, v10=v101, &
                   uox=uoce(:,1), vox=voce(:,1), rthraten=rthraten(:,1:nz,1), ysu_topdown_pblmix=.false., &
                   ctopo=ctopo(:,1), ctopo2=ctopo(:,1), flag_bep=.false., its=1, ite=1, kte=nz, kme=kme, &
                   xkzqx=xkzq, errmsg=errmsg, errflg=errflg)
   call check(errflg == 0, 'direct bl_ysu_run succeeds')

   call check(same(exq(1,1:nz,1), xkzq(1,:)), &
              'exch_q(k) is bl_ysu_run''s xkzqx(k), the interface between layers k and k+1 (bitwise)')
   call check(exq(1,nz,1) == 0.0 .and. all(exq(1,1:nz-1,1) > 0.0), &
              'exch_q is positive at interior interfaces and zero at the top')
   call check(exq(1,kme,1) == sentinel, 'exch_q is not written at kme = kte+1')
   do c = 2, 3
      if (c == 2) label = 'exch_q absent: every other output is unchanged (bitwise)'
      if (c == 3) label = 'exch_q an unallocated array: every other output is unchanged (bitwise)'
      call check(same(ru(:,:,:,c), ru(:,:,:,1)) .and. same(rv_(:,:,:,c), rv_(:,:,:,1)) .and. &
                 same(rth(:,:,:,c), rth(:,:,:,1)) .and. same(rqv(:,:,:,c), rqv(:,:,:,1)) .and. &
                 same(rqc(:,:,:,c), rqc(:,:,:,1)) .and. &
                 same(exh(:,:,:,c), exh(:,:,:,1)) .and. same(exm(:,:,:,c), exm(:,:,:,1)) .and. &
                 same(hpbl(:,:,c), hpbl(:,:,1)) .and. same(wstar(:,:,c), wstar(:,:,1)) .and. &
                 all(kpbl(:,:,c) == kpbl(:,:,1)), trim(label))
   end do

   call finish('test_ysu_wrapper_export')

contains

   subroutine run_wrapper(c, exch_q)
      ! exch_q is passed on to the wrapper's optional explicit-shape dummy as the driver does: an
      ! allocatable, absent (case 2) or not allocated (case 3)
      integer, intent(in) :: c
      real(kind_phys), dimension(:,:,:), allocatable, intent(inout), optional :: exch_q
      ru(:,:,:,c) = sentinel
      rv_(:,:,:,c) = sentinel
      rth(:,:,:,c) = sentinel
      rqv(:,:,:,c) = sentinel
      rqc(:,:,:,c) = sentinel
      rqi(:,:,:,c) = sentinel
      exh(:,:,:,c) = sentinel
      exm(:,:,:,c) = sentinel
      u10(:,:,c) = 0.8*u(1,1,1)
      v10(:,:,c) = 0.8*v(1,1,1)
      call ysu(u3d=u, v3d=v, t3d=t, qv3d=qv, qc3d=qc, qi3d=qi, p3d=p, p3di=pint, pi3d=pi, &
               rublten=ru(:,:,:,c), rvblten=rv_(:,:,:,c), rthblten=rth(:,:,:,c), rqvblten=rqv(:,:,:,c), &
               rqcblten=rqc(:,:,:,c), rqiblten=rqi(:,:,:,c), flag_qc=.true., flag_qi=.false., &
               cp=cp, g=g, rovcp=rd/cp, rd=rd, rovg=rd/g, ep1=rv/rd-1.0, ep2=rd/rv, karman=karman, xlv=xlv, &
               rv=rv, dz8w=dz, psfc=psfc, znt=znt, ust=ust, hpbl=hpbl(:,:,c), psim=psim, psih=psih, &
               xland=xland, hfx=hfx, qfx=qfx, wspd=wspd, br=br, dt=dt, kpbl2d=kpbl(:,:,c), &
               exch_h=exh(:,:,:,c), exch_m=exm(:,:,:,c), exch_q=exch_q, wstar=wstar(:,:,c), &
               delta=delta(:,:,c), u10=u10(:,:,c), v10=v10(:,:,c), uoce=uoce, voce=voce, rthraten=rthraten, &
               ysu_topdown_pblmix=0, ctopo=ctopo, ctopo2=ctopo, idiff=0, flag_bep=.false., errmsg=errmsg, &
               errflg=errflg, ids=1, ide=2, jds=1, jde=2, kds=1, kde=kme, ims=1, ime=1, jms=1, jme=1, &
               kms=1, kme=kme, its=1, ite=1, jts=1, jte=1, kts=1, kte=nz)
      call check(errflg == 0, 'wrapper call succeeds')
   end subroutine run_wrapper

   logical function same(a, b)
      real(kind_phys), intent(in) :: a(..), b(..)
      select rank (a)
      rank (1)
         select rank (b)
         rank (1)
            same = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
         end select
      rank (2)
         select rank (b)
         rank (2)
            same = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
         end select
      rank (3)
         select rank (b)
         rank (3)
            same = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
         end select
      rank default
         same = .false.
      end select
   end function same

end program test_ysu_wrapper_export
