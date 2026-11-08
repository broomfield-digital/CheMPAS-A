! Column reference test for the YSU gas operator (P3 of docs/plans/pbl-transport.md).
!
! Synthetic single columns (stable, convective, and a cloud-topped column with radiative cooling and
! YSU top-down mixing) are run through the real bl_ysu_run with two passive gases (qmix) and the
! patched xkzqx export.
!   A. YSU's own tracer update, rebuilt from xkzqx with YSU's moist-pressure weights and its 2*dt
!      implicit solve, reproduces qmixtnp to roundoff: the export is the operator YSU applies to
!      tracers, and K(k) sits on the interface between layers k and k+1.
!   B. The CheMPAS-A operator (dry layer mass, model dt, pbl_gas_exchange + pbl_diffuse_column)
!      conserves dry tracer mass to roundoff, and converges to YSU's tracer tendency as dt -> 0;
!      the remaining difference is the moist/dry weighting, of order qv.
!   C. At the model dt (30 s) the two time discretizations differ: per mode, YSU's half-applied
!      2*dt step multiplies by (1+lam)/(1+2*lam), which never mixes more than half per step, and
!      backward Euler by 1/(1+lam). Against a near-exact solution (1000 backward-Euler substeps),
!      the CheMPAS-A operator is the more accurate one wherever the column is stiff (lam >= 0.5).
program test_ysu_gas_operator

   use ccpp_kind_types, only: kind_phys
   use mpas_kind_types, only: RKIND
   use unit_assert, only: check, finish
   use bl_ysu, only: bl_ysu_run
   use mpas_atmphys_pbl_gas, only: pbl_gas_exchange
   use mpas_atmphys_pbl_diffusion, only: pbl_diffuse_column
   use mpas_atmphys_tracer_column, only: relative_mass_residual

   implicit none

   integer, parameter :: nz = 40, nmix = 2
   real(kind_phys), parameter :: g = 9.81, rd = 287.0, rv = 461.6, cp = 3.5*287.0, xlv = 2.5e6
   real(kind_phys), parameter :: karman = 0.4, p0 = 1.0e5
   real(kind_phys) :: hpbl_stable, hpbl_convective
   real(kind_phys) :: kzq_plain(nz), kzq_topdown(nz)

   call run_case('stable', theta_mixed=0.0, hfx=-30.0, qfx=0.0, ust=0.2, br=0.3, psim=2.0, psih=2.0, &
                 cloudy=.false., topdown=.false., hpbl_out=hpbl_stable)
   call run_case('convective', theta_mixed=1200.0, hfx=250.0, qfx=1.0e-4, ust=0.4, br=-1.0, psim=-1.0, &
                 psih=-1.5, cloudy=.false., topdown=.false., hpbl_out=hpbl_convective)
   call check(hpbl_convective > hpbl_stable, 'convective PBL is deeper than the stable one')

   call run_case('cloud-topped, top-down off', theta_mixed=900.0, hfx=20.0, qfx=3.0e-5, ust=0.3, br=-0.1, &
                 psim=-0.2, psih=-0.3, cloudy=.true., topdown=.false., kzq_out=kzq_plain)
   call run_case('cloud-topped, top-down on', theta_mixed=900.0, hfx=20.0, qfx=3.0e-5, ust=0.3, br=-0.1, &
                 psim=-0.2, psih=-0.3, cloudy=.true., topdown=.true., kzq_out=kzq_topdown)
   call check(any(abs(kzq_topdown - kzq_plain) > 1.0e-6*maxval(kzq_plain)), &
              'top-down mixing changes the exported operator')

   call finish('test_ysu_gas_operator')

contains

   subroutine run_case(name, theta_mixed, hfx, qfx, ust, br, psim, psih, cloudy, topdown, hpbl_out, kzq_out)
      character(len=*), intent(in) :: name
      real(kind_phys), intent(in) :: theta_mixed            ! depth of a well-mixed layer [m]; 0 = stable
      real(kind_phys), intent(in) :: hfx, qfx, ust, br, psim, psih
      logical, intent(in) :: cloudy, topdown
      real(kind_phys), intent(out), optional :: hpbl_out
      real(kind_phys), intent(out), optional :: kzq_out(nz)

      real(kind_phys), dimension(1,nz) :: u, v, t, qv, qc, qi, pi2d, p2d, dz, rthraten
      real(kind_phys), dimension(1,nz) :: utnp, vtnp, ttnp, qvtnp, qctnp, qitnp, exch_h, exch_m, xkzq
      real(kind_phys), dimension(1,nz+1) :: p2di
      real(kind_phys), dimension(1,nz,nmix) :: qmix, r_ysu
      real(kind_phys), dimension(1) :: psfc, znt, ust1, hpbl, xland, hfx1, qfx1, wspd, br1, psim1, psih1
      real(kind_phys), dimension(1) :: wstar, delta, u10, v10, uox, vox
      integer, dimension(1) :: kpbl
      character(len=256) :: errmsg
      integer :: errflg, k, n, ip
      real(kind_phys) :: zi(nz+1), zc(nz), theta(nz), thv, pii(nz+1), kap, del(nz), m_d(nz)
      real(kind_phys) :: r_ref(nz), r_ours(nz), gex(nz+1), q(1,nz), worst_a, worst_b, scale, resid, dts(2)
      real(kind_phys) :: diff_b(2), q_ex(1,nz), r_exact(nz), err_ours, err_ysu, lam
      integer :: s
      logical :: ok

      kap = rd/cp
      zi(1) = 0.0
      do k = 1, nz
         dz(1,k) = 50.0 + 10.0*(k-1)
         zi(k+1) = zi(k) + dz(1,k)
         zc(k) = 0.5*(zi(k) + zi(k+1))
      end do
      do k = 1, nz
         if (zc(k) < theta_mixed) then
            theta(k) = 300.0
         else
            theta(k) = 300.0 + 4.0e-3*(zc(k) - theta_mixed)
         end if
         qv(1,k) = 8.0e-3*exp(-zc(k)/2500.0)
         u(1,k) = 5.0 + 5.0e-3*zc(k)
         v(1,k) = 2.0
         qc(1,k) = 0.0
         qi(1,k) = 0.0
         rthraten(1,k) = 0.0
      end do
      if (cloudy) then
         do k = 1, nz
            if (zc(k) > 500.0 .and. zc(k) < theta_mixed + 300.0) qc(1,k) = 3.0e-4
         end do
         ! cloud-top radiative cooling inside the PBL (YSU sums it over the layers below kpbl)
         do k = 1, nz
            if (zc(k) > 700.0 .and. zc(k) < theta_mixed) rthraten(1,k) = -2.0e-4
         end do
      end if

      ! hydrostatic Exner at interfaces, then pressures and temperature
      pii(1) = 1.0
      do k = 1, nz
         thv = theta(k)*(1.0 + 0.608*qv(1,k))
         pii(k+1) = pii(k) - g*dz(1,k)/(cp*thv)
      end do
      do k = 1, nz+1
         p2di(1,k) = p0*pii(k)**(1.0/kap)
      end do
      do k = 1, nz
         pi2d(1,k) = 0.5*(pii(k) + pii(k+1))
         p2d(1,k) = p0*pi2d(1,k)**(1.0/kap)
         t(1,k) = theta(k)*pi2d(1,k)
      end do

      do k = 1, nz
         qmix(1,k,1) = merge(1.0e-9_kind_phys, 0.0_kind_phys, k >= 3 .and. k <= 5)
         qmix(1,k,2) = 2.0e-9*(1.0 - zc(k)/1.0e4)
      end do
      psfc = p0
      znt = 0.1
      ust1 = ust
      xland = 1.0
      hfx1 = hfx
      qfx1 = qfx
      wspd = sqrt(u(1,1)**2 + v(1,1)**2)
      br1 = br
      psim1 = psim
      psih1 = psih
      u10 = 0.8*u(1,1)
      v10 = 0.8*v(1,1)
      uox = 0.0
      vox = 0.0
      xkzq = -1.0

      dts = [30.0_kind_phys, 1.0e-3_kind_phys]
      err_ours = 0.0
      err_ysu = 0.0
      lam = 0.0
      do ip = 1, 2
         call bl_ysu_run(ux=u, vx=v, tx=t, qvx=qv, qcx=qc, qix=qi, nmix=nmix, qmix=qmix, p2d=p2d, &
                         p2di=p2di, pi2d=pi2d, f_qc=.true., f_qi=.false., utnp=utnp, vtnp=vtnp, ttnp=ttnp, &
                         qvtnp=qvtnp, qctnp=qctnp, qitnp=qitnp, qmixtnp=r_ysu, cp=cp, g=g, rovcp=rd/cp, rd=rd, &
                         rovg=rd/g, ep1=rv/rd-1.0, ep2=rd/rv, karman=karman, xlv=xlv, rv=rv, dz8w2d=dz, &
                         psfcpa=psfc, znt=znt, ust=ust1, hpbl=hpbl, psim=psim1, psih=psih1, xland=xland, &
                         hfx=hfx1, qfx=qfx1, wspd=wspd, br=br1, dt=dts(ip), kpbl1d=kpbl, exch_hx=exch_h, &
                         exch_mx=exch_m, wstar=wstar, delta=delta, u10=u10, v10=v10, uox=uox, vox=vox, &
                         rthraten=rthraten, ysu_topdown_pblmix=topdown, flag_bep=.false., its=1, ite=1, &
                         kte=nz, kme=nz+1, xkzqx=xkzq, errmsg=errmsg, errflg=errflg)
         if (ip == 1) then
            call check(errflg == 0, trim(name)//': bl_ysu_run succeeds')
            call check(all(xkzq(1,1:nz-1) > 0.0) .and. xkzq(1,nz) == 0.0, &
                       trim(name)//': export is positive at interior interfaces and zero at the top')
            if (present(hpbl_out)) hpbl_out = hpbl(1)
            if (present(kzq_out)) kzq_out = xkzq(1,:)
         end if

         ! A: rebuild YSU's tracer update from the export (moist weights, 2*dt, ad(1) = 1)
         do k = 1, nz
            del(k) = p2di(1,k) - p2di(1,k+1)
         end do
         worst_a = 0.0
         worst_b = 0.0
         do n = 1, nmix
            call ysu_tracer_reference(dts(ip), del, p2d(1,:), zi, xkzq(1,:), qmix(1,:,n), r_ref)
            scale = maxval(abs(r_ysu(1,:,n)))
            worst_a = max(worst_a, maxval(abs(r_ref - r_ysu(1,:,n)))/scale)

            ! B: CheMPAS-A operator on dry layer mass with the model dt
            m_d = del/g/(1.0 + qv(1,:))
            call pbl_gas_exchange(xkzq(1,:), m_d, dz(1,:), gex, ok)
            q(1,:) = qmix(1,:,n)
            call pbl_diffuse_column(dts(ip), m_d, gex, q)
            r_ours = (q(1,:) - qmix(1,:,n))/dts(ip)
            resid = relative_mass_residual(m_d, qmix(1,:,n), q(1,:))
            worst_b = max(worst_b, maxval(abs(r_ours - r_ysu(1,:,n)))/scale)
            if (ip == 1) then
               call check(ok .and. resid <= 1.0e-12, trim(name)//': CheMPAS-A operator conserves dry mass')
               ! C: near-exact solution of the same dry operator
               q_ex(1,:) = qmix(1,:,n)
               do s = 1, 1000
                  call pbl_diffuse_column(dts(ip)/1000.0, m_d, gex, q_ex)
               end do
               r_exact = (q_ex(1,:) - qmix(1,:,n))/dts(ip)
               err_ours = max(err_ours, maxval(abs(r_ours - r_exact))/scale)
               err_ysu = max(err_ysu, maxval(abs(r_ysu(1,:,n) - r_exact))/scale)
               do k = 1, nz
                  lam = max(lam, dts(ip)*(gex(k) + gex(k+1))/m_d(k))
               end do
            end if
         end do
         diff_b(ip) = worst_b
         if (ip == 1) call check(worst_a <= 1.0e-10, trim(name)//': export reproduces YSU''s tracer update')
      end do
      write(*, '(a,a,a,es9.2,a,es9.2,a,es9.2)') '  ', trim(name), ': hpbl ', hpbl(1), ' m; max K ', &
         maxval(xkzq), ' m2/s; max dt*G/m_d at 30 s ', lam
      write(*, '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)') '    |R_ours-R_ysu|/max|R_ysu|: dt=30 s ', diff_b(1), &
         ', dt=1e-3 s ', diff_b(2), ';  error vs near-exact at 30 s: CheMPAS-A ', err_ours, ', YSU ', err_ysu
      call check(diff_b(2) <= 2.0*maxval(qv) + 1.0e-3, &
                 trim(name)//': as dt -> 0 the operators agree to the moist/dry weighting (order qv)')
      if (lam >= 0.5) call check(err_ours <= err_ysu, &
                 trim(name)//': at dt = 30 s the CheMPAS-A step is closer to the exact solution than YSU''s')
   end subroutine run_case

   subroutine ysu_tracer_reference(dt, del, p2d, zi, xkzq, q0, r)
      ! YSU's qmix update (bl_ysu.F90, moisture matrix with sfk = vlk = 1), solved independently.
      real(kind_phys), intent(in) :: dt, del(nz), p2d(nz), zi(nz+1), xkzq(nz), q0(nz)
      real(kind_phys), intent(out) :: r(nz)
      real(kind_phys) :: al(nz), ad(nz), au(nz), za(nz), dza, dt2, dsig, dsdz2, a(nz), c(nz), w, x(nz), d(nz)
      integer :: k

      dt2 = 2.0*dt
      do k = 1, nz
         za(k) = 0.5*(zi(k) + zi(k+1))
      end do
      ad = 0.0
      al = 0.0
      au = 0.0
      ad(1) = 1.0
      do k = 1, nz-1
         dza = za(k+1) - za(k)
         dsig = p2d(k) - p2d(k+1)
         dsdz2 = dsig*xkzq(k)/dza/dza
         au(k) = -dt2/del(k)*dsdz2
         al(k) = -dt2/del(k+1)*dsdz2
         ad(k) = ad(k) - au(k)
         ad(k+1) = 1.0 - al(k)
      end do
      ! rows: al(k-1)*x(k-1) + ad(k)*x(k) + au(k)*x(k+1) = q0(k)
      a(1) = 0.0
      c(nz) = 0.0
      do k = 2, nz
         a(k) = al(k-1)
      end do
      do k = 1, nz-1
         c(k) = au(k)
      end do
      d(1) = q0(1)/ad(1)
      x(1) = c(1)/ad(1)
      do k = 2, nz
         w = ad(k) - a(k)*x(k-1)
         if (k < nz) x(k) = c(k)/w
         d(k) = (q0(k) - a(k)*d(k-1))/w
      end do
      r(nz) = d(nz)
      do k = nz-1, 1, -1
         r(k) = d(k) - x(k)*r(k+1)
      end do
      r = (r - q0)/dt2
   end subroutine ysu_tracer_reference

end program test_ysu_gas_operator
