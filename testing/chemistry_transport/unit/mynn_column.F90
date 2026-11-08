! Synthetic single-column driver for the MYNN unit tests (P4 of docs/plans/pbl-transport.md): one
! column through the MPAS wrapper mynn_bl_driver (physics_wrf/module_bl_mynn.F), called the way the PBL
! driver calls it, with the MPAS default MYNN options except those a test sets. The column is 40 layers
! (50 m thick at the surface, 10 m thicker per layer) in hydrostatic balance.
module mynn_column

   use mpas_kind_types, only: RKIND
   use bl_mynn, only: bl_mynn_init
   use module_bl_mynn, only: mynn_bl_driver

   implicit none
   private
   public :: column_state, column_output, make_column, run_mynn, nz, kme

   integer, parameter :: nz = 40, kme = nz + 1
   real(RKIND), parameter :: g = 9.81_RKIND, rd = 287.0_RKIND, rv = 461.6_RKIND, cp = 1004.5_RKIND
   real(RKIND), parameter :: p0 = 1.0e5_RKIND

   type column_state
      real(RKIND), dimension(nz) :: dz, zc, th, tt, p, exner, rho, qv, qc, u, v, rthraten
      real(RKIND) :: ps, ts, qsfc, ust, ch, hfx, qfx, rmol, wspd, znt
   end type column_state

   type column_output
      real(RKIND), dimension(nz) :: exch_h, rqoz, rnwfa, edmf_a, edmf_w
      real(RKIND) :: pblh, maxmf
      integer :: kpbl, errflg
   end type column_output

contains

   subroutine make_column(regime, moist, col)
      ! regime: 'stable' (theta rising 6 K/km, downward heat flux) or 'convective' (1 km mixed layer
      ! under 4 K/km, 250 W m-2 upward heat flux). moist = .false. gives qv = 0, so moist and dry
      ! density coincide.
      character(len=*), intent(in) :: regime
      logical, intent(in) :: moist
      type(column_state), intent(out) :: col
      real(RKIND) :: pii(nz+1), kap, thv, z
      integer :: k

      kap = rd/cp
      z = 0.0_RKIND
      do k = 1, nz
         col%dz(k) = 50.0_RKIND + 10.0_RKIND*(k-1)
         col%zc(k) = z + 0.5_RKIND*col%dz(k)
         z = z + col%dz(k)
      end do
      select case (regime)
      case ('stable')
         col%th = 300.0_RKIND + 6.0e-3_RKIND*col%zc
      case ('convective')
         col%th = merge(300.0_RKIND, 300.0_RKIND + 4.0e-3_RKIND*(col%zc - 1000.0_RKIND), col%zc < 1000.0_RKIND)
      case default
         error stop 'make_column: unknown regime'
      end select
      if (moist) then
         col%qv = 8.0e-3_RKIND*exp(-col%zc/2500.0_RKIND)
      else
         col%qv = 0.0_RKIND
      end if
      col%qc = 0.0_RKIND
      col%rthraten = 0.0_RKIND
      col%u = 5.0_RKIND + 5.0e-3_RKIND*col%zc
      col%v = 2.0_RKIND
      col%ps = p0
      pii(1) = 1.0_RKIND
      do k = 1, nz
         thv = col%th(k)*(1.0_RKIND + 0.608_RKIND*col%qv(k))
         pii(k+1) = pii(k) - g*col%dz(k)/(cp*thv)
         col%exner(k) = 0.5_RKIND*(pii(k) + pii(k+1))
         col%p(k) = p0*col%exner(k)**(1.0_RKIND/kap)
         col%tt(k) = col%th(k)*col%exner(k)
         col%rho(k) = col%p(k)/(rd*col%tt(k)*(1.0_RKIND + 0.608_RKIND*col%qv(k)))
      end do
      col%znt = 0.1_RKIND
      col%wspd = sqrt(col%u(1)**2 + col%v(1)**2)
      col%qsfc = col%qv(1)
      col%ch = 0.01_RKIND
      if (regime == 'stable') then
         col%hfx = -30.0_RKIND
         col%ust = 0.2_RKIND
         col%rmol = 0.04_RKIND
         col%ts = col%tt(1) - 3.0_RKIND
      else
         col%hfx = 250.0_RKIND
         col%ust = 0.4_RKIND
         col%rmol = -0.04_RKIND
         col%ts = col%tt(1) + 5.0_RKIND
      end if
      col%qfx = merge(1.0e-4_RKIND, 0.0_RKIND, moist)
   end subroutine make_column

   subroutine run_mynn(col, edmf, edmf_dd, output, delt, out, qoz, rqoz, nwfa, rnwfa)
      ! One MYNN step (initflag = 1: MYNN initializes its turbulence from the column first). nwfa, when
      ! given, switches on MYNN's mixing of water-friendly aerosol (MYNN reports that tendency only with
      ! ice-friendly aerosol on as well, so nifa = nwfa). The pinned MYNN mixes ozone but reports a zero
      ! ozone tendency, so qoz cannot serve as a reference.
      type(column_state), intent(in) :: col
      integer, intent(in) :: edmf, edmf_dd, output
      real(RKIND), intent(in) :: delt
      type(column_output), intent(out) :: out
      real(RKIND), intent(in), optional :: qoz(nz), nwfa(nz)
      real(RKIND), intent(out), optional :: rqoz(nz), rnwfa(nz)

      real(RKIND), dimension(1,kme,1) :: dz, u, w, v, th, tt, p, exner, rho, qv, rthraten, qc, qi, qs, qoz3, &
                                         nc, ni, nifa, nwfa3, nbca, cldfra_bl, qc_bl, qi_bl, el_pbl, qke, &
                                         qke_adv, cov, qsq, tsq, sh3d, sm3d, rublten, rvblten, rthblten, &
                                         rqvblten, rqcblten, rqiblten, rqsblten, rqozblten, rncblten, &
                                         rniblten, rnifablten, rnwfablten, rnbcablten, edmf_a, edmf_w, &
                                         edmf_qt, edmf_thl, edmf_ent, edmf_qc, sub_thl, sub_sqv, det_thl, &
                                         det_sqv, exch_h, exch_m, dqke, qwt, qshear, qbuoy, qdiss
      real(RKIND), dimension(1,1) :: dx, xland, ps, ts, qsfc, ust, ch, hfx, qfx, rmol, wspd, uoce, voce, znt, &
                                     pblh, maxwidth, maxmf, ztop_plume
      integer, dimension(1,1) :: kpbl, ktop_plume
      character(len=256) :: errmsg
      integer :: errflg

      call bl_mynn_init(cp, 1846.4_RKIND, 2106.0_RKIND, 4190.0_RKIND, rv/rd - 1.0_RKIND, rd/rv, g, 0.4_RKIND, p0, &
                        rd, rv, 0.6112_RKIND, 17.67_RKIND, 29.65_RKIND, 273.15_RKIND, 3.50e5_RKIND, &
                        2.85e6_RKIND, 2.50e6_RKIND, errmsg, errflg)
      if (errflg /= 0) error stop 'bl_mynn_init failed'

      call fill(dz, col%dz)
      call fill(u, col%u)
      call fill(v, col%v)
      w = 0.0_RKIND
      call fill(th, col%th)
      call fill(tt, col%tt)
      call fill(p, col%p)
      call fill(exner, col%exner)
      call fill(rho, col%rho)
      call fill(qv, col%qv)
      call fill(qc, col%qc)
      call fill(rthraten, col%rthraten)
      qi = 0.0_RKIND
      qs = 0.0_RKIND
      nc = 0.0_RKIND
      ni = 0.0_RKIND
      nbca = 0.0_RKIND
      qoz3 = 0.0_RKIND
      nwfa3 = 0.0_RKIND
      if (present(qoz)) call fill(qoz3, qoz)
      if (present(nwfa)) call fill(nwfa3, nwfa)
      nifa = nwfa3
      cldfra_bl = 0.0_RKIND
      qc_bl = 0.0_RKIND
      qi_bl = 0.0_RKIND
      el_pbl = 0.0_RKIND
      qke = 0.0_RKIND
      qke_adv = 0.0_RKIND
      cov = 0.0_RKIND
      qsq = 0.0_RKIND
      tsq = 0.0_RKIND
      sh3d = 0.0_RKIND
      sm3d = 0.0_RKIND
      rublten = 0.0_RKIND
      rvblten = 0.0_RKIND
      rthblten = 0.0_RKIND
      rqvblten = 0.0_RKIND
      rqcblten = 0.0_RKIND
      rqiblten = 0.0_RKIND
      rqsblten = 0.0_RKIND
      rqozblten = 0.0_RKIND
      rncblten = 0.0_RKIND
      rniblten = 0.0_RKIND
      rnifablten = 0.0_RKIND
      rnwfablten = 0.0_RKIND
      rnbcablten = 0.0_RKIND
      edmf_a = 0.0_RKIND
      edmf_w = 0.0_RKIND
      edmf_qt = 0.0_RKIND
      edmf_thl = 0.0_RKIND
      edmf_ent = 0.0_RKIND
      edmf_qc = 0.0_RKIND
      sub_thl = 0.0_RKIND
      sub_sqv = 0.0_RKIND
      det_thl = 0.0_RKIND
      det_sqv = 0.0_RKIND
      dx = 12000.0_RKIND
      xland = 1.0_RKIND
      ps = col%ps
      ts = col%ts
      qsfc = col%qsfc
      ust = col%ust
      ch = col%ch
      hfx = col%hfx
      qfx = col%qfx
      rmol = col%rmol
      wspd = col%wspd
      uoce = 0.0_RKIND
      voce = 0.0_RKIND
      znt = col%znt
      pblh = 0.0_RKIND
      kpbl = 1
      ktop_plume = 0

      call mynn_bl_driver( &
           f_qc=.true., f_qi=.false., f_qs=.false., f_qoz=present(qoz), f_nc=.false., f_ni=.false., &
           f_nifa=present(nwfa), f_nwfa=present(nwfa), f_nbca=.false., icloud_bl=1, delt=delt, dx=dx, &
           xland=xland, ps=ps, ts=ts, qsfc=qsfc, ust=ust, ch=ch, hfx=hfx, qfx=qfx, rmol=rmol, wspd=wspd, &
           znt=znt, uoce=uoce, voce=voce, dz=dz, u=u, v=v, w=w, th=th, tt=tt, p=p, exner=exner, rho=rho, &
           qv=qv, qc=qc, qi=qi, qs=qs, nc=nc, ni=ni, nifa=nifa, nwfa=nwfa3, nbca=nbca, qoz=qoz3, &
           rthraten=rthraten, pblh=pblh, kpbl=kpbl, cldfra_bl=cldfra_bl, qc_bl=qc_bl, qi_bl=qi_bl, &
           maxwidth=maxwidth, maxmf=maxmf, ktop_plume=ktop_plume, ztop_plume=ztop_plume, dqke=dqke, &
           qke_adv=qke_adv, tsq=tsq, qsq=qsq, cov=cov, el_pbl=el_pbl, rublten=rublten, rvblten=rvblten, &
           rthblten=rthblten, rqvblten=rqvblten, rqcblten=rqcblten, rqiblten=rqiblten, rqsblten=rqsblten, &
           rncblten=rncblten, rniblten=rniblten, rnifablten=rnifablten, rnwfablten=rnwfablten, &
           rnbcablten=rnbcablten, rqozblten=rqozblten, edmf_a=edmf_a, edmf_w=edmf_w, edmf_qt=edmf_qt, &
           edmf_thl=edmf_thl, edmf_ent=edmf_ent, edmf_qc=edmf_qc, sub_thl=sub_thl, sub_sqv=sub_sqv, &
           det_thl=det_thl, det_sqv=det_sqv, exch_h=exch_h, exch_m=exch_m, qke=qke, qwt=qwt, &
           qshear=qshear, qbuoy=qbuoy, qdiss=qdiss, sh3d=sh3d, sm3d=sm3d, spp_pbl=0, &
           do_restart=.false., do_DAcycling=.false., initflag=1, bl_mynn_tkeadvect=.false., &
           bl_mynn_tkebudget=0, bl_mynn_cloudpdf=2, bl_mynn_mixlength=2, bl_mynn_closure=2.5_RKIND, &
           bl_mynn_stfunc=1, bl_mynn_topdown=0, bl_mynn_scaleaware=1, bl_mynn_dheat_opt=1, &
           bl_mynn_edmf=edmf, bl_mynn_edmf_dd=edmf_dd, bl_mynn_edmf_mom=0, bl_mynn_edmf_tke=0, &
           bl_mynn_output=output, bl_mynn_mixscalars=1, bl_mynn_cloudmix=1, bl_mynn_mixqt=0, &
           ids=1, ide=2, jds=1, jde=2, kds=1, kde=kme, ims=1, ime=1, jms=1, jme=1, kms=1, kme=kme, &
           its=1, ite=1, jts=1, jte=1, kts=1, kte=nz, errmsg=errmsg, errflg=errflg)

      out%errflg = errflg
      out%exch_h = exch_h(1,1:nz,1)
      out%rqoz = rqozblten(1,1:nz,1)
      out%rnwfa = rnwfablten(1,1:nz,1)
      out%edmf_a = edmf_a(1,1:nz,1)
      out%edmf_w = edmf_w(1,1:nz,1)
      out%pblh = pblh(1,1)
      out%maxmf = maxmf(1,1)
      out%kpbl = kpbl(1,1)
      if (present(rqoz)) rqoz = out%rqoz
      if (present(rnwfa)) rnwfa = out%rnwfa
   end subroutine run_mynn

   subroutine fill(a, column)
      real(RKIND), intent(out) :: a(1,kme,1)
      real(RKIND), intent(in) :: column(nz)
      a(1,1:nz,1) = column
      a(1,kme,1) = column(nz)
   end subroutine fill

end module mynn_column
