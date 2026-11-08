! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_pbl_edmf
 use mpas_kind_types,only: RKIND

 implicit none
 private
 public:: pbl_edmf_updraft_fluxes,   &
          pbl_edmf_updraft_plume,    &
          pbl_edmf_downdraft_fluxes, &
          pbl_edmf_downdraft_plume,  &
          pbl_edmf_column


!Closed-column mass-flux (EDMF) transport of MICM gases by MYNN's updraft plumes (P4b) and downdrafts (P4c of
!docs/plans/pbl-transport.md). MPAS level k = 1 is at the ground; face k is the lower interface of layer k.
!
!MYNN exports, per plume i, the area-weighted velocity aw(k,i) = UPA*UPW*Psig_w on the top face of layer k
!and the fraction eps(k,i) = ENT*dz with which it mixes the plume through layer k. The dry plume mass flux on
!face k+1 is M(k+1,i) = rhoz_d(k+1)*aw(k,i), with MYNN's face density interpolation. MYNN mixes scalars into a
!plume with eps whatever its mass flux does, which can drain a layer of gas it never had where a plume
!accelerates. Here each plume is mass-consistent: in layer k,
!
!   E(k) = max(eps*M(k)/(1-eps), M(k+1)-M(k)),   D(k) = M(k) + E(k) - M(k+1) >= 0,
!   M(k+1)*c(k) = M(k)*c(k-1) + E(k)*q(k) - D(k)*c(k),  i.e.  c(k) = (M(k)*c(k-1) + E(k)*q(k))/(M(k)+E(k)),
!
!so c(k) is MYNN's (1-eps)*c(k-1) + eps*q(k) wherever MYNN's own fraction keeps D >= 0, and otherwise entrains
!just enough environmental air. A plume is rooted in layer 1 (M(1) = 0, so c(1) = q(1)); where it does not leave
!layer k (M(k+1) = 0) all of it detrains there and E(k) = 0.
!
!Downdrafts mirror this. MYNN exports, per downdraft j, aw(k,j) = DOWNA*(-DOWNW) on the bottom face of layer k
!(where its sd_aw(k) stores it) and eps(k,j); N(k,j) = rhoz_d(k)*aw(k,j) is the downward dry mass flux on face k,
!N(1) = 0 (MYNN keeps downdrafts off the ground). In layer k a downdraft enters through face k+1 and leaves
!through face k:
!
!   E(k) = max(eps*N(k+1)/(1-eps), N(k)-N(k+1)),   D(k) = N(k+1) + E(k) - N(k) >= 0,
!   c(k) = (N(k+1)*c(k+1) + E(k)*q(k))/(N(k+1)+E(k))   (concentration crossing face k downward),
!
!rooted in its start layer (N(k+1) = 0, so c(k) = q(k)), and detraining wholly in the layer it ends in.
!
!The plume concentrations are rebuilt from the current gases at every substep. The environment compensates each
!plume, upwind: downward from the layer above a face for updrafts, upward from the layer below for downdrafts.
!The upward flux through face k is
!
!   F(k) = sum_i M_i(k)*c_i(k-1) - Mtot(k)*q(k) - sum_j N_j(k)*c_j(k) + Ntot(k)*q(k-1),   F(1) = F(nz+1) = 0,
!
!and m_d(k)*dq(k)/dt = F(k) - F(k+1) = sum over plumes [D*c(k) - E*q(k)] - (Mtot(k) + Ntot(k+1))*q(k)
!+ Mtot(k+1)*q(k+1) + Ntot(k)*q(k-1). An explicit step of length dts is therefore a convex combination of q(k),
!its neighbours and the plume concentrations when dts*(Eu(k) + Mtot(k) + Ed(k) + Ntot(k+1)) <= m_d(k):
!monotone, nonnegative, uniform profiles unchanged and a zero species exactly zero. The step is split into the
!fewest equal substeps that satisfy that bound in every layer. The flux form conserves sum(m_d*q) to roundoff.


 real(kind=RKIND),parameter:: max_outgoing = 1.e4_RKIND   ! dt*max(outgoing/m_d) above this: bogus operator


 contains


!=================================================================================================================
 pure subroutine pbl_edmf_updraft_fluxes(aw,eps,rho_d,dz,M,E,ok)
!=================================================================================================================
!dry plume mass flux on faces and mass-consistent entrainment from MYNN's plume export. ok is false, and M and E
!are zero, unless every value is finite, aw >= 0, 0 <= eps < 1, aw(nz,:) = 0 (no flux through the model top),
!rho_d > 0 and dz > 0.
 use,intrinsic:: ieee_arithmetic,only: ieee_is_finite
 real(kind=RKIND),intent(in),dimension(:,:):: aw       ! (nz,np) [m s-1], top face of layer k
 real(kind=RKIND),intent(in),dimension(:,:):: eps      ! (nz,np) mixing fraction through layer k
 real(kind=RKIND),intent(in),dimension(:):: rho_d      ! (nz) dry air density [kg m-3]
 real(kind=RKIND),intent(in),dimension(:):: dz         ! (nz) layer thickness [m]
 real(kind=RKIND),intent(out),dimension(:,:):: M       ! (nz+1,np) dry mass flux on face k [kg m-2 s-1]
 real(kind=RKIND),intent(out),dimension(:,:):: E       ! (nz,np) entrainment in layer k [kg m-2 s-1]
 logical,intent(out):: ok
 real(kind=RKIND):: rhoz
 integer:: i,k,nz

 nz = size(rho_d)
 M(:,:) = 0._RKIND
 E(:,:) = 0._RKIND
 ok = all(ieee_is_finite(aw)) .and. all(ieee_is_finite(eps)) .and. all(ieee_is_finite(rho_d)) .and. &
      all(ieee_is_finite(dz))
 if(ok) ok = all(aw >= 0._RKIND) .and. all(eps >= 0._RKIND) .and. all(eps < 1._RKIND) .and. &
             all(aw(nz,:) == 0._RKIND) .and. all(rho_d > 0._RKIND) .and. all(dz > 0._RKIND)
 if(.not.ok) return

 do k = 1, nz-1
    rhoz = (rho_d(k+1)*dz(k) + rho_d(k)*dz(k+1))/(dz(k) + dz(k+1))
    M(k+1,:) = rhoz*aw(k,:)
 enddo
 do i = 1, size(aw,2)
    do k = 1, nz
       if(M(k+1,i) > 0._RKIND) E(k,i) = max(eps(k,i)*M(k,i)/(1._RKIND - eps(k,i)), M(k+1,i) - M(k,i))
    enddo
 enddo

 end subroutine pbl_edmf_updraft_fluxes

!=================================================================================================================
 pure subroutine pbl_edmf_updraft_plume(M,E,q,c)
!=================================================================================================================
!gas concentration c(:,k) of one plume as it leaves layer k through its top face (zero where it does not).
 real(kind=RKIND),intent(in),dimension(:):: M          ! (nz+1) plume mass flux on faces
 real(kind=RKIND),intent(in),dimension(:):: E          ! (nz) entrainment
 real(kind=RKIND),intent(in),dimension(:,:):: q        ! (nt,nz) environment
 real(kind=RKIND),intent(out),dimension(:,:):: c       ! (nt,nz)
 integer:: k

 c(:,:) = 0._RKIND
 do k = 1, size(E)
    if(M(k+1) > 0._RKIND) then
       if(M(k) > 0._RKIND) then
          c(:,k) = (M(k)*c(:,k-1) + E(k)*q(:,k))/(M(k) + E(k))
       else
          c(:,k) = q(:,k)                               ! rooted in layer k
       endif
    endif
 enddo

 end subroutine pbl_edmf_updraft_plume

!=================================================================================================================
 pure subroutine pbl_edmf_downdraft_fluxes(aw,eps,rho_d,dz,N,E,ok)
!=================================================================================================================
!downward dry mass flux on faces and mass-consistent entrainment from MYNN's downdraft export. ok is false, and
!N and E are zero, unless every value is finite, aw >= 0, 0 <= eps < 1, aw(1,:) = 0 (no flux into the ground),
!rho_d > 0 and dz > 0.
 use,intrinsic:: ieee_arithmetic,only: ieee_is_finite
 real(kind=RKIND),intent(in),dimension(:,:):: aw       ! (nz,nd) [m s-1], bottom face of layer k, downward
 real(kind=RKIND),intent(in),dimension(:,:):: eps      ! (nz,nd) mixing fraction through layer k
 real(kind=RKIND),intent(in),dimension(:):: rho_d      ! (nz) dry air density [kg m-3]
 real(kind=RKIND),intent(in),dimension(:):: dz         ! (nz) layer thickness [m]
 real(kind=RKIND),intent(out),dimension(:,:):: N       ! (nz+1,nd) downward dry mass flux on face k
 real(kind=RKIND),intent(out),dimension(:,:):: E       ! (nz,nd) entrainment in layer k [kg m-2 s-1]
 logical,intent(out):: ok
 real(kind=RKIND):: rhoz
 integer:: j,k,nz

 nz = size(rho_d)
 N(:,:) = 0._RKIND
 E(:,:) = 0._RKIND
 ok = all(ieee_is_finite(aw)) .and. all(ieee_is_finite(eps)) .and. all(ieee_is_finite(rho_d)) .and. &
      all(ieee_is_finite(dz))
 if(ok) ok = all(aw >= 0._RKIND) .and. all(eps >= 0._RKIND) .and. all(eps < 1._RKIND) .and. &
             all(aw(1,:) == 0._RKIND) .and. all(rho_d > 0._RKIND) .and. all(dz > 0._RKIND)
 if(.not.ok) return

 do k = 2, nz
    rhoz = (rho_d(k)*dz(k-1) + rho_d(k-1)*dz(k))/(dz(k-1) + dz(k))
    N(k,:) = rhoz*aw(k,:)
 enddo
 do j = 1, size(aw,2)
    do k = 1, nz
       if(N(k,j) > 0._RKIND) E(k,j) = max(eps(k,j)*N(k+1,j)/(1._RKIND - eps(k,j)), N(k,j) - N(k+1,j))
    enddo
 enddo

 end subroutine pbl_edmf_downdraft_fluxes

!=================================================================================================================
 pure subroutine pbl_edmf_downdraft_plume(N,E,q,c)
!=================================================================================================================
!gas concentration c(:,k) of one downdraft as it leaves layer k through its bottom face (zero where it does not).
 real(kind=RKIND),intent(in),dimension(:):: N          ! (nz+1) downward mass flux on faces
 real(kind=RKIND),intent(in),dimension(:):: E          ! (nz) entrainment
 real(kind=RKIND),intent(in),dimension(:,:):: q        ! (nt,nz) environment
 real(kind=RKIND),intent(out),dimension(:,:):: c       ! (nt,nz)
 integer:: k

 c(:,:) = 0._RKIND
 do k = size(E), 1, -1
    if(N(k) > 0._RKIND) then
       if(N(k+1) > 0._RKIND) then
          c(:,k) = (N(k+1)*c(:,k+1) + E(k)*q(:,k))/(N(k+1) + E(k))
       else
          c(:,k) = q(:,k)                               ! rooted in layer k
       endif
    endif
 enddo

 end subroutine pbl_edmf_downdraft_plume

!=================================================================================================================
 pure subroutine pbl_edmf_column(dt,m_d,Mu,Eu,Nd,Ed,q,nsub,ok)
!=================================================================================================================
!one model step of updraft and downdraft transport of every gas in a closed column, in nsub equal substeps. ok is
!false, and q is unchanged, when dt*max(outgoing/m_d) exceeds 1e4 or is not finite, or m_d <= 0. Without
!downdrafts, pass Nd and Ed with zero columns or values.
 use,intrinsic:: ieee_arithmetic,only: ieee_is_finite
 real(kind=RKIND),intent(in):: dt                      ! model time step [s]
 real(kind=RKIND),intent(in),dimension(:):: m_d        ! (nz) dry layer mass [kg m-2]
 real(kind=RKIND),intent(in),dimension(:,:):: Mu       ! (nz+1,np) from pbl_edmf_updraft_fluxes
 real(kind=RKIND),intent(in),dimension(:,:):: Eu       ! (nz,np)
 real(kind=RKIND),intent(in),dimension(:,:):: Nd       ! (nz+1,nd) from pbl_edmf_downdraft_fluxes
 real(kind=RKIND),intent(in),dimension(:,:):: Ed       ! (nz,nd)
 real(kind=RKIND),intent(inout),dimension(:,:):: q     ! (nt,nz) [kg/kg]
 integer,intent(out):: nsub
 logical,intent(out):: ok

 real(kind=RKIND),dimension(size(q,1),size(q,2)):: c
 real(kind=RKIND),dimension(size(q,1),size(q,2)+1):: F
 real(kind=RKIND),dimension(size(m_d)+1):: mtot,ntot
 real(kind=RKIND):: courant,dts
 integer:: i,k,nz,isub

 nz = size(m_d)
 nsub = 0
 ok = all(m_d > 0._RKIND)
 if(.not.ok) return
 mtot(:) = sum(Mu,dim=2)
 ntot(:) = sum(Nd,dim=2)
 courant = dt*maxval((sum(Eu,dim=2) + mtot(1:nz) + sum(Ed,dim=2) + ntot(2:nz+1))/m_d)
 ok = ieee_is_finite(courant) .and. courant <= max_outgoing
 if(.not.ok) return
 nsub = max(1,ceiling(courant))
 dts = dt/nsub

 do isub = 1, nsub
    F(:,:) = 0._RKIND
    do i = 1, size(Mu,2)
       if(all(Mu(:,i) == 0._RKIND)) cycle
       call pbl_edmf_updraft_plume(Mu(:,i),Eu(:,i),q,c)
       do k = 1, nz-1
          if(Mu(k+1,i) > 0._RKIND) F(:,k+1) = F(:,k+1) + Mu(k+1,i)*c(:,k)
       enddo
    enddo
    do i = 1, size(Nd,2)
       if(all(Nd(:,i) == 0._RKIND)) cycle
       call pbl_edmf_downdraft_plume(Nd(:,i),Ed(:,i),q,c)
       do k = 2, nz
          if(Nd(k,i) > 0._RKIND) F(:,k) = F(:,k) - Nd(k,i)*c(:,k)
       enddo
    enddo
    do k = 2, nz
       F(:,k) = F(:,k) - mtot(k)*q(:,k) + ntot(k)*q(:,k-1)
    enddo
    do k = 1, nz
       q(:,k) = q(:,k) + dts*(F(:,k) - F(:,k+1))/m_d(k)
    enddo
 enddo

 end subroutine pbl_edmf_column

!=================================================================================================================
 end module mpas_atmphys_pbl_edmf
!=================================================================================================================
