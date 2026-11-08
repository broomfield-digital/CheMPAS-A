! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_pbl_diffusion
 use mpas_kind_types,only: RKIND

 implicit none
 private
 public:: pbl_diffuse_column


!Closed-column implicit (backward-Euler) vertical diffusion of MICM gases (P2 of docs/plans/pbl-transport.md).
!For layer k (MPAS k = 1 at the ground) with dry layer mass m_d(k) = zz*rho_zz*dz [kg m-2]:
!
!   m_d(k)*(q1(k)-q0(k))/dt = G(k+1)*(q1(k+1)-q1(k)) - G(k)*(q1(k)-q1(k-1))
!
!where G(k) = rho*Kz/dz [kg m-2 s-1] at the interface below layer k. The ground (G(1)) and top (G(nz+1))
!interfaces are closed whatever values the caller supplies: no emission, deposition or escape.
!
!The step is solved for the increment d = q1 - q0, whose right-hand side is the flux divergence of q0, with one
!Thomas factorization of the M-matrix for every species. The interface fluxes of q1 = q0 + d are then applied
!to q0 in flux form. So, to roundoff at any stiffness: sum(m_d*q) telescopes exactly; a uniform profile has
!zero right-hand side and stays exactly unchanged; values stay within their initial range (monotone, no new
!extrema, nonnegative); a zero species stays exactly zero (there is no positive floor).
!
!Preconditions (the caller checks them): m_d > 0, G >= 0 and finite; size(G) = size(m_d)+1 = size(q,2)+1.


 contains


!=================================================================================================================
 pure subroutine pbl_diffuse_column(dt,m_d,G,q)
!=================================================================================================================
 real(kind=RKIND),intent(in):: dt                       ! time step [s]
 real(kind=RKIND),intent(in),dimension(:):: m_d         ! (nz) dry layer mass [kg m-2]
 real(kind=RKIND),intent(in),dimension(:):: G           ! (nz+1) interface exchange [kg m-2 s-1]
 real(kind=RKIND),intent(inout),dimension(:,:):: q      ! (nt,nz) mixing ratio [kg/kg dry air]

 real(kind=RKIND),dimension(size(m_d)):: a,b,c,cp
 real(kind=RKIND),dimension(size(q,1),size(q,2)):: dp,d
 real(kind=RKIND),dimension(size(q,1),size(q,2)+1):: F
 real(kind=RKIND):: g_below,g_above,w
 integer:: k,nz

 nz = size(m_d)

!tridiagonal row k for the increment: a(k)*d(k-1) + b(k)*d(k) + c(k)*d(k+1) = dt*div(F(q0))(k)/m_d(k).
 do k = 1, nz
    g_below = 0._RKIND
    g_above = 0._RKIND
    if(k > 1 ) g_below = G(k)
    if(k < nz) g_above = G(k+1)
    a(k) = -dt*g_below/m_d(k)
    c(k) = -dt*g_above/m_d(k)
    b(k) = 1._RKIND - a(k) - c(k)
 enddo

!upward-positive interface fluxes of the old state [kg m-2 s-1]; closed ground and top.
 F(:,1)    = 0._RKIND
 F(:,nz+1) = 0._RKIND
 do k = 2, nz
    F(:,k) = G(k)*(q(:,k-1)-q(:,k))
 enddo

!Thomas sweep for the increment; the matrix is diagonally dominant, so no pivoting is needed.
 cp(1)   = c(1)/b(1)
 dp(:,1) = dt*(F(:,1)-F(:,2))/m_d(1)/b(1)
 do k = 2, nz
    w       = b(k) - a(k)*cp(k-1)
    cp(k)   = c(k)/w
    dp(:,k) = (dt*(F(:,k)-F(:,k+1))/m_d(k) - a(k)*dp(:,k-1))/w
 enddo
 d(:,nz) = dp(:,nz)
 do k = nz-1, 1, -1
    d(:,k) = dp(:,k) - cp(k)*d(:,k+1)
 enddo

!flux-form finish: the interface fluxes of q1 = q0 + d, applied to q0, so the column mass telescopes exactly.
 do k = 2, nz
    F(:,k) = F(:,k) + G(k)*(d(:,k-1)-d(:,k))
 enddo
 do k = 1, nz
    q(:,k) = q(:,k) + dt*(F(:,k)-F(:,k+1))/m_d(k)
 enddo

 end subroutine pbl_diffuse_column

!=================================================================================================================
 end module mpas_atmphys_pbl_diffusion
!=================================================================================================================
