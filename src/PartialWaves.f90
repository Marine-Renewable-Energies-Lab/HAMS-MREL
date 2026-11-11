!  ------------------------------------------------------------------------------------------------------
!
!    Program HAMS-MREL for the diffraction and radiation of waves
!    by 3D structures.
!
!  License:
!
!    This routine is part of HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Vaibhav Raghavan
!
!  ------------------------------------------------------------------------------------------------------
MODULE PartialWaves
   ! Partial-wave decomposition of a first-order field outside a vertical cylinder r = R about centre, in water of finite
   ! depth h (Phase 6a; reused in Phase 8). Outside R, a field that satisfies Laplace, the free-surface condition
   ! -omega^2 phi + g phi_z = 0, the bottom condition and the radiation condition is (time factor e^{-i omega t})
   !
   !   phi(r, theta, z) = sum_{m=-M..M} e^{i m theta} [ a(m, 0) H_m(k0 r) Z_0(z) + sum_{n=1..N} a(m, n) K_m(kn r) Z_n(z) ]
   !
   ! with theta measured from +x about centre, H_m = H_m^(1) = J_m + i Y_m, and the vertical eigenfunctions
   !   Z_0 = cosh k0(z+h) / cosh k0h,   Z_n = cos kn(z+h) / cos knh,   Z(0) = 1,
   ! where omega^2/g = k0 tanh k0h = -kn tan knh (the roots DISPERSION puts in WVN). H^(1) radiates outwards with
   ! e^{-i omega t}: H_m^(1)(kr) ~ e^{i(kr - ...)}, and the finite-depth Green function of HAMS has the far field
   ! P0 (i J_0 - Y_0) = i P0 H_0^(1)(kr) (FinGreen3D.f90, EIGEN). The incident wave is e^{+i k x} (incident_phi).
   !
   ! The field is sampled on the cylinder at N_theta = 4M + 4 equally spaced angles and N_z = 2(N + 1) Gauss-Legendre depths
   ! on [-h, 0]. A DFT in theta gives the angular modes at each depth; at each mode, a least-squares fit of Z_0..Z_N to the
   ! N_z depths (weighted by the square roots of the Gauss weights, so that it tends to the orthogonal projection on
   ! [-h, 0]) gives the mode amplitudes at r = R, which are divided by H_m(k0 R) and K_m(kn R).

   USE SecondOrderTypes, ONLY : cs_points, CS_WALL
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: pw_set, pw_sample_points, pw_decompose, pw_eval_surface, pw_bessel_k, pw_evanescent_from_wvn

   REAL(8), PARAMETER :: PI = 3.14159265358979323846D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   ! K_m(kn R) below this underflows in the coefficients (K_0(x) ~ e^{-x}); the decomposition stops rather than divide by it
   REAL(8), PARAMETER :: KR_MAX = 600D0

   TYPE :: pw_set
      ! The modes of one field. a(m, 0) multiplies H_m(k0 r) Z_0(z), a(m, n) multiplies K_m(kn r) Z_n(z), n = 1..nev.
      INTEGER :: mmax = 0, nev = 0
      REAL(8) :: k0 = 0D0, h = 0D0, radius = 0D0, centre(2) = 0D0
      REAL(8), ALLOCATABLE :: kev(:)          ! (nev)
      COMPLEX(8), ALLOCATABLE :: a(:, :)      ! (-mmax:mmax, 0:nev)
   END TYPE pw_set

CONTAINS

   SUBROUTINE pw_sample_points(centre, radius, h, mmax, nev, pts)
      ! The sample points on the cylinder r = radius about centre: N_theta = 4 mmax + 4 angles theta_j = 2 pi j / N_theta
      ! (j = 0..N_theta - 1) and N_z = 2 (nev + 1) Gauss-Legendre depths z_l on [-h, 0] (l = 1..N_z, from the bottom up).
      ! Point i = j + 1 + N_theta (l - 1): theta runs fastest. Kind CS_WALL, normal outward (away from centre), w the
      ! quadrature area R dtheta w_l, psize the larger of the arc spacing and h / N_z, body 0 (the circle is about all bodies).
      REAL(8), INTENT(IN) :: centre(2), radius, h
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(cs_points), INTENT(OUT) :: pts
      REAL(8), ALLOCATABLE :: z(:), wz(:)
      REAL(8) :: th, dth
      INTEGER :: nth, nz, j, l, i

      nth = 4*mmax + 4
      nz = 2*(nev + 1)
      CALL depth_nodes(h, nz, z, wz)
      dth = 2D0*PI/nth
      pts%n = nth*nz
      ALLOCATE (pts%x(3, pts%n), pts%nrm(3, pts%n), pts%w(pts%n), pts%psize(pts%n), pts%kind(pts%n), pts%body(pts%n))
      DO l = 1, nz
         DO j = 0, nth - 1
            i = j + 1 + nth*(l - 1)
            th = dth*j
            pts%x(:, i) = [centre(1) + radius*COS(th), centre(2) + radius*SIN(th), z(l)]
            pts%nrm(:, i) = [COS(th), SIN(th), 0D0]
            pts%w(i) = radius*dth*wz(l)
         END DO
      END DO
      pts%psize = MAX(radius*dth, h/nz)
      pts%kind = CS_WALL
      pts%body = 0
   END SUBROUTINE pw_sample_points

   SUBROUTINE pw_decompose(phi_samples, centre, radius, h, k0, kev, mmax, nev, pw)
      ! The modes of the field phi_samples, given at the points of pw_sample_points(centre, radius, h, mmax, nev) in their
      ! order. k0 is the propagating wave number and kev(1:nev) the evanescent ones (WVN(2:nev+1), pw_evanescent_from_wvn).
      COMPLEX(8), INTENT(IN) :: phi_samples(:)
      REAL(8), INTENT(IN) :: centre(2), radius, h, k0, kev(:)
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(pw_set), INTENT(OUT) :: pw
      REAL(8), ALLOCATABLE :: z(:), wz(:), amat(:, :), bmat(:, :), work(:), kvals(:)
      COMPLEX(8), ALLOCATABLE :: c(:, :)      ! (-mmax:mmax, nz), the angular modes at each depth
      COMPLEX(8) :: hk, s
      REAL(8) :: th, sq
      INTEGER :: nth, nz, ncol, nrhs, j, l, m, n, info, lwork
      EXTERNAL :: DGELS

      nth = 4*mmax + 4
      nz = 2*(nev + 1)
      ncol = nev + 1
      IF (SIZE(phi_samples) /= nth*nz) ERROR STOP "pw_decompose: phi_samples is not the size of pw_sample_points"
      IF (SIZE(kev) < nev) ERROR STOP "pw_decompose: fewer evanescent wave numbers than nev"
      IF (nev > 0) THEN
         IF (kev(nev)*radius > KR_MAX) ERROR STOP "pw_decompose: kn R too large for K_m(kn R) (lower nev or the radius)"
      END IF

      pw%mmax = mmax
      pw%nev = nev
      pw%k0 = k0
      pw%h = h
      pw%radius = radius
      pw%centre = centre
      ALLOCATE (pw%kev(nev), pw%a(-mmax:mmax, 0:nev))
      pw%kev = kev(1:nev)

      ! DFT in theta at each depth: c(m, l) = (1/N_theta) sum_j phi(theta_j, z_l) e^{-i m theta_j}
      ALLOCATE (c(-mmax:mmax, nz))
      DO l = 1, nz
         DO m = -mmax, mmax
            s = (0D0, 0D0)
            DO j = 0, nth - 1
               th = 2D0*PI*MODULO(m*j, nth)/nth
               s = s + phi_samples(j + 1 + nth*(l - 1))*CMPLX(COS(th), -SIN(th), KIND=8)
            END DO
            c(m, l) = s/nth
         END DO
      END DO

      ! Weighted least squares at every m at once: sqrt(w_l) Z_n(z_l) b(m, n) = sqrt(w_l) c(m, l), real and imaginary parts
      ! as separate right-hand sides (the matrix is real)
      CALL depth_nodes(h, nz, z, wz)
      nrhs = 2*(2*mmax + 1)
      ALLOCATE (amat(nz, ncol), bmat(nz, nrhs))
      DO l = 1, nz
         sq = SQRT(wz(l))
         amat(l, 1) = sq*z0_mode(k0, h, z(l))
         DO n = 1, nev
            amat(l, n + 1) = sq*COS(kev(n)*(z(l) + h))/COS(kev(n)*h)
         END DO
         DO m = -mmax, mmax
            bmat(l, 2*(m + mmax) + 1) = sq*REAL(c(m, l))
            bmat(l, 2*(m + mmax) + 2) = sq*AIMAG(c(m, l))
         END DO
      END DO
      lwork = -1
      ALLOCATE (work(1))
      CALL DGELS('N', nz, ncol, nrhs, amat, nz, bmat, nz, work, lwork, info)
      lwork = MAX(1, INT(work(1)))
      DEALLOCATE (work)
      ALLOCATE (work(lwork))
      CALL DGELS('N', nz, ncol, nrhs, amat, nz, bmat, nz, work, lwork, info)
      IF (info /= 0) ERROR STOP "pw_decompose: the least-squares fit of the vertical modes failed (DGELS)"

      ! Divide by the radial functions at R: H_{-m} = (-1)^m H_m, K_{-m} = K_m
      ALLOCATE (kvals(0:mmax))
      DO m = 0, mmax
         hk = CMPLX(BESSEL_JN(m, k0*radius), BESSEL_YN(m, k0*radius), KIND=8)
         pw%a(m, 0) = CMPLX(bmat(1, 2*(m + mmax) + 1), bmat(1, 2*(m + mmax) + 2), KIND=8)/hk
         IF (m > 0) pw%a(-m, 0) = CMPLX(bmat(1, 2*(mmax - m) + 1), bmat(1, 2*(mmax - m) + 2), KIND=8) &
                                  /(MERGE(-1D0, 1D0, MODULO(m, 2) == 1)*hk)
      END DO
      DO n = 1, nev
         CALL pw_bessel_k(mmax, kev(n)*radius, kvals)
         DO m = -mmax, mmax
            pw%a(m, n) = CMPLX(bmat(n + 1, 2*(m + mmax) + 1), bmat(n + 1, 2*(m + mmax) + 2), KIND=8)/kvals(ABS(m))
         END DO
      END DO
   END SUBROUTINE pw_decompose

   PURE SUBROUTINE pw_eval_surface(pw, x, phi, grad, phizz)
      ! phi, grad phi and phi_zz on z = 0 at the horizontal position x(1:2) (x(3) is not used), from the modes of pw. For
      ! r >= R only, where the expansion holds. phi_z = Z'(0) phi per mode: k0 tanh k0h and -kn tan knh, both omega^2/g;
      ! phi_zz = Z''(0) phi per mode: k0^2 and -kn^2.
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(OUT) :: phi, grad(3), phizz
      COMPLEX(8) :: hk(0:pw%mmax + 1), rad(-pw%mmax:pw%mmax), drad(-pw%mmax:pw%mmax), e, sgn, dr, dth
      COMPLEX(8) :: phiz
      REAL(8) :: kv(0:pw%mmax + 1), dxy(2), r, th, xk, dz0, dzn
      INTEGER :: m, n, am, mmax

      mmax = pw%mmax
      dxy = x(1:2) - pw%centre
      r = NORM2(dxy)
      th = ATAN2(dxy(2), dxy(1))

      phi = (0D0, 0D0)
      dr = (0D0, 0D0)
      dth = (0D0, 0D0)
      phiz = (0D0, 0D0)
      phizz = (0D0, 0D0)

      ! Propagating mode: H_m and H_m' = (H_{m-1} - H_{m+1}) / 2, H_{-1} = -H_1
      xk = pw%k0*r
      DO m = 0, mmax + 1
         hk(m) = CMPLX(BESSEL_JN(m, xk), BESSEL_YN(m, xk), KIND=8)
      END DO
      dz0 = pw%k0*TANH(pw%k0*pw%h)
      DO m = -mmax, mmax
         am = ABS(m)
         sgn = MERGE(-1D0, 1D0, m < 0 .AND. MODULO(am, 2) == 1)
         rad(m) = sgn*hk(am)
         IF (am == 0) THEN
            drad(m) = -hk(1)*pw%k0
         ELSE
            drad(m) = sgn*0.5D0*(hk(am - 1) - hk(am + 1))*pw%k0
         END IF
      END DO
      DO m = -mmax, mmax
         e = EXP(CI*(m*th))*pw%a(m, 0)
         phi = phi + e*rad(m)
         dr = dr + e*drad(m)
         dth = dth + CI*m*e*rad(m)
         phiz = phiz + dz0*e*rad(m)
         phizz = phizz + pw%k0**2*e*rad(m)
      END DO

      ! Evanescent modes: K_m and K_m' = -(K_{m-1} + K_{m+1}) / 2, K_{-1} = K_1
      DO n = 1, pw%nev
         xk = pw%kev(n)*r
         CALL pw_bessel_k(mmax + 1, xk, kv)
         dzn = -pw%kev(n)*TAN(pw%kev(n)*pw%h)
         DO m = -mmax, mmax
            am = ABS(m)
            e = EXP(CI*(m*th))*pw%a(m, n)
            IF (am == 0) THEN
               drad(m) = -kv(1)*pw%kev(n)
            ELSE
               drad(m) = -0.5D0*(kv(am - 1) + kv(am + 1))*pw%kev(n)
            END IF
            phi = phi + e*kv(am)
            dr = dr + e*drad(m)
            dth = dth + CI*m*e*kv(am)
            phiz = phiz + dzn*e*kv(am)
            phizz = phizz - pw%kev(n)**2*e*kv(am)
         END DO
      END DO

      grad(1) = COS(th)*dr - SIN(th)*dth/r
      grad(2) = SIN(th)*dr + COS(th)*dth/r
      grad(3) = phiz
   END SUBROUTINE pw_eval_surface

   PURE SUBROUTINE pw_bessel_k(mmax, x, k)
      ! Modified Bessel functions of the second kind K_0(x)..K_mmax(x), x > 0. K_0 and K_1 by their power series for x <= 2
      ! (Abramowitz and Stegun 9.6.13 and 9.6.11) and by Steed's continued fraction CF2 with Temme's normalisation for x > 2
      ! (Temme 1975, J. Comput. Phys. 19, 324; Numerical Recipes, 2nd ed., section 6.7), then the upward recurrence
      ! K_{m+1} = K_{m-1} + (2m/x) K_m, which is stable for K. Relative accuracy about 1e-15 (pw_bessel_k_values).
      ! K underflows to zero for x above about 700.
      INTEGER, INTENT(IN) :: mmax
      REAL(8), INTENT(IN) :: x
      REAL(8), INTENT(OUT) :: k(0:mmax)
      REAL(8), PARAMETER :: EULER = 0.57721566490153286061D0, EPS = 1D-17
      REAL(8) :: t, lg, psi1, psi2, term0, term1, s0, s1, si1, k0, k1
      REAL(8) :: b, d, hh, delh, q1, q2, qnew, q, cc, a, a1, s, dels
      INTEGER :: i, m

      IF (x <= 2D0) THEN
         t = 0.25D0*x*x
         lg = LOG(0.5D0*x)
         psi1 = -EULER
         psi2 = 1D0 - EULER
         term0 = 1D0
         term1 = 1D0
         s0 = psi1 - lg
         s1 = psi1 + psi2
         si1 = 1D0
         DO i = 1, 100
            psi1 = psi2
            psi2 = psi1 + 1D0/(i + 1)
            term0 = term0*t/(REAL(i, 8)*i)
            term1 = term1*t/(REAL(i, 8)*(i + 1))
            s0 = s0 + term0*(psi1 - lg)
            s1 = s1 + term1*(psi1 + psi2)
            si1 = si1 + term1
            IF (term0 < EPS*ABS(s0) .AND. term1 < EPS) EXIT
         END DO
         k0 = s0
         k1 = 1D0/x + lg*0.5D0*x*si1 - 0.25D0*x*s1
      ELSE
         b = 2D0*(1D0 + x)
         d = 1D0/b
         hh = d
         delh = d
         q1 = 0D0
         q2 = 1D0
         a1 = 0.25D0
         q = a1
         cc = a1
         a = -a1
         s = 1D0 + q*delh
         DO i = 2, 10000
            a = a - 2D0*(i - 1)
            cc = -a*cc/i
            qnew = (q1 - b*q2)/a
            q1 = q2
            q2 = qnew
            q = q + cc*qnew
            b = b + 2D0
            d = 1D0/(b + a*d)
            delh = (b*d - 1D0)*delh
            hh = hh + delh
            dels = q*delh
            s = s + dels
            IF (ABS(dels/s) < EPS) EXIT
         END DO
         hh = a1*hh
         k0 = SQRT(PI/(2D0*x))*EXP(-x)/s
         k1 = k0*(x + 0.5D0 - hh)/x
      END IF

      k(0) = k0
      IF (mmax >= 1) k(1) = k1
      DO m = 1, mmax - 1
         k(m + 1) = k(m - 1) + (2D0*m/x)*k(m)
      END DO
   END SUBROUTINE pw_bessel_k

   SUBROUTINE pw_evanescent_from_wvn(nev, kev)
      ! The first nev evanescent wave numbers of the current frequency, WVN(2:nev+1) (DISPERSION, WaveDyn_mod): the roots of
      ! -omega^2/g = kn tan(kn h) that the Green function of the solve used. WVN holds NK - 1 of them.
      USE WaveDyn_mod, ONLY : WVN, NK
      INTEGER, INTENT(IN) :: nev
      REAL(8), INTENT(OUT) :: kev(nev)
      CHARACTER(len=160) :: msg

      IF (nev > NK - 1) THEN
         WRITE (msg, '(a,i0,a,i0,a)') "pw_evanescent_from_wvn: ", nev, " evanescent modes asked, but WVN holds only ", NK - 1, &
            " roots"
         ERROR STOP TRIM(msg)
      END IF
      kev = WVN(2:nev + 1)
   END SUBROUTINE pw_evanescent_from_wvn

   PURE REAL(8) FUNCTION z0_mode(k0, h, z)
      ! cosh k0(z+h) / cosh k0h without overflow at large k0 h
      REAL(8), INTENT(IN) :: k0, h, z
      z0_mode = (EXP(k0*z) + EXP(-k0*(z + 2D0*h)))/(1D0 + EXP(-2D0*k0*h))
   END FUNCTION z0_mode

   SUBROUTINE depth_nodes(h, nz, z, wz)
      ! Gauss-Legendre nodes and weights on [-h, 0], from the bottom up: Newton on P_nz from the Chebyshev-like first guess.
      REAL(8), INTENT(IN) :: h
      INTEGER, INTENT(IN) :: nz
      REAL(8), ALLOCATABLE, INTENT(OUT) :: z(:), wz(:)
      REAL(8) :: t, p0, p1, p2, dp, dt
      INTEGER :: i, j, it

      ALLOCATE (z(nz), wz(nz))
      DO i = 1, (nz + 1)/2
         t = COS(PI*(i - 0.25D0)/(nz + 0.5D0))
         DO it = 1, 100
            p0 = 1D0
            p1 = t
            DO j = 2, nz
               p2 = ((2*j - 1)*t*p1 - (j - 1)*p0)/j
               p0 = p1
               p1 = p2
            END DO
            dp = nz*(t*p1 - p0)/(t*t - 1D0)
            dt = p1/dp
            t = t - dt
            IF (ABS(dt) < 1D-15) EXIT
         END DO
         ! recompute at the converged node
         p0 = 1D0
         p1 = t
         DO j = 2, nz
            p2 = ((2*j - 1)*t*p1 - (j - 1)*p0)/j
            p0 = p1
            p1 = p2
         END DO
         dp = nz*(t*p1 - p0)/(t*t - 1D0)
         ! t descends from near +1: node i from the top is z(nz + 1 - i), its mirror from the bottom z(i)
         z(nz + 1 - i) = 0.5D0*h*(t - 1D0)
         z(i) = 0.5D0*h*(-t - 1D0)
         wz(nz + 1 - i) = h/((1D0 - t*t)*dp*dp)
         wz(i) = wz(nz + 1 - i)
      END DO
   END SUBROUTINE depth_nodes

END MODULE PartialWaves
