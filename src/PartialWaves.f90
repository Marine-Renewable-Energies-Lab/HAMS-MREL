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

   PUBLIC :: pw_set, pw_radial, pw_sample_points, pw_decompose, pw_eval_surface, pw_eval_point, pw_radial_at, pw_bessel_k
   PUBLIC :: pw_evanescent_from_wvn, PW_KR_MAX, pw_bessel_k_scaled, pw_bessel_k_scaled_c, pw_hankel_scaled

   REAL(8), PARAMETER :: PI = 3.14159265358979323846D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   ! K_m(kn R) below this underflows in the coefficients (K_0(x) ~ e^{-x}); the decomposition stops rather than divide by it
   REAL(8), PARAMETER :: KR_MAX = 600D0, PW_KR_MAX = KR_MAX

   TYPE :: pw_set
      ! The modes of one field. a(m, 0) multiplies H_m(k0 r) Z_0(z), a(m, n) multiplies K_m(kn r) Z_n(z), n = 1..nev.
      INTEGER :: mmax = 0, nev = 0
      REAL(8) :: k0 = 0D0, h = 0D0, radius = 0D0, centre(2) = 0D0
      REAL(8), ALLOCATABLE :: kev(:)          ! (nev)
      COMPLEX(8), ALLOCATABLE :: a(:, :)      ! (-mmax:mmax, 0:nev)
   END TYPE pw_set

   TYPE :: pw_radial
      ! The radial functions of a pw_set at one radius r (pw_radial_at): H_m(k0 r) and K_m(kn r) for m = 0..mmax + 1
      REAL(8) :: r = -1D0
      INTEGER :: mmax = -1, nev = -1
      COMPLEX(8), ALLOCATABLE :: hk(:)        ! (0:mmax+1)
      REAL(8), ALLOCATABLE :: kv(:, :)        ! (0:mmax+1, nev)
   END TYPE pw_radial

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

   PURE SUBROUTINE pw_radial_at(pw, r, rad)
      ! The radial functions of the modes of pw at r: H_m(k0 r) and K_m(kn r), m = 0..mmax + 1 (the derivatives need
      ! mmax + 1). They depend on r only, so a caller that evaluates many points at one radius (a ring, or the samples of
      ! the circle) computes them once and passes rad to pw_eval_surface or pw_eval_point.
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: r
      TYPE(pw_radial), INTENT(OUT) :: rad
      INTEGER :: m, n

      rad%r = r
      rad%mmax = pw%mmax
      rad%nev = pw%nev
      ALLOCATE (rad%hk(0:pw%mmax + 1), rad%kv(0:pw%mmax + 1, pw%nev))
      DO m = 0, pw%mmax + 1
         rad%hk(m) = CMPLX(BESSEL_JN(m, pw%k0*r), BESSEL_YN(m, pw%k0*r), KIND=8)
      END DO
      DO n = 1, pw%nev
         CALL pw_bessel_k(pw%mmax + 1, pw%kev(n)*r, rad%kv(:, n))
      END DO
   END SUBROUTINE pw_radial_at

   PURE SUBROUTINE pw_eval_point(pw, x, phi, rad)
      ! phi at the point x from the modes of pw, at any depth -h <= x(3) <= 0 and horizontally at r >= R, where the
      ! expansion holds: sum_m e^{i m theta} [ a(m, 0) H_m(k0 r) Z_0(z) + sum_n a(m, n) K_m(kn r) Z_n(z) ]. At z = 0 it is the
      ! phi of pw_eval_surface. For the mode check on the circle samples, which lie at depth. rad as for pw_eval_surface.
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(OUT) :: phi
      TYPE(pw_radial), INTENT(IN), OPTIONAL :: rad
      TYPE(pw_radial) :: own
      COMPLEX(8) :: s, f, df, sgn
      REAL(8) :: dxy(2), r, th, zf
      INTEGER :: m, n, am

      dxy = x(1:2) - pw%centre
      r = NORM2(dxy)
      th = ATAN2(dxy(2), dxy(1))
      CALL radial_values(pw, r, own, rad)

      phi = (0D0, 0D0)
      DO n = 0, pw%nev
         IF (n == 0) THEN
            zf = z0_mode(pw%k0, pw%h, x(3))
         ELSE
            zf = COS(pw%kev(n)*(x(3) + pw%h))/COS(pw%kev(n)*pw%h)
         END IF
         s = (0D0, 0D0)
         DO m = -pw%mmax, pw%mmax
            am = ABS(m)
            IF (n == 0) THEN
               sgn = MERGE(-1D0, 1D0, m < 0 .AND. MODULO(am, 2) == 1)
               CALL hankel_pair(own%hk, am, sgn, pw%k0, f, df)
            ELSE
               CALL bessel_k_pair(own%kv(:, n), am, pw%kev(n), f, df)
            END IF
            s = s + EXP(CI*(m*th))*pw%a(m, n)*f
         END DO
         phi = phi + zf*s
      END DO
   END SUBROUTINE pw_eval_point

   PURE SUBROUTINE pw_eval_surface(pw, x, phi, grad, phizz, rad)
      ! phi, grad phi and phi_zz on z = 0 at the horizontal position x(1:2) (x(3) is not used), from the modes of pw. For
      ! r >= R only, where the expansion holds. phi_z = Z'(0) phi per mode: k0 tanh k0h and -kn tan knh, both omega^2/g;
      ! phi_zz = Z''(0) phi per mode: k0^2 and -kn^2. rad, if given and made by pw_radial_at for pw at this r, holds the
      ! radial functions; otherwise they are computed here.
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(OUT) :: phi, grad(3), phizz
      TYPE(pw_radial), INTENT(IN), OPTIONAL :: rad
      TYPE(pw_radial) :: own
      COMPLEX(8) :: rd, drd, e, sgn, dr, dth, phiz
      REAL(8) :: dxy(2), r, th, dz0, dzn
      INTEGER :: m, n, am

      dxy = x(1:2) - pw%centre
      r = NORM2(dxy)
      th = ATAN2(dxy(2), dxy(1))
      ! The radial functions: those of rad if they are at this r, else computed here (own is then a copy, so that no
      ! branch below refers to rad, which may be absent)
      CALL radial_values(pw, r, own, rad)

      phi = (0D0, 0D0)
      dr = (0D0, 0D0)
      dth = (0D0, 0D0)
      phiz = (0D0, 0D0)
      phizz = (0D0, 0D0)

      ! Propagating mode: H_m and H_m' = (H_{m-1} - H_{m+1}) / 2, H_{-1} = -H_1
      dz0 = pw%k0*TANH(pw%k0*pw%h)
      DO m = -pw%mmax, pw%mmax
         am = ABS(m)
         sgn = MERGE(-1D0, 1D0, m < 0 .AND. MODULO(am, 2) == 1)
         CALL hankel_pair(own%hk, am, sgn, pw%k0, rd, drd)
         e = EXP(CI*(m*th))*pw%a(m, 0)
         phi = phi + e*rd
         dr = dr + e*drd
         dth = dth + CI*m*e*rd
         phiz = phiz + dz0*e*rd
         phizz = phizz + pw%k0**2*e*rd
      END DO

      ! Evanescent modes: K_m and K_m' = -(K_{m-1} + K_{m+1}) / 2, K_{-1} = K_1
      DO n = 1, pw%nev
         dzn = -pw%kev(n)*TAN(pw%kev(n)*pw%h)
         DO m = -pw%mmax, pw%mmax
            am = ABS(m)
            CALL bessel_k_pair(own%kv(:, n), am, pw%kev(n), rd, drd)
            e = EXP(CI*(m*th))*pw%a(m, n)
            phi = phi + e*rd
            dr = dr + e*drd
            dth = dth + CI*m*e*rd
            phiz = phiz + dzn*e*rd
            phizz = phizz - pw%kev(n)**2*e*rd
         END DO
      END DO

      grad(1) = COS(th)*dr - SIN(th)*dth/r
      grad(2) = SIN(th)*dr + COS(th)*dth/r
      grad(3) = phiz
   END SUBROUTINE pw_eval_surface

   PURE SUBROUTINE radial_values(pw, r, own, rad)
      ! own = rad if rad is present and holds the radial functions of pw at r, else pw_radial_at(pw, r)
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: r
      TYPE(pw_radial), INTENT(OUT) :: own
      TYPE(pw_radial), INTENT(IN), OPTIONAL :: rad
      LOGICAL :: cached

      cached = .FALSE.
      IF (PRESENT(rad)) cached = radial_fits(pw, r, rad)
      IF (cached) THEN
         own = rad
      ELSE
         CALL pw_radial_at(pw, r, own)
      END IF
   END SUBROUTINE radial_values

   PURE LOGICAL FUNCTION radial_fits(pw, r, rad)
      ! rad holds the radial functions of pw at r (to rounding)
      TYPE(pw_set), INTENT(IN) :: pw
      REAL(8), INTENT(IN) :: r
      TYPE(pw_radial), INTENT(IN) :: rad

      radial_fits = .FALSE.
      IF (.NOT. (ALLOCATED(rad%hk) .AND. ALLOCATED(rad%kv))) RETURN
      radial_fits = rad%mmax == pw%mmax .AND. rad%nev == pw%nev .AND. ABS(rad%r - r) <= 1D-12*MAX(r, 1D0)
   END FUNCTION radial_fits

   PURE SUBROUTINE hankel_pair(hk, am, sgn, k0, f, df)
      ! H_m(k0 r) and dH_m(k0 r)/dr for m = +-am (sgn = (-1)^am for m < 0, else 1), from hk = H_0..H_{am+1} at k0 r
      COMPLEX(8), INTENT(IN) :: hk(0:), sgn
      INTEGER, INTENT(IN) :: am
      REAL(8), INTENT(IN) :: k0
      COMPLEX(8), INTENT(OUT) :: f, df

      f = sgn*hk(am)
      IF (am == 0) THEN
         df = -hk(1)*k0
      ELSE
         df = sgn*0.5D0*(hk(am - 1) - hk(am + 1))*k0
      END IF
   END SUBROUTINE hankel_pair

   PURE SUBROUTINE bessel_k_pair(kv, am, kn, f, df)
      ! K_m(kn r) and dK_m(kn r)/dr, K_{-m} = K_m, from kv = K_0..K_{am+1} at kn r
      REAL(8), INTENT(IN) :: kv(0:), kn
      INTEGER, INTENT(IN) :: am
      COMPLEX(8), INTENT(OUT) :: f, df

      f = kv(am)
      IF (am == 0) THEN
         df = -kv(1)*kn
      ELSE
         df = -0.5D0*(kv(am - 1) + kv(am + 1))*kn
      END IF
   END SUBROUTINE bessel_k_pair

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


   PURE SUBROUTINE pw_bessel_k_scaled(mmax, x, k)
      ! e^x K_m(x), m = 0..mmax, x > 0: pw_bessel_k_scaled_c on the real axis. It stays O(x^-1/2) where K_m itself
      ! underflows (x above about 700), as the outer radial integral of QTFFreeSurface needs far out (pw_bessel_k_far).
      INTEGER, INTENT(IN) :: mmax
      REAL(8), INTENT(IN) :: x
      REAL(8), INTENT(OUT) :: k(0:mmax)
      COMPLEX(8) :: kc(0:mmax)
      INTEGER :: m

      CALL pw_bessel_k_scaled_c(mmax, CMPLX(x, 0D0, KIND=8), kc)
      DO m = 0, mmax
         k(m) = REAL(kc(m))
      END DO
   END SUBROUTINE pw_bessel_k_scaled

   PURE SUBROUTINE pw_bessel_k_scaled_c(mmax, w, k)
      ! e^w K_m(w), m = 0..mmax, for complex w /= 0 with Re w >= 0 (the closed right half-plane, where the Hankel functions
      ! of the tail contours need it: H^(1)_m(z) = (2/pi) i^-(m+1) K_m(-i z)). The algorithm of pw_bessel_k in complex
      ! arithmetic: for |w| <= 2 the power series of K_0 and K_1 (Abramowitz and Stegun 9.6.13 and 9.6.11, principal log),
      ! times e^w; for |w| > 2 Steed's continued fraction CF2 with Temme's normalisation (Temme 1975; Numerical Recipes 2nd
      ! ed. section 6.7), which gives e^w K_0 directly. Then the upward recurrence K_{m+1} = K_{m-1} + (2m/w) K_m, stable
      ! for K (the dominant solution in m), with the common factor e^w. Relative accuracy about 1e-15 (pw_bessel_k_far).
      INTEGER, INTENT(IN) :: mmax
      COMPLEX(8), INTENT(IN) :: w
      COMPLEX(8), INTENT(OUT) :: k(0:mmax)
      REAL(8), PARAMETER :: EULER = 0.57721566490153286061D0, EPS = 1D-17
      COMPLEX(8) :: t, lg, term0, term1, s0, s1, si1, k0, k1, b, d, hh, delh, q1, q2, qnew, q, s, dels, ew
      REAL(8) :: psi1, psi2, cc, a, a1
      INTEGER :: i, m

      IF (ABS(w) <= 2D0) THEN
         t = 0.25D0*w*w
         lg = LOG(0.5D0*w)
         psi1 = -EULER
         psi2 = 1D0 - EULER
         term0 = (1D0, 0D0)
         term1 = (1D0, 0D0)
         s0 = psi1 - lg
         s1 = CMPLX(psi1 + psi2, 0D0, KIND=8)
         si1 = (1D0, 0D0)
         DO i = 1, 200
            psi1 = psi2
            psi2 = psi1 + 1D0/(i + 1)
            term0 = term0*t/(REAL(i, 8)*i)
            term1 = term1*t/(REAL(i, 8)*(i + 1))
            s0 = s0 + term0*(psi1 - lg)
            s1 = s1 + term1*(psi1 + psi2)
            si1 = si1 + term1
            IF (ABS(term0) < EPS*ABS(s0) .AND. ABS(term1) < EPS*ABS(si1)) EXIT
         END DO
         ew = EXP(w)
         k0 = s0*ew
         k1 = (1D0/w + lg*0.5D0*w*si1 - 0.25D0*w*s1)*ew
      ELSE
         b = 2D0*(1D0 + w)
         d = 1D0/b
         hh = d
         delh = d
         q1 = (0D0, 0D0)
         q2 = (1D0, 0D0)
         a1 = 0.25D0
         q = a1
         cc = a1
         a = -a1
         s = 1D0 + q*delh
         DO i = 2, 100000
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
            IF (ABS(dels) < EPS*ABS(s)) EXIT
         END DO
         hh = a1*hh
         k0 = SQRT(PI/(2D0*w))/s
         k1 = k0*(w + 0.5D0 - hh)/w
      END IF

      k(0) = k0
      IF (mmax >= 1) k(1) = k1
      DO m = 1, mmax - 1
         k(m + 1) = k(m - 1) + (2D0*m/w)*k(m)
      END DO
   END SUBROUTINE pw_bessel_k_scaled_c

   PURE SUBROUTINE pw_hankel_scaled(kind, mmax, z, hs)
      ! The Hankel functions with their exponential factor taken out, at complex z with Re z > 0, m = 0..mmax:
      !   kind 1: H^(1)_m(z) e^{-iz},   kind 2: H^(2)_m(z) e^{+iz},
      ! both O(|z|^-1/2) for large |z| in either half-plane, so that a product of them along a rotated contour carries
      ! its exponential e^{i kappa r} explicitly (the tail of QTFFreeSurface). Each kind is computed where it decays and
      ! reflected (H^(2)_m(conj z) = conj H^(1)_m(z) for integer m) into the other half-plane:
      !   H^(1), Im z >= 0: (2/pi) i^-(m+1) e^w K_m(w), w = -iz (Re w >= 0; Abramowitz and Stegun 9.6.4), by
      !     pw_bessel_k_scaled_c;
      !   H^(2), Im z > 0: 2 J_m(z) e^{iz} - H^(1)_m(z) e^{-iz} e^{2iz}, with J_m(z) e^{iz} by Miller's backward recurrence
      !     normalised by e^{-iz} = J_0 + 2 sum_n (-i)^n J_n (the Jacobi-Anger series at theta = pi), whose terms are of the
      !     size of the sum, so that there is no cancellation however large Im z is;
      !   Im z < 0: the conjugates of these at conj z.
      ! Checked against scipy (Amos) to 1e-12 in both half-planes, |z| 0.7 to 3000 (pw_hankel_scaled_values).
      INTEGER, INTENT(IN) :: kind, mmax
      COMPLEX(8), INTENT(IN) :: z
      COMPLEX(8), INTENT(OUT) :: hs(0:mmax)
      COMPLEX(8) :: zu
      LOGICAL :: lower

      lower = AIMAG(z) < 0D0
      zu = z
      IF (lower) zu = CONJG(z)
      ! In the upper half-plane, kind 1 at z is the decaying one; kind 2 at z in the lower half-plane is conj(kind 1 at zu)
      IF ((kind == 1) .NEQV. lower) THEN
         CALL hankel1_upper(mmax, zu, hs)
      ELSE
         CALL hankel2_upper(mmax, zu, hs)
      END IF
      IF (lower) hs = CONJG(hs)
   END SUBROUTINE pw_hankel_scaled

   PURE SUBROUTINE hankel1_upper(mmax, z, hs)
      ! H^(1)_m(z) e^{-iz} for Im z >= 0: (2/pi) i^-(m+1) e^w K_m(w), w = -iz
      INTEGER, INTENT(IN) :: mmax
      COMPLEX(8), INTENT(IN) :: z
      COMPLEX(8), INTENT(OUT) :: hs(0:mmax)
      COMPLEX(8) :: ph
      INTEGER :: m

      CALL pw_bessel_k_scaled_c(mmax, -CI*z, hs)
      ph = -CI*(2D0/PI)
      DO m = 0, mmax
         hs(m) = ph*hs(m)
         ph = -CI*ph
      END DO
   END SUBROUTINE hankel1_upper

   PURE SUBROUTINE hankel2_upper(mmax, z, hs)
      ! H^(2)_m(z) e^{iz} for Im z >= 0: 2 J_m(z) e^{iz} - H^(1)_m(z) e^{-iz} e^{2iz}. J_m(z) e^{iz} by Miller's backward
      ! recurrence J_{n-1} = (2n/z) J_n - J_{n+1} from n0 = mmax + |z| + 40 + 8 |z|^(1/3) (even), normalised by
      ! J_0 + 2 sum_{n>=1} (-i)^n J_n = e^{-iz}; rescaled on the way down so that nothing overflows.
      INTEGER, INTENT(IN) :: mmax
      COMPLEX(8), INTENT(IN) :: z
      COMPLEX(8), INTENT(OUT) :: hs(0:mmax)
      REAL(8), PARAMETER :: BIG = 1D250, SMALL = 1D-250
      COMPLEX(8) :: jp, jc, jm, nrm, ph, e2, h1(0:mmax)
      INTEGER :: n0, n

      n0 = mmax + INT(ABS(z)) + 40 + INT(8D0*ABS(z)**(1D0/3D0))
      IF (MODULO(n0, 2) == 1) n0 = n0 + 1
      jp = (0D0, 0D0)        ! J_{n+1}
      jc = (1D-30, 0D0)      ! J_n, n = n0
      ! (-i)^n0
      ph = (-CI)**MODULO(n0, 4)
      nrm = 2D0*ph*jc
      hs = (0D0, 0D0)
      IF (n0 <= mmax) hs(n0) = jc
      DO n = n0, 1, -1
         jm = (2D0*n/z)*jc - jp
         jp = jc
         jc = jm
         ! jc is J_{n-1}
         ph = CI*ph         ! (-i)^(n-1) = (-i)^n * i
         IF (n - 1 >= 1) THEN
            nrm = nrm + 2D0*ph*jc
         ELSE
            nrm = nrm + jc
         END IF
         IF (n - 1 <= mmax) hs(n - 1) = jc
         IF (ABS(jc) > BIG) THEN
            jc = jc*SMALL
            jp = jp*SMALL
            nrm = nrm*SMALL
            hs = hs*SMALL
         END IF
      END DO
      CALL hankel1_upper(mmax, z, h1)
      e2 = EXP(2D0*CI*z)
      hs = 2D0*hs/nrm - h1*e2
   END SUBROUTINE hankel2_upper

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
