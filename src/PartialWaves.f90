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
   PUBLIC :: pw_green_modes, pw_green_mmax, pw_source_modes, pw_from_panels

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
            ! A zero amplitude adds nothing: skipping it keeps out radial functions that overflow at high order (the
            ! orders that pw_green_modes sets to zero), which would give 0 * Inf = NaN
            IF (pw%a(m, n) == (0D0, 0D0)) CYCLE
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
         IF (pw%a(m, 0) == (0D0, 0D0)) CYCLE       ! as in pw_eval_point
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
            IF (pw%a(m, n) == (0D0, 0D0)) CYCLE    ! as in pw_eval_point
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

   SUBROUTINE pw_green_modes(x, centre, radius, h, k0, kev, mmax, nev, pw, mcut)
      ! The partial waves of the finite-depth Green function G(x, xi) of HAMS as a field of xi outside the circle r = radius
      ! about centre, for a source x inside it (r_x < radius): pw_eval_point(pw, xi) is the full G (1/r and 1/r1 included)
      ! for r_xi >= radius, -h <= xi_z <= 0 (Phase 7b, the outer contraction of the free-surface forcing).
      !
      ! G as FINGREEN3D sums it (EIGEN; John 1950), rho the horizontal distance of x and xi:
      !   G = i P0 H_0(k0 rho) + sum_n Pn K_0(kn rho),
      !   P0 = (pi / N0) cosh k0(z_x+h) cosh k0(xi_z+h),   N0 = (h/2) (1 + sinh 2k0h / (2k0h))   (WAVE_DEPTH_FACTORS),
      !   Pn = 4 / (h (1 + sin 2knh / (2knh))) cos kn(z_x+h) cos kn(xi_z+h).
      ! The series holds for rho > 0 and contains the Rankine parts 1/r + 1/r1 (EIGEN subtracts them afterwards). On the
      ! Z(0) = 1 eigenfunctions of pw_set, P0 = c0 Z_0(z_x) Z_0(xi_z) and Pn = cn Z_n(z_x) Z_n(xi_z) with
      !   c0 = pi cosh^2 k0h / N0 = 2 pi k0 / (k0h sech^2 k0h + tanh k0h),   cn = 4 cos^2 knh / (h + sin 2knh / (2kn)).
      ! Graf's addition theorem for r_xi > r_x, alpha = theta_xi - theta_x (Abramowitz and Stegun 9.1.79 for H_0; for K_0 the
      ! same form with I_m, A&S 9.1.80 at nu = 0 with z -> iz, all terms of sign +, as the test green_modes_match_fingreen
      ! confirms):
      !   H_0(k rho) = sum_m J_m(k r_x) H_m(k r_xi) e^{i m alpha},   K_0(k rho) = sum_m I_m(k r_x) K_m(k r_xi) e^{i m alpha},
      ! so a(m, 0) = i c0 Z_0(z_x) J_m(k0 r_x) e^{-i m theta_x},  a(m, n) = cn Z_n(z_x) I_m(kn r_x) e^{-i m theta_x}
      ! (J_{-m} = (-1)^m J_m, I_{-m} = I_m). I_m is computed as e^{-x} I_m(x) (bessel_i_scaled) and multiplied by e^x, x =
      ! kn r_x < kn R <= PW_KR_MAX, which does not overflow. The terms decay like (r_x / R)^m (pw_green_mmax).
      !
      ! Orders m whose radial function at R, H_{m+1}(k0 R) or K_{m+1}(kn R), exceeds 1e300 (the plain pw_set amplitude
      ! times it would overflow) are set to exactly zero; pw_eval_point and pw_eval_surface skip zero amplitudes. Their
      ! terms are of order (r_x / R)^m / m relative to the leading ones: for k0 R = 0.22 (the test case) the cut is near
      ! m = 118 for H, where (0.9)^m / m is 4e-8. A long wave cuts earlier: k0 R = 0.157 at m = 109, where at r_x = 0.97 R
      ! the dropped tail is about 3e-4 (green_modes_cut_reported). mcut, if present, returns the highest order represented
      ! in every mode (all a(m, n) with |m| > mcut may be zero; mcut = mmax if nothing is cut), so that a caller can
      ! compare it with the order it needs (pw_green_mmax) and stop or warn: the expansion is then only good to about
      ! (r_x / R)^(mcut+1) / (mcut+1).
      REAL(8), INTENT(IN) :: x(3), centre(2), radius, h, k0, kev(:)
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(pw_set), INTENT(OUT) :: pw
      INTEGER, INTENT(OUT), OPTIONAL :: mcut
      REAL(8), PARAMETER :: BIG = 1D300
      REAL(8) :: dxy(2), rx, thx, e1, e2, c0, cn, xn, jm(0:mmax), is(0:mmax), kv(0:mmax + 1), ym, ym1, yp, sgn
      COMPLEX(8) :: e(0:mmax)
      INTEGER :: m, n, m0, mrep

      dxy = x(1:2) - centre
      rx = NORM2(dxy)
      thx = ATAN2(dxy(2), dxy(1))
      IF (rx >= radius) ERROR STOP "pw_green_modes: the source is not inside the circle"
      IF (SIZE(kev) < nev) ERROR STOP "pw_green_modes: fewer evanescent wave numbers than nev"
      IF (nev > 0) THEN
         IF (kev(nev)*radius > KR_MAX) ERROR STOP "pw_green_modes: kn R too large for K_m(kn R) (lower nev or the radius)"
      END IF

      pw%mmax = mmax
      pw%nev = nev
      pw%k0 = k0
      pw%h = h
      pw%radius = radius
      pw%centre = centre
      ALLOCATE (pw%kev(nev), pw%a(-mmax:mmax, 0:nev))
      pw%kev = kev(1:nev)
      pw%a = (0D0, 0D0)
      DO m = 0, mmax
         e(m) = EXP(-CI*(m*thx))
      END DO

      ! Propagating mode. c0 with sech and tanh from e^{-k0h}, which do not overflow at large k0h
      e1 = EXP(-k0*h)
      e2 = e1*e1
      c0 = 2D0*PI*k0/(k0*h*(2D0*e1/(1D0 + e2))**2 + (1D0 - e2)/(1D0 + e2))*z0_mode(k0, h, x(3))
      jm = BESSEL_JN(0, mmax, k0*rx)
      ! The cut: the largest m with |Y_{m'}(k0 R)| <= BIG for m' <= m + 1, Y by the upward recurrence (A&S 9.1.27), stable
      ! for Y; |J| < 1 and |Y| grows with m beyond k0 R, so |H_m| <= BIG then
      ym = BESSEL_Y0(k0*radius)
      ym1 = BESSEL_Y1(k0*radius)
      m0 = -1
      DO m = 0, mmax
         IF (ABS(ym) > BIG .OR. ABS(ym1) > BIG) EXIT
         m0 = m
         yp = (2D0*(m + 1)/(k0*radius))*ym1 - ym
         ym = ym1
         ym1 = yp
      END DO
      mrep = m0
      DO m = 0, m0
         pw%a(m, 0) = CI*c0*jm(m)*e(m)
         sgn = MERGE(-1D0, 1D0, MODULO(m, 2) == 1)
         IF (m > 0) pw%a(-m, 0) = CI*c0*sgn*jm(m)*CONJG(e(m))
      END DO

      ! Evanescent modes: cn Z_n(z_x) = 4 cos knh cos kn(z_x+h) / (h + sin 2knh / (2kn))
      DO n = 1, nev
         cn = 4D0*COS(kev(n)*h)*COS(kev(n)*(x(3) + h))/(h + SIN(2D0*kev(n)*h)/(2D0*kev(n)))
         xn = kev(n)*rx
         CALL bessel_i_scaled(mmax, xn, is)
         CALL pw_bessel_k(mmax + 1, kev(n)*radius, kv)
         DO m = 0, mmax
            IF (.NOT. (kv(m + 1) <= BIG)) THEN
               mrep = MIN(mrep, m - 1)
               EXIT
            END IF
            pw%a(m, n) = cn*is(m)*EXP(xn)*e(m)
            IF (m > 0) pw%a(-m, n) = cn*is(m)*EXP(xn)*CONJG(e(m))
         END DO
      END DO
      IF (PRESENT(mcut)) mcut = mrep
   END SUBROUTINE pw_green_modes

   ! ------------------------------------------------------------------ the circle modes from the hull integral (item 3)
   !
   ! Derivation (own, from the textbook formulas cited). The centroid-route field of a state outside the hull is the sum
   ! over the hull panels j (images included) of (1 / 4 pi) A_j [G(x, xi_j) q_j - n_j . grad_xi G(x, xi_j) phi_j]
   ! (field_influence_row with near_panel = .FALSE., field_potential of WaterlineField). At finite depth G is the series
   ! of pw_green_modes, G = i c0 Z_0(z) Z_0(zeta) H_0(k0 d) + sum_n cn Z_n(z) Z_n(zeta) K_0(kn d), d the horizontal
   ! distance of x and xi (John 1950). Let the source xi be at (rho', theta') about the circle centre, the field point at
   ! (r, theta), with rho' < r.
   !   Propagating terms: Graf's addition theorem (Abramowitz and Stegun 9.1.79 with nu = 0; Watson 1944, section 11.3)
   !     H_0(k d) = sum_m H_m(k r) J_m(k rho') e^{i m (theta - theta')},   rho' < r.
   !   Evanescent terms: put k = i kappa into the same sum (it holds for complex k by analytic continuation, the
   !   condition being |rho'| < |r|). With I_m(x) = i^-m J_m(i x) (A&S 9.6.3) and K_m(x) = (pi i / 2) i^m H_m(i x)
   !   (A&S 9.6.4), each term H_m(i kappa r) J_m(i kappa rho') = (2 / (pi i)) i^-m K_m(kappa r) i^m I_m(kappa rho')
   !   = (2 / (pi i)) K_m(kappa r) I_m(kappa rho'), and the left side is (2 / (pi i)) K_0(kappa d). The powers of i
   !   cancel, so there is no (-1)^m:
   !     K_0(kappa d) = sum_m K_m(kappa r) I_m(kappa rho') e^{i m (theta - theta')},   rho' < r.
   !   (Check at theta = theta': d = r - rho', and every term is positive, so the sum exceeds K_0(kappa r) I_0, as
   !   K_0(r - rho') does; at theta - theta' = pi the terms alternate, d = r + rho'. graf_matches_direct_sum pins both
   !   forms and shows that the opposite signs miss.)
   ! So G(x, xi) = sum_m e^{i m theta} [b(m, 0) H_m(k0 r) Z_0(z) + sum_n b(m, n) K_m(kn r) Z_n(z)], the pw_set form, with
   !   b(m, 0) = i c0 Z_0(zeta) J_m(k0 rho') e^{-i m theta'},   b(m, n) = cn Z_n(zeta) I_m(kn rho') e^{-i m theta'}.
   ! The panel sum makes a(m, n) = (1 / 4 pi) sum_j A_j [q_j b_j(m, n) - phi_j n_j . grad_xi b_j(m, n)].
   ! The gradient in the source point: with f_m = J_m(k rho') e^{-i m theta'} (the conjugate of the regular wave
   ! J_m e^{i m theta'}, J real) and (d_x + i d_y) [J_m e^{i m theta}] = -k J_{m+1} e^{i (m+1) theta},
   ! (d_x - i d_y) [J_m e^{i m theta}] = k J_{m-1} e^{i (m-1) theta} (from A&S 9.1.27),
   !   d_x f_m = (k / 2) (f_{m-1} - f_{m+1}),   d_y f_m = -(i k / 2) (f_{m-1} + f_{m+1}),
   ! and with g_m = I_m(kappa rho') e^{-i m theta'} and the same rules for I with +kappa I_{m+1}, +kappa I_{m-1}
   ! (A&S 9.6.26),
   !   d_x g_m = (kappa / 2) (g_{m-1} + g_{m+1}),   d_y g_m = (i kappa / 2) (g_{m+1} - g_{m-1}),
   ! for every integer m (J_{-m} = (-1)^m J_m, I_{-m} = I_m), so nothing is singular at rho' = 0. The vertical derivative
   ! is that of the depth factor: Z_0' = k0 sinh k0(z+h) / cosh k0h, d/dz cos kn(z+h) = -kn sin kn(z+h).

   PURE SUBROUTINE pw_source_modes(xs, nrm, centre, h, k0, kev, mmax, nev, b, db)
      ! The partial waves about centre of the finite-depth Green function G(xs, .) of a source xs inside the circle (b), and
      ! of n . grad_xs G(xs, .) (db), in the pw_set form (derivation above): pw_eval_point of a pw_set with a = b gives G
      ! at r > |xs - centre|, Rankine parts included, as pw_green_modes does. kev(1:nev) the evanescent wave numbers.
      ! No cut here: pw_from_panels zeroes the orders whose radial functions at R overflow.
      REAL(8), INTENT(IN) :: xs(3), nrm(3), centre(2), h, k0, kev(:)
      INTEGER, INTENT(IN) :: mmax, nev
      COMPLEX(8), INTENT(OUT) :: b(-mmax:mmax, 0:nev), db(-mmax:mmax, 0:nev)
      INTEGER :: n

      CALL source_mode(0, xs, nrm, centre, h, k0, kev, mmax, b(:, 0), db(:, 0))
      DO n = 1, nev
         CALL source_mode(n, xs, nrm, centre, h, k0, kev, mmax, b(:, n), db(:, n))
      END DO
   END SUBROUTINE pw_source_modes

   PURE SUBROUTINE source_mode(n, xs, nrm, centre, h, k0, kev, mmax, b, db)
      ! Mode n of pw_source_modes: n = 0 the propagating mode (k0), n >= 1 the evanescent mode of kev(n). b(-mmax:mmax), db.
      ! I_m(kn rho') is e^{-x} I_m(x) times e^x (bessel_i_scaled), x = kn rho' < kn R <= PW_KR_MAX: no overflow.
      INTEGER, INTENT(IN) :: n, mmax
      REAL(8), INTENT(IN) :: xs(3), nrm(3), centre(2), h, k0, kev(:)
      COMPLEX(8), INTENT(OUT) :: b(-mmax:mmax), db(-mmax:mmax)
      COMPLEX(8) :: e(-mmax - 1:mmax + 1), f(-mmax - 1:mmax + 1), c
      REAL(8) :: dxy(2), rs, ths, e1, e2, zf, dzf, x, kn, cn, v(0:mmax + 1)
      INTEGER :: m

      dxy = xs(1:2) - centre
      rs = NORM2(dxy)
      ths = 0D0
      IF (rs > 0D0) ths = ATAN2(dxy(2), dxy(1))
      DO m = 0, mmax + 1
         e(m) = EXP(-CI*(m*ths))
         e(-m) = CONJG(e(m))
      END DO

      IF (n == 0) THEN
         ! i c0 Z_0(zeta) J_m(k0 rho') e^{-i m theta'}: c0 as in pw_green_modes (sech and tanh from e^{-k0h}); Z_0 and
         ! Z_0' without overflow at large k0 h
         e1 = EXP(-k0*h)
         e2 = e1*e1
         c = CI*2D0*PI*k0/(k0*h*(2D0*e1/(1D0 + e2))**2 + (1D0 - e2)/(1D0 + e2))
         zf = z0_mode(k0, h, xs(3))
         dzf = k0*(EXP(k0*xs(3)) - EXP(-k0*(xs(3) + 2D0*h)))/(1D0 + e2)
         x = k0*rs
         IF (x > 0D0) THEN
            v = BESSEL_JN(0, mmax + 1, x)
         ELSE
            v = 0D0
            v(0) = 1D0
         END IF
         ! f_m = J_m e^{-i m theta'}, J_{-m} = (-1)^m J_m
         DO m = 0, mmax + 1
            f(m) = v(m)*e(m)
            f(-m) = MERGE(-1D0, 1D0, MODULO(m, 2) == 1)*v(m)*e(-m)
         END DO
         DO m = -mmax, mmax
            b(m) = c*zf*f(m)
            db(m) = c*(zf*0.5D0*k0*(nrm(1)*(f(m - 1) - f(m + 1)) - CI*nrm(2)*(f(m - 1) + f(m + 1))) + dzf*nrm(3)*f(m))
         END DO
      ELSE
         ! cn Z_n(zeta) I_m(kn rho') e^{-i m theta'}, cn Z_n(zeta) = 4 cos knh cos kn(zeta+h) / (h + sin 2knh / (2kn)), as
         ! in pw_green_modes
         kn = kev(n)
         cn = 4D0*COS(kn*h)/(h + SIN(2D0*kn*h)/(2D0*kn))
         zf = cn*COS(kn*(xs(3) + h))
         dzf = -cn*kn*SIN(kn*(xs(3) + h))
         x = kn*rs
         CALL bessel_i_scaled(mmax + 1, x, v)
         v = v*EXP(x)
         ! g_m = I_m e^{-i m theta'}, I_{-m} = I_m
         DO m = 0, mmax + 1
            f(m) = v(m)*e(m)
            f(-m) = v(m)*e(-m)
         END DO
         DO m = -mmax, mmax
            b(m) = zf*f(m)
            db(m) = zf*0.5D0*kn*(nrm(1)*(f(m - 1) + f(m + 1)) + CI*nrm(2)*(f(m + 1) - f(m - 1))) + dzf*nrm(3)*f(m)
         END DO
      END IF
   END SUBROUTINE source_mode

   SUBROUTINE pw_from_panels(xs, nrm, area, q, phi, centre, radius, h, k0, kev, mmax, nev, pw, mcut)
      ! The partial waves outside the circle r = radius about centre of the fields of the panel sources xs(3, np) (images
      ! already placed: the caller mirrors the base panels, positions and normals), normals nrm(3, np) and areas area(np),
      ! one field per column of q(np, ns) = d phi_P / dn and phi(np, ns) = phi_P at the panels:
      !   phi_P(x) = (1 / 4 pi) sum_j A_j [G(x, xs_j) q_j - n_j . grad_xs G(x, xs_j) phi_j],
      ! the centroid route of fs_sample_phi (centroid_only) at finite depth, so that the result is the pw_set its samples
      ! fit, without the samples (item 3 of the field-pass speed-ups; derivation above pw_source_modes). k0 and
      ! kev(1:nev) the wave numbers of the frequency (pw_evanescent_from_wvn). Every source must lie inside the circle
      ! (horizontally), else the run stops. Orders whose radial function at R overflows (H_{m+1}(k0 R) or K_{m+1}(kn R)
      ! above 1e300, as in pw_green_modes) are set to zero; mcut, if present, is the highest order kept in every mode.
      ! Each vertical mode n is summed by one thread over the panels in their order: the result does not depend on the
      ! number of threads.
      USE HAMS_mod, ONLY: NTHREAD
      REAL(8), INTENT(IN) :: xs(:, :), nrm(:, :), area(:), centre(2), radius, h, k0, kev(:)
      COMPLEX(8), INTENT(IN) :: q(:, :), phi(:, :)
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(pw_set), INTENT(OUT) :: pw(:)
      INTEGER, INTENT(OUT), OPTIONAL :: mcut
      REAL(8), PARAMETER :: BIG = 1D300
      COMPLEX(8), ALLOCATABLE :: acc(:, :, :), b(:), db(:)
      REAL(8) :: kv(0:mmax + 1), ym, ym1, yp, rmax
      INTEGER :: np, ns, j, n, m, ist, jmax, mtop(0:nev), mrep
      CHARACTER(LEN=200) :: msg

      np = SIZE(xs, 2)
      ns = SIZE(q, 2)
      IF (SIZE(nrm, 2) /= np .OR. SIZE(area) /= np .OR. SIZE(q, 1) /= np .OR. SIZE(phi, 1) /= np .OR. &
          SIZE(phi, 2) /= ns .OR. SIZE(pw) /= ns) ERROR STOP "pw_from_panels: the panel and state arrays do not match"
      IF (.NOT. h > 0D0) ERROR STOP "pw_from_panels: finite depth only"
      IF (SIZE(kev) < nev) ERROR STOP "pw_from_panels: fewer evanescent wave numbers than nev"
      IF (nev > 0) THEN
         IF (kev(nev)*radius > KR_MAX) ERROR STOP "pw_from_panels: kn R too large for K_m(kn R) (lower nev or the radius)"
      END IF
      ! Graf's sums hold for a source inside the circle only; the mirror images of a half or quarter mesh included
      rmax = -1D0
      jmax = 0
      DO j = 1, np
         IF (NORM2(xs(1:2, j) - centre) > rmax) THEN
            rmax = NORM2(xs(1:2, j) - centre)
            jmax = j
         END IF
      END DO
      IF (.NOT. rmax < radius) THEN
         WRITE (msg, '(a,i0,a,3f12.4,a,f0.4,a,f0.4,a)') ' pw_from_panels: panel source ', jmax, ' (images included) at', &
            xs(:, jmax), ' lies ', rmax, ' m from the centre, not inside the circle of radius ', radius, ' m.'
         WRITE (*, '(a)') TRIM(msg)
         ERROR STOP 1
      END IF

      ! The cut of pw_green_modes: the highest m per mode with its radial functions at R up to m + 1 below BIG
      ym = BESSEL_Y0(k0*radius)
      ym1 = BESSEL_Y1(k0*radius)
      mtop(0) = -1
      DO m = 0, mmax
         IF (ABS(ym) > BIG .OR. ABS(ym1) > BIG) EXIT
         mtop(0) = m
         yp = (2D0*(m + 1)/(k0*radius))*ym1 - ym
         ym = ym1
         ym1 = yp
      END DO
      DO n = 1, nev
         CALL pw_bessel_k(mmax + 1, kev(n)*radius, kv)
         mtop(n) = mmax
         DO m = 0, mmax
            IF (.NOT. (kv(m + 1) <= BIG)) THEN
               mtop(n) = m - 1
               EXIT
            END IF
         END DO
      END DO
      mrep = MINVAL(mtop)

      ALLOCATE (acc(-mmax:mmax, ns, 0:nev))
      !$OMP PARALLEL NUM_THREADS(NTHREAD) PRIVATE(b, db, j, ist, m)
      ALLOCATE (b(-mmax:mmax), db(-mmax:mmax))
      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO n = 0, nev
         acc(:, :, n) = (0D0, 0D0)
         DO j = 1, np
            CALL source_mode(n, xs(:, j), nrm(:, j), centre, h, k0, kev, mmax, b, db)
            DO ist = 1, ns
               acc(:, ist, n) = acc(:, ist, n) + area(j)*(q(j, ist)*b - phi(j, ist)*db)
            END DO
         END DO
         DO m = -mmax, mmax
            IF (ABS(m) > mtop(n)) acc(m, :, n) = (0D0, 0D0)
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (b, db)
      !$OMP END PARALLEL

      DO ist = 1, ns
         pw(ist)%mmax = mmax
         pw(ist)%nev = nev
         pw(ist)%k0 = k0
         pw(ist)%h = h
         pw(ist)%radius = radius
         pw(ist)%centre = centre
         ALLOCATE (pw(ist)%kev(nev), pw(ist)%a(-mmax:mmax, 0:nev))
         pw(ist)%kev = kev(1:nev)
         DO n = 0, nev
            pw(ist)%a(:, n) = acc(:, ist, n)/(4D0*PI)
         END DO
      END DO
      IF (PRESENT(mcut)) mcut = mrep
   END SUBROUTINE pw_from_panels

   PURE INTEGER FUNCTION pw_green_mmax(r_x, radius, tol)
      ! The angular order for pw_green_modes: the smallest m with (r_x / radius)^m <= tol, at least 4. The m-th terms of
      ! Graf's sums at r >= radius, J_m(k r_x) H_m(k r) and I_m(k r_x) K_m(k r), are of order (r_x / r)^m / m for large m
      ! (A&S 9.3.1). Not capped: it grows like log(tol) / log(r_x / radius) as the source nears the circle.
      REAL(8), INTENT(IN) :: r_x, radius, tol
      REAL(8) :: q
      INTEGER :: m

      q = r_x/radius
      IF (q >= 1D0) ERROR STOP "pw_green_mmax: r_x is not inside the circle"
      IF (q <= 0D0 .OR. tol >= 1D0) THEN
         pw_green_mmax = 4
         RETURN
      END IF
      m = MAX(1, CEILING(LOG(tol)/LOG(q)))
      ! Rounding of the logarithms: step to the smallest m
      DO WHILE (m > 1)
         IF (q**(m - 1) > tol) EXIT
         m = m - 1
      END DO
      DO WHILE (q**m > tol)
         m = m + 1
      END DO
      pw_green_mmax = MAX(4, m)
   END FUNCTION pw_green_mmax

   PURE SUBROUTINE bessel_i_scaled(mmax, x, is)
      ! e^{-x} I_m(x), m = 0..mmax, x >= 0, by Miller's backward recurrence I_{m-1} = I_{m+1} + (2m/x) I_m (A&S 9.6.26),
      ! for which I_m is the minimal solution, normalised by e^{-x} (I_0 + 2 sum_{m>=1} I_m) = 1 (A&S 9.6.34 at theta = 0).
      ! The start n0 lies beyond mmax and beyond the orders that carry the sum (I_m / I_0 ~ e^{-m^2 / 2x} for large x);
      ! the values are rescaled on the way down so that nothing overflows.
      INTEGER, INTENT(IN) :: mmax
      REAL(8), INTENT(IN) :: x
      REAL(8), INTENT(OUT) :: is(0:mmax)
      REAL(8), PARAMETER :: BIG = 1D250, SMALL = 1D-250
      REAL(8) :: ip, ic, im, nrm
      INTEGER :: n0, n

      is = 0D0
      IF (x < 1D-300) THEN
         is(0) = 1D0
         RETURN
      END IF
      n0 = mmax + 40 + INT(SQRT(80D0*x)) + INT(x)
      ip = 0D0          ! I_{n+1}
      ic = 1D-30        ! I_n, n = n0
      nrm = 2D0*ic
      DO n = n0, 1, -1
         im = (2D0*n/x)*ic + ip
         ip = ic
         ic = im
         ! ic is I_{n-1}
         IF (n - 1 >= 1) THEN
            nrm = nrm + 2D0*ic
         ELSE
            nrm = nrm + ic
         END IF
         IF (n - 1 <= mmax) is(n - 1) = ic
         IF (ic > BIG) THEN
            ic = ic*SMALL
            ip = ip*SMALL
            nrm = nrm*SMALL
            is = is*SMALL
         END IF
      END DO
      is = is/nrm
   END SUBROUTINE bessel_i_scaled

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
