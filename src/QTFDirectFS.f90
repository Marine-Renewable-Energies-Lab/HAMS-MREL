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
MODULE QTFDirectFS
   ! Free-surface part of the right side of the direct second-order solve (Phase 7b): the forcing Q_D of phi_D on z = 0
   ! enters the boundary integral equation of QTFDirect as
   !    r_FS(x) = (s_D / g) int_{z=0} G(x, xi) Q_D(xi) dS_xi
   ! at every row point x of the solve (the full hull panels and, with IRSP = 1, the full lid points).
   !
   ! The sign s_D (DIRECT_FS_SIGN). Green's second identity over the fluid volume, bounded by the hull S, the free surface
   ! F (z = 0 outside the waterlines), the bottom and a far cylinder, with n the normal out of the fluid on all of them and
   ! lap G = -4 pi delta (G = 1/r + regular), gives at x on a smooth part of S (half the delta, a principal value)
   !    -2 pi phi(x) = int_S (phi dG/dn - G dphi/dn) dS + int_F (phi dG/dn - G dphi/dn) dS
   ! (the bottom term vanishes, both have d/dz = 0 there; the far cylinder by the radiation condition, the limit the outer
   ! integral of QTFFreeSurface takes on its rotated tail contours), that is
   !    2 pi phi + int_S phi dG/dn dS = int_S G dphi/dn dS - int_F (phi dG/dn - G dphi/dn) dS,
   ! the equation of QTFDirect (ASSB_LEFT, ASSB_RBC: 2 pi phi + int phi dG/dn = int G dphi/dn) plus the F term. On F the
   ! normal out of the fluid is +z. G satisfies -Omega^2 G + g G_z = 0 there (V = Omega^2 / g), and phi_D satisfies
   ! -Omega^2 phi_D + g phi_D,z = Q_D (QTFFreeSurface header), so
   !    phi_D G_z - G phi_D,z = phi_D (Omega^2 / g) G - G (Omega^2 phi_D + Q_D) / g = -G Q_D / g,
   ! and -int_F (phi_D dG/dn - G dphi_D/dn) dS = +(1/g) int_F G Q_D dS. So s_D = +1: the free-surface part is added to the
   ! right side with the sign of the body part int_S G b dS. At a lid point (inside the body, outside the fluid) the
   ! identity has 0 in place of -2 pi phi, and the F term is the same: the lid rows get the same r_FS. For the first-order
   ! problem Q = 0 and the F term vanishes, which is why ASSB_RBC has none. Pinned by T_dir = T_pot (full): test
   ! direct_fs_equals_indirect_full, 2.8e-3 of max |T_pot| with s_D = +1, 0.78 with -1.
   ! g is the gravity of the Green function, g = Omega^2 / V (Const_mod G): Q_D must be formed with the same g (fs%g).
   !
   ! The inner part (the inner free-surface mesh of build_fs_inner, nfs panels, areas w_j, centroids xi_j):
   !    r_in(x_i) = sum_j gfs(i, j) w_j Q_D(xi_j),   gfs(i, j) = (1 / w_j) int_{panel j} G(x_i, xi) dS
   ! with Q_D constant per panel (the one-point rule of the inner integral of QTFFreeSurface, fs_inner_pair). The panel
   ! integral of G, as CALGREEN does it on the hull: a panel within NEAR_PANEL_SIZES (50) panel sizes of x_i takes the
   ! Rankine part 1/r + 1/r1 in closed form (rankine_moments of QTFDirect at x_i and at its mirror (x, y, -z): 1/r1 is
   ! 1/r of the mirrored point), and the rest G - 1/r - 1/r1 (FINGREEN3D with flag 1, which leaves out 1/r, less 1/r1) by
   ! an n x n Gauss-Legendre rule on the panel (fs_remainder_order; n = 1 is the centroid); a panel further away takes G
   ! at its centroid (FINGREEN3D, flag 0). The rest is not smooth on the scale of the panel when x_i is close: next to
   ! the waterline it carries the logarithm of the free-surface wave part and in finite depth the bottom images, so one
   ! point leaves errors of the order of (a / d)^2 of it (a the panel radius, d the distance): 1.7e-3 of r_in on the
   ! coarse default inner mesh of the 448-panel cylinder, away from the waterline as at it (test fs_rhs_inner_only with
   ! one point). With n from d / a (fs_remainder_order) it is 2.9e-7, about the accuracy of FINGREEN3D itself.
   !
   ! The outer part (beyond the inner mesh: the annulus, the tail and the sliver of fs_outer_moments): G(x_i, .) on z = 0
   ! outside the circle R is a sum of partial waves (pw_green_modes, Graf's addition theorem), so
   !    r_out(x_i) = sum_{m, n} a(m, n, i) mom(m, n)
   ! with mom = int B_mn Q_D dS from fs_outer_moments and the amplitudes in the form of the contraction there (the
   ! set_amplitudes form of QTFFreeSurface): B_m0 = H_m(k0 r) e^{i m theta} with a(m, 0, i) = a_pw(m, 0), and
   ! B_mn = K_m(kn r) / K_m(kn R) e^{i m theta} with a(m, n, i) = a_pw(m, n) K_m(kn R), a_pw the amplitudes of
   ! pw_green_modes(x_i). mom(m, n) pairs with a(m, n), same index (mom uses the coefficient q_{+m} of
   ! e^{-i m theta} of Q_D, so that sum a mom is the integral; test outer_moments_contract_green).
   !
   ! Orders. The m-th Graf term at r >= R is at most of order (r_x / R)^m / m (A&S 9.3.1; I_m(a) / I_m(b) <= (a / b)^m
   ! for a < b, as I_m(x) / x^m increases), so mg = pw_green_mmax(r_max, R, DIRECT_FS_TOL), r_max the largest horizontal
   ! distance of a row point from the centre. The evanescent term n is of order e^{-kn (R - r_x)} (I_m(kn r_x) K_m(kn R),
   ! A&S 9.7.1 and 9.7.2), and kn > (n - 1/2) pi / h (the roots of -Omega^2 / g = kn tan(kn h)), so the terms beyond n
   ! sum to at most e^{-kn d} / (1 - e^{-pi d / h}), d = R - r_max. nev is the smallest n with that below DIRECT_FS_TOL:
   !    nev = ceiling(L h / (pi d) + 1/2),   L = ln(1 / (DIRECT_FS_TOL (1 - e^{-pi d / h})))
   ! It grows like 1 / d as the rows near the circle (a point at 0.9 R needed nev = 200). The rows
   ! nearest the circle (the waterline rows of a wall hull, or the lid points next to the waterline) set both. nev is
   ! capped where kn R exceeds PW_KR_MAX (pw_green_modes stops beyond), with a warning of the error left.
   ! pw_green_modes zeroes the orders whose H_{m+1}(k0 R) or K_{m+1}(kn R) overflow (long waves on small circles) and
   ! returns the highest order kept, mcut (the same for all rows: it depends on k0, kn and R only). fs_outer_moments is
   ! not finite at those orders, so only |m| <= mcut is contracted (never Inf x 0). If mcut < mg, the outer part drops
   ! terms m > mcut, each at most q^m / m, q = r_max / R, so at most sum_{m > mcut} q^m / m <= q^(mcut + 1) /
   ! ((mcut + 1) (1 - q)) (geometric bound) of the leading terms (fs%cut_error). Above DIRECT_FS_TOL a warning states
   ! that bound and the run goes on: the loss is bounded and stated, and a stop at one low Omega would
   ! end a long run. The bound is not small for every input: R close to the hull (q near 1) with a long wave (small
   ! k0 R, low mcut) raises it (a measured case gives 3.4e-3 for R = 1.05 r_x, Omega = 0.01 rad/s in 15 m of
   ! water); the warning names the cure, a larger qtf_fs_radius.
   ! The sliver of fs_outer_moments lies between R cos(pi / np) and R, where the Graf series still converges for
   ! r_max < R cos(pi / np) at the slower rate (r_max / (R cos(pi / np)))^m; it covers a small fraction of the free
   ! surface, so the order is kept at that of R.
   !
   ! r_far. The outer moments drop the evanescent B_mn beyond r_far, where they are e^{-kev(1) (r_far - R)} of their value
   ! at R. direct_fs_rfar raises r_far so that kev(1) (r_far - R) >= RFAR_EV_DECAY (40): the dropped part is below e^-40.
   !
   ! Symmetry. The rows are the full points: full panel (a - 1) nbase + i is image a of base panel i, full lid row
   ! (a - 1) nl + l is base lid point l mirrored by image a, as direct_rhs_b2_block has them; the inner mesh is the whole
   ! free surface. The block right side of block ip at base row i is sum_a RXY(ip, a) r(image a of i) (RXY_MULTI for
   ! several bodies), the layout of ASSB_RHS_FIELD_IRR.
   !
   ! The Green function is that of the solver globals (V, WVN, NK, H; finite depth only, as Phase 6a): the caller sets
   ! them to Omega, as direct_at_omega requires.

   USE Const_mod, ONLY: RXY, RXY_MULTI
   USE HAMS_mod, ONLY: NTHREAD
   USE WaveDyn_mod, ONLY: H, NK, V, WVN, W1
   USE PanelMesh_mod, ONLY: ISX, ISY, iXYZ_GLOBAL_MULTI_COMB_P
   USE Inerfs_mod, ONLY: iXYZ_P
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE SecondOrderTypes, ONLY: panel_mesh, cs_points
   USE PartialWaves, ONLY: pw_set, pw_green_modes, pw_green_mmax, pw_bessel_k, PW_KR_MAX
   ! direct_fs_data is defined in QTFDirect, so that direct_at_omega can take it, and passed on from here
   USE QTFDirect, ONLY: direct_lid_points, rankine_moments, direct_fs_data
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: direct_fs_data, direct_fs_setup, direct_fs_rhs, direct_fs_rfar, fs_remainder_order
   PUBLIC :: direct_fs_orders, direct_fs_array_mb, direct_fs_estimate
   PUBLIC :: DIRECT_FS_SIGN, DIRECT_FS_TOL

   ! The sign s_D of r_FS = s_D / g int G Q_D dS (header)
   INTEGER, PARAMETER :: DIRECT_FS_SIGN = 1
   ! The relative size of the dropped Graf terms that sets mg and nev (header)
   REAL(8), PARAMETER :: DIRECT_FS_TOL = 1D-8
   ! As CALGREEN and QTFDirect: a panel closer than this many panel sizes gets the Rankine part in closed form
   REAL(8), PARAMETER :: NEAR_PANEL_SIZES = 50D0
   ! The rest G - 1/r - 1/r1 on a near panel: Gauss-Legendre order for an error of about REM_EPS (fs_remainder_order),
   ! at most REM_MAX points per direction
   REAL(8), PARAMETER :: REM_EPS = 1D-8
   INTEGER, PARAMETER :: REM_MAX = 8
   ! direct_fs_rfar: the decay of the slowest evanescent basis function at r_far
   REAL(8), PARAMETER :: RFAR_EV_DECAY = 40D0
   REAL(8), PARAMETER :: PI = 3.14159265358979323846D0

CONTAINS

   SUBROUTINE direct_fs_setup(mesh, fs_pts, fs_panels, centre, radius, h_in, omega, is_multi, fs)
      ! The inner Green matrix and the Graf amplitudes of every row at Omega (header). fs_pts and fs_panels are the inner
      ! free surface of build_fs_inner (points and panels, in the same order), centre and radius its circle, h_in the depth
      ! (finite). The solver globals must be at omega (W1, V, WVN, H); with IRSP = 1 the lid points of direct_lid_points
      ! give the lid rows. Stops if the globals are not at omega or h_in, in deep water, if fs_pts and fs_panels differ,
      ! or if a row point lies outside the circle. Writes a warning if the Graf orders or the evanescent modes had to be
      ! cut (header). The pointers of fs to the 6a store are left unset (the caller of direct_at_omega sets them).
      TYPE(panel_mesh), INTENT(IN) :: mesh, fs_panels
      TYPE(cs_points), INTENT(IN) :: fs_pts
      REAL(8), INTENT(IN) :: centre(2), radius, h_in, omega
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(direct_fs_data), INTENT(OUT) :: fs
      REAL(8), ALLOCATABLE :: xr(:, :), xv(:, :, :), prad(:), wvn_l(:), kr(:, :)
      INTEGER, ALLOCATABLE :: nv(:)
      REAL(8) :: d, rr
      INTEGER :: i, j, k, nfs, nrow, mcut_all, nev_rule, n, ierr
      CHARACTER(LEN=300) :: line
      EXTERNAL :: DISPERSION

      IF (.NOT. (H > 0D0) .OR. ABS(H - h_in) > 1D-9*h_in) THEN
         WRITE (*, '(a)') ' direct_fs_setup: the free-surface part needs finite depth, and the depth of the solver.'
         ERROR STOP 1
      END IF
      IF (.NOT. (omega > 0D0) .OR. ABS(W1 - omega) > 1D-12*omega) THEN
         WRITE (*, '(a,ES14.6,a,ES14.6,a)') ' direct_fs_setup: Omega =', omega, ' rad/s, but the solver is at', W1, &
            ' rad/s; the Green function must be that of Omega.'
         ERROR STOP 1
      END IF
      nfs = fs_pts%n
      IF (fs_panels%n /= nfs) THEN
         WRITE (*, '(a)') ' direct_fs_setup: the inner points and the inner panels do not match.'
         ERROR STOP 1
      END IF

      fs%omega = omega
      fs%g = omega*omega/V
      fs%h = h_in
      fs%radius = radius
      fs%centre = centre
      fs%multi = is_multi
      fs%nbase = mesh%nbase
      fs%nsys = mesh%nsys
      fs%n = mesh%n
      fs%nl = direct_lid_points(is_multi)
      fs%nfs = nfs
      nrow = mesh%n + mesh%nsys*fs%nl
      fs%nrow = nrow

      ! The row points
      ALLOCATE (xr(3, nrow))
      DO i = 1, mesh%n
         xr(:, i) = mesh%xc(:, i)
      END DO
      DO i = 1, mesh%nsys*fs%nl
         xr(:, mesh%n + i) = lid_point(i, fs%nl, is_multi)
      END DO
      fs%r_max = 0D0
      DO i = 1, nrow
         fs%r_max = MAX(fs%r_max, HYPOT(xr(1, i) - centre(1), xr(2, i) - centre(2)))
      END DO
      IF (.NOT. fs%r_max < radius) THEN
         WRITE (line, '(a,f0.4,a,f0.4,a)') ' direct_fs_setup: a row point lies ', fs%r_max, &
            ' m from the centre, not inside the circle R = ', radius, ' m.'
         WRITE (*, '(a)') TRIM(line)
         ERROR STOP 1
      END IF

      ! The inner panels: vertices, and the radius about the centroid (for the order of the rest)
      ALLOCATE (xv(3, 4, nfs), nv(nfs), prad(nfs))
      DO j = 1, nfs
         nv(j) = fs_panels%ncn(j)
         prad(j) = 0D0
         DO k = 1, nv(j)
            xv(:, k, j) = fs_panels%xnode(:, fs_panels%ncon(k, j))
            prad(j) = MAX(prad(j), NORM2(xv(:, k, j) - fs_pts%x(:, j)))
         END DO
         IF (nv(j) == 3) xv(:, 4, j) = xv(:, 3, j)
      END DO
      fs%w = fs_pts%w

      ! The inner Green matrix, OpenMP over the rows
      ALLOCATE (fs%gfs(nrow, nfs), STAT=ierr)
      IF (ierr /= 0) CALL alloc_stop('the inner Green matrix', 16D0*REAL(nrow, 8)*REAL(nfs, 8)/1D6)
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) SCHEDULE(DYNAMIC, 4) PRIVATE(j)
      DO i = 1, nrow
         DO j = 1, nfs
            fs%gfs(i, j) = panel_green(xr(:, i), xv(:, :, j), nv(j), fs_panels%nrm(:, j), fs_pts%x(:, j), fs_pts%w(j), &
                                       fs_pts%psize(j), prad(j))
         END DO
      END DO
      !$OMP END PARALLEL DO

      ! The orders (header)
      d = radius - fs%r_max
      CALL direct_fs_orders(fs%r_max, radius, h_in, fs%mg, nev_rule)
      ALLOCATE (wvn_l(nev_rule + 1))
      CALL DISPERSION(wvn_l, nev_rule + 1, omega, h_in)
      fs%k0 = wvn_l(1)
      n = nev_rule
      DO WHILE (n > 0)
         IF (wvn_l(n + 1)*radius <= PW_KR_MAX) EXIT
         n = n - 1
      END DO
      fs%nev = n
      fs%kev = wvn_l(2:n + 1)
      IF (n < nev_rule) THEN
         rr = 0D0
         IF (n > 0) rr = EXP(-fs%kev(n)*d)
         IF (n == 0) rr = 1D0
         WRITE (line, '(a,ES12.4,a,i0,a,i0,a,ES9.2,a)') ' QTF direct FS: warning: at Omega =', omega, &
            ' rad/s the rows nearest the circle need ', nev_rule, ' evanescent modes, but kn R stays below the limit '// &
            'of the partial waves only up to ', n, '; the outer part is good to about ', rr, &
            ' there. A larger qtf_fs_radius lowers it.'
         WRITE (*, '(a)') TRIM(line)
      END IF

      ! K_m(kn R) of the contraction form
      ALLOCATE (kr(0:fs%mg, MAX(fs%nev, 1)))
      DO n = 1, fs%nev
         CALL pw_bessel_k(fs%mg, fs%kev(n)*radius, kr(:, n))
      END DO

      ALLOCATE (fs%a(-fs%mg:fs%mg, 0:fs%nev, nrow), STAT=ierr)
      IF (ierr /= 0) CALL alloc_stop('the Graf amplitudes', direct_fs_array_mb(nrow, 0, fs%mg, fs%nev))
      mcut_all = fs%mg
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) SCHEDULE(DYNAMIC, 4) REDUCTION(MIN:mcut_all)
      DO i = 1, nrow
         CALL row_amplitudes(i, mcut_all)
      END DO
      !$OMP END PARALLEL DO
      fs%mcut = mcut_all
      ! Orders above the smallest cut of all rows are zeroed in every row (one cut for the contraction)
      fs%a(-fs%mg:-fs%mcut - 1, :, :) = (0D0, 0D0)
      fs%a(fs%mcut + 1:fs%mg, :, :) = (0D0, 0D0)
      fs%cut_error = 0D0
      IF (fs%mcut < fs%mg) THEN
         ! The bound of the dropped tail (header): the first dropped term times 1 / (1 - q)
         fs%cut_error = (fs%r_max/radius)**(fs%mcut + 1)/(REAL(fs%mcut + 1, 8)*(1D0 - fs%r_max/radius))
         IF (fs%cut_error > DIRECT_FS_TOL) THEN
            WRITE (line, '(a,ES12.4,a,i0,a,i0,a,ES9.2,a)') ' QTF direct FS: warning: at Omega =', omega, &
               ' rad/s the partial waves of G keep the orders up to ', fs%mcut, ' (H_m(k0 R) or K_m(kn R) overflow '// &
               'beyond), fewer than the ', fs%mg, ' the rows nearest the circle need; the outer part is good to about ', &
               fs%cut_error, ' there. A larger qtf_fs_radius lowers it.'
            WRITE (*, '(a)') TRIM(line)
         END IF
      END IF

   CONTAINS

      SUBROUTINE alloc_stop(what, mb)
         CHARACTER(LEN=*), INTENT(IN) :: what
         REAL(8), INTENT(IN) :: mb

         WRITE (*, '(a,a,a,F0.1,a,i0,a,i0,a,F0.3,a)') ' direct_fs_setup: not enough memory for ', what, ' (', mb, ' MB; ', &
            nrow, ' rows, orders ', fs%mg, ', R =', radius, ' m). A larger qtf_fs_radius lowers the orders and the '// &
            'evanescent modes; so does a coarser hull mesh.'
         ERROR STOP 1
      END SUBROUTINE alloc_stop

      SUBROUTINE row_amplitudes(row, mc)
         ! a(:, :, row) in the contraction form, and the cut of pw_green_modes into mc
         INTEGER, INTENT(IN) :: row
         INTEGER, INTENT(INOUT) :: mc
         TYPE(pw_set) :: pw
         INTEGER :: m, nn, mcut

         CALL pw_green_modes(xr(:, row), centre, radius, h_in, fs%k0, fs%kev, fs%mg, fs%nev, pw, mcut)
         mc = MIN(mc, mcut)
         fs%a(:, :, row) = (0D0, 0D0)
         DO m = -mcut, mcut
            fs%a(m, 0, row) = pw%a(m, 0)
         END DO
         DO nn = 1, fs%nev
            DO m = -mcut, mcut
               fs%a(m, nn, row) = pw%a(m, nn)*kr(ABS(m), nn)
            END DO
         END DO
      END SUBROUTINE row_amplitudes

   END SUBROUTINE direct_fs_setup

   SUBROUTINE direct_fs_rhs(fs, q, mom, rhs)
      ! Adds the free-surface part r_FS = s_D / g [sum_j gfs(i, j) w_j q_j + sum_mn a(m, n, i) mom(m, n)] of ncol right
      ! sides to rhs:
      !   q    (nfs, ncol)              Q_D at the inner points (fs_forcing_points), formed with the g of fs (fs%g)
      !   mom  (-mg:mg, 0:nev, ncol)    the outer moments of fs_outer_moments with k0, kev, mg and nev of fs
      !   rhs  (nbase + nl, ncol, nsys) the layout of ASSB_RHS_FIELD_IRR (header); the block form is added
      ! Only |m| <= mcut of mom is used; a moment that is not finite there stops the run.
      TYPE(direct_fs_data), INTENT(IN) :: fs
      COMPLEX(8), INTENT(IN) :: q(:, :), mom(-fs%mg:, 0:, :)
      COMPLEX(8), INTENT(INOUT) :: rhs(:, :, :)
      COMPLEX(8), ALLOCATABLE :: wq(:, :), mw(:, :, :), r(:, :)
      INTEGER :: ncol, nb, nl, ns, nmode, ip, im, i, col, j
      REAL(8) :: rs
      EXTERNAL :: ZGEMM

      ncol = SIZE(q, 2)
      nb = fs%nbase
      nl = fs%nl
      ns = fs%nsys
      IF (SIZE(q, 1) /= fs%nfs .OR. SIZE(mom, 1) /= 2*fs%mg + 1 .OR. SIZE(mom, 2) /= fs%nev + 1 .OR. &
          SIZE(mom, 3) /= ncol .OR. SIZE(rhs, 1) /= nb + nl .OR. SIZE(rhs, 2) /= ncol .OR. SIZE(rhs, 3) /= ns) THEN
         WRITE (*, '(a)') ' direct_fs_rhs: q, mom or rhs does not have the size of the free-surface data.'
         ERROR STOP 1
      END IF
      IF (ncol == 0) RETURN
      nmode = (2*fs%mg + 1)*(fs%nev + 1)

      ! The orders the amplitudes represent, the rest zero (they may not be finite there)
      ALLOCATE (mw(-fs%mg:fs%mg, 0:fs%nev, ncol))
      mw = (0D0, 0D0)
      mw(-fs%mcut:fs%mcut, :, :) = mom(-fs%mcut:fs%mcut, :, :)
      IF (.NOT. ALL(IEEE_IS_FINITE(REAL(mw))) .OR. .NOT. ALL(IEEE_IS_FINITE(AIMAG(mw)))) THEN
         WRITE (*, '(a,ES12.4,a,i0,a)') ' direct_fs_rhs: at Omega =', fs%omega, ' rad/s an outer moment of order |m| <= ', &
            fs%mcut, ' is not finite.'
         ERROR STOP 1
      END IF
      ALLOCATE (wq(fs%nfs, ncol), r(fs%nrow, ncol))
      DO col = 1, ncol
         DO j = 1, fs%nfs
            wq(j, col) = fs%w(j)*q(j, col)
         END DO
      END DO
      r = (0D0, 0D0)
      IF (fs%nfs > 0) CALL ZGEMM('N', 'N', fs%nrow, ncol, fs%nfs, (1D0, 0D0), fs%gfs, fs%nrow, wq, fs%nfs, (0D0, 0D0), r, &
                                 fs%nrow)
      CALL ZGEMM('T', 'N', fs%nrow, ncol, nmode, (1D0, 0D0), fs%a, nmode, mw, nmode, (1D0, 0D0), r, fs%nrow)
      r = (DIRECT_FS_SIGN/fs%g)*r

      ! The block form, as ASSB_RBC (ASSB_RBC_IRR)
      DO ip = 1, ns
         DO im = 1, ns
            rs = RXY(ip, im)
            IF (fs%multi) rs = RXY_MULTI(ip, im)
            DO i = 1, nb
               rhs(i, :, ip) = rhs(i, :, ip) + rs*r((im - 1)*nb + i, :)
            END DO
            DO i = 1, nl
               rhs(nb + i, :, ip) = rhs(nb + i, :, ip) + rs*r(fs%n + (im - 1)*nl + i, :)
            END DO
         END DO
      END DO
   END SUBROUTINE direct_fs_rhs

   SUBROUTINE direct_fs_orders(r_max, radius, h_in, mg, nev)
      ! The Graf order M and the number of evanescent modes of the direct free-surface set-up for the row points up to
      ! r_max from the centre and the circle R (header): M from pw_green_mmax, N_ev = ceil(L h / (pi d) + 1/2) with
      ! d = R - r_max and L = ln(1 / (tol (1 - e^(-pi d / h)))). Neither depends on Omega; the cut at kn R > PW_KR_MAX
      ! of direct_fs_setup can only lower N_ev, so this is the upper bound.
      REAL(8), INTENT(IN) :: r_max, radius, h_in
      INTEGER, INTENT(OUT) :: mg, nev
      REAL(8) :: d, el

      d = radius - r_max
      mg = pw_green_mmax(r_max, radius, DIRECT_FS_TOL)
      el = LOG(1D0/(DIRECT_FS_TOL*(1D0 - EXP(-PI*d/h_in))))
      nev = MAX(1, CEILING(el*h_in/(PI*d) + 0.5D0))
   END SUBROUTINE direct_fs_orders

   PURE REAL(8) FUNCTION direct_fs_array_mb(nrow, nfs, mg, nev)
      ! The size in MB of the per-Omega arrays: gfs(nrow, nfs) and a(-mg:mg, 0:nev, nrow), both complex(8)
      INTEGER, INTENT(IN) :: nrow, nfs, mg, nev

      direct_fs_array_mb = 16D0*REAL(nrow, 8)*(REAL(nfs, 8) + REAL(2*mg + 1, 8)*REAL(nev + 1, 8))/1D6
   END FUNCTION direct_fs_array_mb

   SUBROUTINE direct_fs_estimate(mesh, nfs, centre, radius, h_in, is_multi, mg, nev, mb)
      ! The orders and the size of the per-Omega arrays of direct_fs_setup before the frequency loop (they do not depend
      ! on Omega): the rows are the hull panel centres and the lid points, as in direct_fs_setup.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: nfs
      REAL(8), INTENT(IN) :: centre(2), radius, h_in
      LOGICAL, INTENT(IN) :: is_multi
      INTEGER, INTENT(OUT) :: mg, nev
      REAL(8), INTENT(OUT) :: mb
      REAL(8) :: r_max, xl(3)
      INTEGER :: i, nl

      nl = direct_lid_points(is_multi)
      r_max = 0D0
      DO i = 1, mesh%n
         r_max = MAX(r_max, HYPOT(mesh%xc(1, i) - centre(1), mesh%xc(2, i) - centre(2)))
      END DO
      DO i = 1, mesh%nsys*nl
         xl = lid_point(i, nl, is_multi)
         r_max = MAX(r_max, HYPOT(xl(1) - centre(1), xl(2) - centre(2)))
      END DO
      IF (.NOT. r_max < radius) THEN
         WRITE (*, '(a,f0.4,a,f0.4,a)') ' direct_fs_estimate: a row point lies ', r_max, &
            ' m from the centre, not inside the circle R = ', radius, ' m.'
         ERROR STOP 1
      END IF
      CALL direct_fs_orders(r_max, radius, h_in, mg, nev)
      mb = direct_fs_array_mb(mesh%n + mesh%nsys*nl, nfs, mg, nev)
   END SUBROUTINE direct_fs_estimate

   REAL(8) FUNCTION direct_fs_rfar(fs, r_far)
      ! r_far raised, if needed, so that the slowest evanescent basis function has decayed by e^-40 there (header):
      ! max(r_far, R + RFAR_EV_DECAY / kev(1)). Without evanescent modes r_far itself.
      TYPE(direct_fs_data), INTENT(IN) :: fs
      REAL(8), INTENT(IN) :: r_far

      direct_fs_rfar = r_far
      IF (fs%nev > 0) direct_fs_rfar = MAX(r_far, fs%radius + RFAR_EV_DECAY/fs%kev(1))
   END FUNCTION direct_fs_rfar

   PURE INTEGER FUNCTION fs_remainder_order(d, a)
      ! The Gauss-Legendre order n (n x n points) for the rest G - 1/r - 1/r1 on a near panel of radius a about its
      ! centroid at the distance d of the row point. The rest is singular where 1/r1 is, at the mirrored row point, at the
      ! same distance d from a panel on z = 0, so the rule converges as for a function analytic in an ellipse of
      ! semi-axis sum rho = t + sqrt(t^2 + 1), t = d / a: error ~ rho^(-2n) (Davis and Rabinowitz, Methods of Numerical
      ! Integration, 2nd ed., 4.6). n = ceiling(ln(1 / REM_EPS) / (2 ln rho)), from 1 to REM_MAX.
      REAL(8), INTENT(IN) :: d, a
      REAL(8) :: t, rho

      IF (.NOT. (a > 0D0)) THEN
         fs_remainder_order = 1
         RETURN
      END IF
      t = d/a
      rho = t + SQRT(t*t + 1D0)
      fs_remainder_order = MIN(REM_MAX, MAX(1, CEILING(LOG(1D0/REM_EPS)/(2D0*LOG(rho)))))
   END FUNCTION fs_remainder_order

   COMPLEX(8) FUNCTION panel_green(x, xv, nv, nrm, xc, w, psz, prad)
      ! (1 / w) int_panel G(x, xi) dS over the inner panel with vertices xv(:, 1:nv) (xv(:, 4) = xv(:, 3) for a
      ! triangle), normal nrm, centroid xc, area w, size psz and radius prad about xc (header)
      REAL(8), INTENT(IN) :: x(3), xv(3, 4), nrm(3), xc(3), w, psz, prad
      INTEGER, INTENT(IN) :: nv
      REAL(8) :: d, xm(3), m0, m0i, g0(3), m1(3, 3), xq(3), wq, jac, r1
      REAL(8) :: gx(REM_MAX), gw(REM_MAX)
      COMPLEX(8) :: grn(4), rest
      INTEGER :: n, iu, iv

      d = NORM2(x - xc)
      IF (d > NEAR_PANEL_SIZES*psz) THEN
         CALL FINGREEN3D(xc(1), x(1), xc(2), x(2), xc(3), x(3), V, WVN, NK, H, grn, 0)
         panel_green = grn(1)
         RETURN
      END IF
      ! Rankine part in closed form: 1/r at x, 1/r1 as 1/r at the mirrored point
      xm(1) = x(1)
      xm(2) = x(2)
      xm(3) = -x(3)
      CALL rankine_moments(xv, nv, nrm, x, xc, m0, g0, m1)
      CALL rankine_moments(xv, nv, nrm, xm, xc, m0i, g0, m1)
      ! The rest G - 1/r - 1/r1 (FINGREEN3D with flag 1 leaves out 1/r only)
      n = fs_remainder_order(d, prad)
      rest = (0D0, 0D0)
      IF (n == 1) THEN
         CALL FINGREEN3D(xc(1), x(1), xc(2), x(2), xc(3), x(3), V, WVN, NK, H, grn, 1)
         r1 = SQRT((xc(1) - x(1))**2 + (xc(2) - x(2))**2 + (xc(3) + x(3))**2)
         rest = (grn(1) - 1D0/r1)*w
      ELSE
         CALL gauss_legendre(n, gx, gw)
         DO iu = 1, n
            DO iv = 1, n
               CALL bilinear_point(xv, gx(iu), gx(iv), xq, jac)
               wq = gw(iu)*gw(iv)*jac
               IF (wq <= 0D0) CYCLE
               CALL FINGREEN3D(xq(1), x(1), xq(2), x(2), xq(3), x(3), V, WVN, NK, H, grn, 1)
               r1 = SQRT((xq(1) - x(1))**2 + (xq(2) - x(2))**2 + (xq(3) + x(3))**2)
               rest = rest + (grn(1) - 1D0/r1)*wq
            END DO
         END DO
      END IF
      panel_green = (m0 + m0i + rest)/w
   END FUNCTION panel_green

   PURE SUBROUTINE bilinear_point(xv, u, s, xq, jac)
      ! The point of parameters (u, s) in [-1, 1]^2 of the bilinear map of the quadrilateral xv(:, 1:4) (a triangle with
      ! xv(:, 4) = xv(:, 3)), and the area element jac = |d xq / du x d xq / ds|
      REAL(8), INTENT(IN) :: xv(3, 4), u, s
      REAL(8), INTENT(OUT) :: xq(3), jac
      REAL(8) :: n1, n2, n3, n4, du(3), ds(3), c(3)

      n1 = 0.25D0*(1D0 - u)*(1D0 - s)
      n2 = 0.25D0*(1D0 + u)*(1D0 - s)
      n3 = 0.25D0*(1D0 + u)*(1D0 + s)
      n4 = 0.25D0*(1D0 - u)*(1D0 + s)
      xq = n1*xv(:, 1) + n2*xv(:, 2) + n3*xv(:, 3) + n4*xv(:, 4)
      du = 0.25D0*((1D0 - s)*(xv(:, 2) - xv(:, 1)) + (1D0 + s)*(xv(:, 3) - xv(:, 4)))
      ds = 0.25D0*((1D0 - u)*(xv(:, 4) - xv(:, 1)) + (1D0 + u)*(xv(:, 3) - xv(:, 2)))
      c(1) = du(2)*ds(3) - du(3)*ds(2)
      c(2) = du(3)*ds(1) - du(1)*ds(3)
      c(3) = du(1)*ds(2) - du(2)*ds(1)
      jac = SQRT(c(1)**2 + c(2)**2 + c(3)**2)
   END SUBROUTINE bilinear_point

   PURE SUBROUTINE gauss_legendre(n, x, w)
      ! Nodes x(1:n) and weights w(1:n) of n-point Gauss-Legendre on [-1, 1], by Newton on P_n from the Chebyshev guesses
      ! (A&S 25.4.29)
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(OUT) :: x(:), w(:)
      REAL(8) :: z, p0, p1, p2, dp
      INTEGER :: i, k, it

      DO i = 1, n
         z = COS(PI*(i - 0.25D0)/(n + 0.5D0))
         DO it = 1, 100
            p0 = 1D0
            p1 = z
            DO k = 2, n
               p2 = ((2*k - 1)*z*p1 - (k - 1)*p0)/k
               p0 = p1
               p1 = p2
            END DO
            dp = n*(z*p1 - p0)/(z*z - 1D0)
            z = z - p1/dp
            IF (ABS(p1/dp) < 1D-15) EXIT
         END DO
         x(i) = z
         w(i) = 2D0/((1D0 - z*z)*dp*dp)
      END DO
   END SUBROUTINE gauss_legendre

   FUNCTION lid_point(lrow, nl, multi) RESULT(xl)
      ! Full lid row lrow = (a - 1) nl + l: base lid point l mirrored by image a (as QTFDirect, direct_rhs_b2_block)
      INTEGER, INTENT(IN) :: lrow, nl
      LOGICAL, INTENT(IN) :: multi
      REAL(8) :: xl(3)
      INTEGER :: l, ia

      ia = (lrow - 1)/nl + 1
      l = lrow - (ia - 1)*nl
      IF (multi) THEN
         xl = iXYZ_GLOBAL_MULTI_COMB_P(l, :)
      ELSE
         xl = iXYZ_P(l, :)
      END IF
      IF (ia == 2) THEN
         IF (ISX == 1 .AND. ISY == 0) THEN
            xl(1) = -xl(1)
         ELSE
            xl(2) = -xl(2)
         END IF
      END IF
   END FUNCTION lid_point

END MODULE QTFDirectFS
