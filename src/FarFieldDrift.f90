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
MODULE FarFieldDrift
   ! Mean wave drift force and yaw moment from the far-field (Kochin) amplitude of the perturbation potential.
   ! Formulae: docs/superpowers/specs/2026-11-04-second-order-waves-design.md, section 5.3. Time factor e^{-i omega t}.
   ! The amplitudes scale with s%amp; the callers pass states with amp = 1 to get loads per unit amplitude.

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, so_options
   USE FirstOrderState, ONLY: depth_factor
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: far_amplitude, group_ratio, far_field_pair, far_field_pairs

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   INTEGER, PARAMETER :: MAX_AZIMUTH = 2**14

CONTAINS

   PURE REAL(8) FUNCTION group_ratio(k, h)
      ! C_g/C_p: 1/2 in deep water (h < 0), otherwise 1/2 (1 + 2kh/sinh 2kh), written to be safe for large kh.
      REAL(8), INTENT(IN) :: k, h
      REAL(8) :: x

      IF (h < 0D0) THEN
         group_ratio = 0.5D0
      ELSE
         x = 2D0*k*h
         IF (x <= 0D0) THEN
            group_ratio = 1D0
         ELSE
            group_ratio = 0.5D0*(1D0 + 2D0*x*EXP(-x)/(1D0 - EXP(-2D0*x)))
         END IF
      END IF
   END FUNCTION group_ratio

   PURE REAL(8) FUNCTION kochin_c(k, nu, h)
      ! Wave-part strength of the Green function: k in deep water, k^2/((k^2-nu^2)h + nu) otherwise.
      REAL(8), INTENT(IN) :: k, nu, h

      IF (h < 0D0) THEN
         kochin_c = k
      ELSE
         kochin_c = k*k/((k*k - nu*nu)*h + nu)
      END IF
   END FUNCTION kochin_c

   SUBROUTINE far_amplitude(mesh, s, theta, a, da)
      ! Far-field amplitude a(theta) of the outgoing perturbation wave, eta_P ~ a e^{ikr}/sqrt(r), and da = da/dtheta.
      ! Computed from phi_p and dphi_p_dn only.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s
      REAL(8), INTENT(IN) :: theta
      COMPLEX(8), INTENT(OUT) :: a, da
      REAL(8) :: c, sn, zf, dzf, xx, yy, zz, nd, dnd
      COMPLEX(8) :: ph, e, dn_e, q, de, ddn_e, hh, dhh, pref
      INTEGER :: i

      c = COS(theta)
      sn = SIN(theta)
      hh = (0D0, 0D0)
      dhh = (0D0, 0D0)
      DO i = 1, mesh%n
         xx = mesh%xc(1, i)
         yy = mesh%xc(2, i)
         zz = mesh%xc(3, i)
         CALL depth_factor(s%k, zz, s%h, zf, dzf)
         ph = EXP(CMPLX(0D0, -s%k*(xx*c + yy*sn), KIND=8))
         nd = mesh%nrm(1, i)*c + mesh%nrm(2, i)*sn
         dnd = -mesh%nrm(1, i)*sn + mesh%nrm(2, i)*c
         e = ph*zf
         dn_e = ph*(-CI*s%k*zf*nd + mesh%nrm(3, i)*dzf)           ! n . grad E
         q = -CI*s%k*(-xx*sn + yy*c)                               ! d(phase)/dtheta
         de = e*q
         ddn_e = dn_e*q + ph*zf*(-CI*s%k)*dnd
         hh = hh + (s%dphi_p_dn(i)*e - s%phi_p(i)*dn_e)*mesh%area(i)
         dhh = dhh + (s%dphi_p_dn(i)*de - s%phi_p(i)*ddn_e)*mesh%area(i)
      END DO
      pref = -(s%omega*kochin_c(s%k, s%nu, s%h)/(2D0*s%g))*SQRT(2D0/(PI*s%k))*EXP(-CI*PI/4D0)
      a = pref*hh
      da = pref*dhh
   END SUBROUTINE far_amplitude

   SUBROUTINE far_field_pair(mesh, s1, s2, opts, f)
      ! Mean load F_12 = (F_x, F_y, M_z) from two states of the same frequency, per unit amplitude pair, with M_z about
      ! the global origin. For s1 = s2 this is the mean drift of one heading. See pair_core for the integration.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      TYPE(so_options), INTENT(IN) :: opts
      COMPLEX(8), INTENT(OUT) :: f(3)
      COMPLEX(8) :: fp(3, 1)
      ! A named array, not the constructor [s1, s2]: an array constructor of a type with allocatable components
      ! crashes at run time under ifx 2024.2 with -heap-arrays 0.
      TYPE(first_order_state) :: pair(2)

      pair(1) = s1
      pair(2) = s2
      CALL pair_core(mesh, pair, [1], [2], opts, fp)
      f = fp(:, 1)
   END SUBROUTINE far_field_pair

   SUBROUTINE far_field_pairs(mesh, states, opts, f)
      ! Mean loads of all pairs of states (one per heading, same frequency): f(:, i, j) = F_ij. Each state's a(theta)
      ! grid is built once and shared by all pairs, and the grids are refined together until every pair converges.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: states(:)
      TYPE(so_options), INTENT(IN) :: opts
      COMPLEX(8), INTENT(OUT) :: f(3, SIZE(states), SIZE(states))
      INTEGER :: ns, i, j, p
      INTEGER, ALLOCATABLE :: ip(:), jp(:)
      COMPLEX(8), ALLOCATABLE :: fp(:, :)

      ns = SIZE(states)
      ALLOCATE (ip(ns*ns), jp(ns*ns), fp(3, ns*ns))
      p = 0
      DO j = 1, ns
         DO i = 1, ns
            p = p + 1
            ip(p) = i
            jp(p) = j
         END DO
      END DO
      CALL pair_core(mesh, states, ip, jp, opts, fp)
      f = RESHAPE(fp, [3, ns, ns])
   END SUBROUTINE far_field_pairs

   SUBROUTINE pair_core(mesh, states, ip, jp, opts, f)
      ! Loads f(:, p) = F_{ip(p), jp(p)} for a list of pairs of states. The azimuth integrals use the trapezoid rule,
      ! doubling the number of points from opts%n_azimuth until, for every pair, the change is below
      ! opts%far_tolerance (absolute after scaling by rho g, or relative), or 2^14 points are reached. a(theta) and
      ! da/dtheta are evaluated once per state and point.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: states(:)
      INTEGER, INTENT(IN) :: ip(:), jp(:)
      TYPE(so_options), INTENT(IN) :: opts
      COMPLEX(8), INTENT(OUT) :: f(:, :)
      COMPLEX(8), ALLOCATABLE :: a(:, :), da(:, :), t(:, :), dt(:, :), pt(:, :), fold(:, :), iv(:, :)
      COMPLEX(8) :: ep, a1b2, da1b2, a2b1, da2b1
      REAL(8) :: kappa, sq, k, scale
      INTEGER :: n, j, m, p, np, ns

      ns = SIZE(states)
      np = SIZE(ip)
      k = states(1)%k
      kappa = 0.5D0*states(1)%rho*states(1)%g*group_ratio(k, states(1)%h)
      sq = SQRT(2D0*PI/k)
      ep = EXP(CI*PI/4D0)
      scale = states(1)%rho*states(1)%g

      ! Terms evaluated at the incident headings
      ALLOCATE (pt(3, np), fold(3, np), iv(3, np))
      DO p = 1, np
         CALL far_amplitude(mesh, states(jp(p)), states(ip(p))%beta, a2b1, da2b1)
         CALL far_amplitude(mesh, states(ip(p)), states(jp(p))%beta, a1b2, da1b2)
         pt(1:2, p) = -sq*([COS(states(ip(p))%beta), SIN(states(ip(p))%beta)]*CONJG(a2b1*ep) &
                           + [COS(states(jp(p))%beta), SIN(states(jp(p))%beta)]*a1b2*ep)
         pt(3, p) = sq*(CI*CONJG(da2b1*ep) - CI*da1b2*ep)
      END DO

      n = MAX(opts%n_azimuth, 4)
      ALLOCATE (a(n, ns), da(n, ns))
      DO m = 1, ns
         DO j = 1, n
            CALL far_amplitude(mesh, states(m), 2D0*PI*REAL(j - 1, 8)/REAL(n, 8), a(j, m), da(j, m))
         END DO
      END DO
      CALL integrate()
      f = assemble()

      DO
         IF (2*n > MAX_AZIMUTH) THEN
            WRITE (*, '(a,i0,a)') 'far_field_pair: azimuth integral not converged with ', n, ' points'
            EXIT
         END IF
         ! Double the grid: old points sit at the odd positions, only the even ones are new
         ALLOCATE (t(2*n, ns), dt(2*n, ns))
         DO m = 1, ns
            DO j = 1, n
               t(2*j - 1, m) = a(j, m)
               dt(2*j - 1, m) = da(j, m)
               CALL far_amplitude(mesh, states(m), 2D0*PI*REAL(2*j - 1, 8)/REAL(2*n, 8), t(2*j, m), dt(2*j, m))
            END DO
         END DO
         CALL MOVE_ALLOC(t, a)
         CALL MOVE_ALLOC(dt, da)
         n = 2*n
         fold = f
         CALL integrate()
         f = assemble()
         IF (ALL(ABS(f - fold) <= opts%far_tolerance*MAX(scale, ABS(f)))) EXIT
      END DO

   CONTAINS

      SUBROUTINE integrate()
         ! Trapezoid rule on [0, 2pi): int a_i a_j* e_theta and int (a_j* da_i - a_i da_j*)/(2i), for every pair
         INTEGER :: pp, jj
         REAL(8) :: th
         COMPLEX(8) :: x1, x2, dx1, dx2

         iv = (0D0, 0D0)
         DO pp = 1, np
            DO jj = 1, n
               th = 2D0*PI*REAL(jj - 1, 8)/REAL(n, 8)
               x1 = a(jj, ip(pp))
               dx1 = da(jj, ip(pp))
               x2 = a(jj, jp(pp))
               dx2 = da(jj, jp(pp))
               iv(1, pp) = iv(1, pp) + x1*CONJG(x2)*COS(th)
               iv(2, pp) = iv(2, pp) + x1*CONJG(x2)*SIN(th)
               iv(3, pp) = iv(3, pp) + (CONJG(x2)*dx1 - x1*CONJG(dx2))/(2D0*CI)
            END DO
         END DO
         iv = iv*(2D0*PI/REAL(n, 8))
      END SUBROUTINE integrate

      FUNCTION assemble() RESULT(r)
         COMPLEX(8) :: r(3, np)

         r(1:2, :) = kappa*(pt(1:2, :) - iv(1:2, :))
         r(3, :) = -(kappa/k)*(iv(3, :) + pt(3, :))
      END FUNCTION assemble

   END SUBROUTINE pair_core

END MODULE FarFieldDrift
