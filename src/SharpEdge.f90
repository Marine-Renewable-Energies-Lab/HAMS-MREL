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
MODULE SharpEdge
   ! Sharp edges of the hull and panel weights for integrands that are singular there (qtf_edge_weights).
   !
   ! Corner flow: near a straight edge where the fluid fills a wedge of angle alpha_f, the potential of the flow across the
   ! edge goes like rho^(pi/alpha_f) cos(pi theta/alpha_f), with rho the distance to the edge, so the velocity goes like
   ! rho^p, p = pi/alpha_f - 1. alpha_f = 2 pi minus the interior angle of the body. At the bottom edge of a truncated
   ! cylinder (interior angle 90 deg, alpha_f = 3 pi/2) p = -1/3.
   !
   ! A one-point rule takes the value of rho^p at the panel centre, which misses the singular part on the panels next to
   ! the edge. The product weight of a panel strip of width h across the edge, centre at rho_c >= h/2,
   !
   !    f = int_{rho_c - h/2}^{rho_c + h/2} rho^p d rho / (h rho_c^p)
   !      = [ (rho_c + h/2)^(p+1) - (rho_c - h/2)^(p+1) ] / ( (p+1) h rho_c^p )
   !
   ! makes the rule exact for rho^p times a function that is constant across the strip. On the row that touches the edge
   ! h = 2 rho_c and f = 2^p/(p+1): 1.1906 for p = -1/3. The second row of the same width would take 1.0086; it is not
   ! weighted (on the sharp-cylinder QTF fixture that moved B2 further from the exact reference, not closer).
   ! A panel on two edges (a box corner) takes the product of the two factors, which is exact for rho_1^p rho_2^p.
   !
   ! Edges: mesh edges shared by exactly two panels of one body, after merging coincident nodes (so that a panel and its
   ! symmetry image join across the plane), whose normals differ by more than edge_deg. A waterline segment belongs to one
   ! panel and is not an edge; an edge in a symmetry plane is one if its normals differ by the angle, like any other.
   ! Only convex edges (alpha_f > pi, p < 0) are weighted: at a concave edge the velocity is not singular.
   ! f multiplies the panel area in the hull term of B2 (b2_body_term in QTFPotential).

   USE SecondOrderTypes, ONLY: panel_mesh
   USE SurfaceGeometry, ONLY: body_sizes, merge_nodes, sort_index
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: sharp_edge_weights, edge_exponent, edge_strip_factor, DEFAULT_EDGE_DEG

   REAL(8), PARAMETER :: DEFAULT_EDGE_DEG = 45D0
   REAL(8), PARAMETER :: PI = 3.141592653589793D0

CONTAINS

   PURE REAL(8) FUNCTION edge_exponent(interior)
      ! Exponent p of the velocity rho^p near an edge of interior angle interior (rad, in the body): p = pi/alpha_f - 1,
      ! alpha_f = 2 pi - interior
      REAL(8), INTENT(IN) :: interior

      edge_exponent = PI/(2D0*PI - interior) - 1D0
   END FUNCTION edge_exponent

   PURE REAL(8) FUNCTION edge_strip_factor(rho_c, h, p)
      ! f of a strip of width h whose centre is at rho_c >= h/2 from the edge, for rho^p (header). p > -1.
      REAL(8), INTENT(IN) :: rho_c, h, p

      edge_strip_factor = ((rho_c + 0.5D0*h)**(p + 1D0) - MAX(rho_c - 0.5D0*h, 0D0)**(p + 1D0)) &
                          /((p + 1D0)*h*rho_c**p)
   END FUNCTION edge_strip_factor

   SUBROUTINE sharp_edge_weights(mesh, edge_deg, f, nedge, nseg, nweighted)
      ! Panel weights f(n) of the header for the panels of mesh (images expanded), and per body the number of sharp edges
      ! (chains of edge segments that meet at a node and turn by less than edge_deg there, so the bottom edge of a cylinder
      ! is one and a box has 12), of sharp edge segments (convex or not) and of panels with f /= 1.
      !   edge_deg   normals of the two panels differ by more than this angle (degrees): DEFAULT_EDGE_DEG
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: edge_deg
      REAL(8), INTENT(OUT) :: f(:)
      INTEGER, INTENT(OUT) :: nedge(:), nseg(:), nweighted(:)

      REAL(8), ALLOCATABLE :: bsize(:), key(:), dir(:, :)
      INTEGER, ALLOCATABLE :: rep(:), ea(:), eb(:), ep(:), ord(:), sa(:), sb(:), sp(:, :), root(:), ends(:, :)
      LOGICAL, ALLOCATABLE :: counted(:)
      REAL(8) :: cosmin, c, s, interior, p, xa(3), xb(3), t(3), d(3), rho
      INTEGER :: n, nn, ne, e, i, j, m, mp, g, h, ns, k, q, a, b, node

      n = mesh%n
      nn = mesh%nnode
      f = 1D0
      nedge = 0
      nseg = 0
      nweighted = 0
      IF (n == 0) RETURN
      cosmin = COS(edge_deg*PI/180D0)
      CALL body_sizes(mesh, bsize)
      CALL merge_nodes(mesh, bsize, rep)

      ! Edges of all panels, by their merged end nodes
      ne = SUM(mesh%ncn)
      ALLOCATE (ea(ne), eb(ne), ep(ne), key(ne))
      e = 0
      DO i = 1, n
         DO m = 1, mesh%ncn(i)
            mp = MOD(m, mesh%ncn(i)) + 1
            e = e + 1
            ea(e) = rep(mesh%ncon(m, i))
            eb(e) = rep(mesh%ncon(mp, i))
            ep(e) = i
            IF (ea(e) == eb(e)) THEN
               key(e) = -1D0
            ELSE
               key(e) = REAL(MIN(ea(e), eb(e)), 8)*REAL(nn + 1, 8) + REAL(MAX(ea(e), eb(e)), 8)
            END IF
         END DO
      END DO
      CALL sort_index(key, ord)

      ! Sharp segments: groups of exactly two panels of one body whose normals differ by more than edge_deg
      ALLOCATE (sa(ne/2 + 1), sb(ne/2 + 1), sp(2, ne/2 + 1))
      ns = 0
      g = 1
      DO WHILE (g <= ne)
         h = g
         DO WHILE (h < ne)
            IF (key(ord(h + 1)) /= key(ord(g))) EXIT
            h = h + 1
         END DO
         IF (key(ord(g)) >= 0D0 .AND. h == g + 1) THEN
            i = ep(ord(g))
            j = ep(ord(h))
            IF (i /= j .AND. mesh%body(i) == mesh%body(j)) THEN
               c = DOT_PRODUCT(mesh%nrm(:, i), mesh%nrm(:, j))/ &
                   (NORM2(mesh%nrm(:, i))*NORM2(mesh%nrm(:, j)))
               IF (c < cosmin) THEN
                  ns = ns + 1
                  sa(ns) = ea(ord(g))
                  sb(ns) = eb(ord(g))
                  sp(:, ns) = [i, j]
               END IF
            END IF
         END IF
         g = h + 1
      END DO

      ! Chains of segments: segments that share a node and turn by less than edge_deg there are one edge
      ALLOCATE (root(ns), dir(3, ns), ends(2, nn), counted(ns))
      ends = 0
      DO k = 1, ns
         root(k) = k
         dir(:, k) = mesh%xnode(:, sb(k)) - mesh%xnode(:, sa(k))
         dir(:, k) = dir(:, k)/NORM2(dir(:, k))
      END DO
      DO k = 1, ns
         DO a = 1, 2
            node = MERGE(sa(k), sb(k), a == 1)
            DO b = 1, 2   ! up to two earlier segments are remembered at a node; more meet only at corners
               q = ends(b, node)
               IF (q == 0) THEN
                  ends(b, node) = k
                  EXIT
               END IF
               IF (ABS(DOT_PRODUCT(dir(:, k), dir(:, q))) > cosmin) CALL join(k, q)
            END DO
         END DO
      END DO
      counted = .FALSE.
      DO k = 1, ns
         nseg(mesh%body(sp(1, k))) = nseg(mesh%body(sp(1, k))) + 1
         q = find(k)
         IF (.NOT. counted(q)) THEN
            counted(q) = .TRUE.
            nedge(mesh%body(sp(1, k))) = nedge(mesh%body(sp(1, k))) + 1
         END IF
      END DO

      ! Weights of the panels on a convex edge
      DO k = 1, ns
         i = sp(1, k)
         j = sp(2, k)
         xa = mesh%xnode(:, sa(k))
         xb = mesh%xnode(:, sb(k))
         c = DOT_PRODUCT(mesh%nrm(:, i), mesh%nrm(:, j))/(NORM2(mesh%nrm(:, i))*NORM2(mesh%nrm(:, j)))
         ! Normals point into the body: at a convex edge each panel's centre is on the body side of the other's plane
         s = DOT_PRODUCT(mesh%xc(:, j) - xa, mesh%nrm(:, i)) + DOT_PRODUCT(mesh%xc(:, i) - xa, mesh%nrm(:, j))
         IF (s <= 0D0) CYCLE
         interior = PI - ACOS(MAX(-1D0, MIN(1D0, c)))
         p = edge_exponent(interior)
         t = (xb - xa)/NORM2(xb - xa)
         DO a = 1, 2
            m = sp(a, k)
            d = mesh%xc(:, m) - xa
            rho = NORM2(d - DOT_PRODUCT(d, t)*t)
            IF (rho <= 0D0) CYCLE
            f(m) = f(m)*edge_strip_factor(rho, 2D0*rho, p)
         END DO
      END DO
      DO i = 1, n
         IF (f(i) /= 1D0) nweighted(mesh%body(i)) = nweighted(mesh%body(i)) + 1
      END DO

   CONTAINS

      INTEGER FUNCTION find(k0)
         INTEGER, INTENT(IN) :: k0
         find = k0
         DO WHILE (root(find) /= find)
            find = root(find)
         END DO
      END FUNCTION find

      SUBROUTINE join(k1, k2)
         INTEGER, INTENT(IN) :: k1, k2
         INTEGER :: r1, r2
         r1 = find(k1)
         r2 = find(k2)
         IF (r1 /= r2) root(MAX(r1, r2)) = MIN(r1, r2)
      END SUBROUTINE join

   END SUBROUTINE sharp_edge_weights

END MODULE SharpEdge
