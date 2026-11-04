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
MODULE SurfaceGeometry
   ! Panel neighbours, waterline extraction and per-body waterplane and volume properties for the near-field drift
   ! integral. Every direction comes from the panel normal mesh%nrm, which points out of the fluid, into the body. The
   ! node order of a panel is used only to walk its edges, because symmetry images keep the base node order on mirrored
   ! coordinates and so have reversed winding.

   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: build_surface_geom, DEFAULT_SHARP_DEG

   REAL(8), PARAMETER :: DEFAULT_SHARP_DEG = 45D0
   REAL(8), PARAMETER :: MERGE_TOL = 1D-6   ! node merging tolerance, relative to the body size
   REAL(8), PARAMETER :: NH_MIN = 1D-3      ! smallest horizontal part of a unit normal accepted on a waterline panel
   REAL(8), PARAMETER :: PI = 3.141592653589793D0

CONTAINS

   SUBROUTINE build_surface_geom(mesh, wl_z_tol, sharp_deg, geom, ok, msg)
      ! Builds neighbours, waterline segments and body properties of a full (images expanded) panel mesh.
      !   wl_z_tol   waterline tolerance, relative to the largest bounding-box dimension of each body
      !   sharp_deg  neighbours whose normals differ by more than this angle are dropped (DEFAULT_SHARP_DEG)
      !   ok         .false. only for a stop: a lid panel in z = 0, or a waterline panel with a nearly horizontal normal
      !   msg        explains a stop; for ok = .true. it holds warnings (waterline gaps, waterline not found), else ''
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: wl_z_tol, sharp_deg
      TYPE(surface_geom), INTENT(OUT) :: geom
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg

      INTEGER :: n, nb, nn, i, m, mp, a, b, e, ne, g, h, iw, ib
      REAL(8), ALLOCATABLE :: bsize(:), ztol(:), key(:)
      INTEGER, ALLOCATABLE :: rep(:), ea(:), eb(:), ep(:), ord(:), segn(:,:), cnt(:), nbad(:), nfree(:), nwl_b(:)
      LOGICAL, ALLOCATABLE :: single(:), onwl(:), nonconforming(:)
      REAL(8) :: nhm, t(2), x1(3), x2(3)
      CHARACTER(LEN=300) :: line

      n = mesh%n
      nb = mesh%nbody
      nn = mesh%nnode
      ok = .TRUE.
      msg = ''
      ALLOCATE (geom%awp(nb), geom%xf(nb), geom%yf(nb), geom%l11(nb), geom%l12(nb), geom%l22(nb), geom%vol(nb), &
                geom%xb(nb), geom%yb(nb), geom%zb(nb), geom%n_wl_gaps(nb))
      geom%awp = 0D0; geom%xf = 0D0; geom%yf = 0D0; geom%l11 = 0D0; geom%l12 = 0D0; geom%l22 = 0D0
      geom%vol = 0D0; geom%xb = 0D0; geom%yb = 0D0; geom%zb = 0D0; geom%n_wl_gaps = 0
      geom%nwl = 0

      CALL body_sizes(mesh, bsize)
      ztol = wl_z_tol*bsize

      ! Lid panels (all nodes on the free surface) are not part of the wetted hull.
      DO i = 1, n
         IF (ALL(ABS(mesh%xnode(3, mesh%ncon(1:mesh%ncn(i), i))) <= ztol(mesh%body(i)))) THEN
            WRITE (line, '(a,i0,a,i0,a)') 'panel ', i, ' of body ', mesh%body(i), &
               ' lies in the free surface z = 0; lid panels are not allowed in the hull mesh'
            msg = TRIM(line)
            ok = .FALSE.
            geom = surface_geom()   ! one defined state after a stop: nwl = 0, nothing allocated
            RETURN
         END IF
      END DO

      CALL merge_nodes(mesh, bsize, rep)

      ! Edges of all panels; an edge belonging to exactly one panel (after merging) is a boundary edge.
      ne = SUM(mesh%ncn)
      ALLOCATE (ea(ne), eb(ne), ep(ne), key(ne), single(ne), onwl(ne))
      e = 0
      DO i = 1, n
         DO m = 1, mesh%ncn(i)
            mp = MOD(m, mesh%ncn(i)) + 1
            e = e + 1
            ea(e) = mesh%ncon(m, i)
            eb(e) = mesh%ncon(mp, i)
            ep(e) = i
            a = rep(ea(e))
            b = rep(eb(e))
            IF (a == b) THEN
               key(e) = -1D0
            ELSE
               ! integer-valued and exact in REAL(8) as long as nnode < 9e7
               key(e) = REAL(MIN(a, b), 8)*REAL(nn + 1, 8) + REAL(MAX(a, b), 8)
            END IF
         END DO
      END DO
      CALL sort_index(key, ord)
      single = .FALSE.
      g = 1
      DO WHILE (g <= ne)
         h = g
         DO WHILE (h < ne)
            IF (key(ord(h + 1)) /= key(ord(g))) EXIT
            h = h + 1
         END DO
         IF (key(ord(g)) >= 0D0 .AND. h == g) single(ord(g)) = .TRUE.
         g = h + 1
      END DO
      ALLOCATE (nfree(nb), nwl_b(nb), nonconforming(n))
      nfree = 0
      nwl_b = 0
      nonconforming = .FALSE.
      DO e = 1, ne
         onwl(e) = .FALSE.
         IF (.NOT. single(e)) CYCLE
         ib = mesh%body(ep(e))
         IF (ABS(mesh%xnode(3, ea(e))) <= ztol(ib) .AND. ABS(mesh%xnode(3, eb(e))) <= ztol(ib)) THEN
            onwl(e) = .TRUE.
         ELSE
            nfree(ib) = nfree(ib) + 1
            nonconforming(ep(e)) = .TRUE.
         END IF
      END DO

      ! Waterline segments, in panel order.
      geom%nwl = COUNT(onwl)
      ALLOCATE (geom%wl_x1(3, geom%nwl), geom%wl_x2(3, geom%nwl), geom%wl_xm(3, geom%nwl), geom%wl_len(geom%nwl), &
                geom%wl_nh(3, geom%nwl), geom%wl_nz(geom%nwl), geom%wl_panel(geom%nwl), geom%wl_body(geom%nwl))
      ALLOCATE (segn(2, geom%nwl))
      iw = 0
      DO e = 1, ne
         IF (.NOT. onwl(e)) CYCLE
         iw = iw + 1
         i = ep(e)
         nhm = SQRT(mesh%nrm(1, i)**2 + mesh%nrm(2, i)**2)
         IF (nhm < NH_MIN*SQRT(SUM(mesh%nrm(:, i)**2))) THEN
            WRITE (line, '(a,i0,a,es10.3,a)') 'panel ', i, ' has a waterline edge but a nearly horizontal normal ' &
               //'(horizontal part ', nhm, '); the waterline direction is undefined'
            msg = TRIM(line)
            ok = .FALSE.
            geom = surface_geom()   ! one defined state after a stop: nwl = 0, nothing allocated
            RETURN
         END IF
         x1 = mesh%xnode(:, ea(e))
         x2 = mesh%xnode(:, eb(e))
         geom%wl_nh(:, iw) = [mesh%nrm(1, i)/nhm, mesh%nrm(2, i)/nhm, 0D0]
         geom%wl_nz(iw) = mesh%nrm(3, i)
         ! Walk the waterline counter-clockwise, which keeps the body on the left.
         t = [geom%wl_nh(2, iw), -geom%wl_nh(1, iw)]
         IF ((x2(1) - x1(1))*t(1) + (x2(2) - x1(2))*t(2) < 0D0) THEN
            geom%wl_x1(:, iw) = x2
            geom%wl_x2(:, iw) = x1
            segn(:, iw) = [rep(eb(e)), rep(ea(e))]
         ELSE
            geom%wl_x1(:, iw) = x1
            geom%wl_x2(:, iw) = x2
            segn(:, iw) = [rep(ea(e)), rep(eb(e))]
         END IF
         geom%wl_xm(:, iw) = 0.5D0*(x1 + x2)
         geom%wl_len(iw) = SQRT(SUM((x2 - x1)**2))
         geom%wl_panel(iw) = i
         geom%wl_body(iw) = mesh%body(i)
         nwl_b(mesh%body(i)) = nwl_b(mesh%body(i)) + 1
      END DO

      CALL waterplane_properties(mesh, geom)
      CALL volume_properties(mesh, geom)
      CALL build_neighbours(mesh, rep, nonconforming, sharp_deg, geom%nbr_start, geom%nbr)

      ! A waterline node should close up in exactly two segments of its body.
      ALLOCATE (cnt(nn), nbad(nb))
      cnt = 0
      nbad = 0
      DO iw = 1, geom%nwl
         cnt(segn(:, iw)) = cnt(segn(:, iw)) + 1
      END DO
      DO iw = 1, geom%nwl
         DO m = 1, 2
            IF (cnt(segn(m, iw)) > 0 .AND. cnt(segn(m, iw)) /= 2) THEN
               nbad(geom%wl_body(iw)) = nbad(geom%wl_body(iw)) + 1
               cnt(segn(m, iw)) = -cnt(segn(m, iw))   ! counted
            END IF
         END DO
      END DO
      DO ib = 1, nb
         geom%n_wl_gaps(ib) = (nbad(ib) + 1)/2
         IF (geom%n_wl_gaps(ib) > 0) THEN
            WRITE (line, '(a,i0,a,i0,a)') 'body ', ib, ': the waterline has ', geom%n_wl_gaps(ib), &
               ' gap(s); some waterline nodes are not shared by exactly two segments'
            CALL add_line(msg, line)
         END IF
         IF (nwl_b(ib) == 0 .AND. nfree(ib) > 0) THEN
            WRITE (line, '(a,i0,a,es9.2,a)') 'body ', ib, ': no waterline found within |z| <= ', wl_z_tol, &
               ' x body size, but the mesh has open edges; check that the hull reaches z = 0'
            CALL add_line(msg, line)
         END IF
      END DO
   END SUBROUTINE build_surface_geom

   SUBROUTINE add_line(msg, line)
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: msg
      CHARACTER(LEN=*), INTENT(IN) :: line

      IF (LEN(msg) > 0) msg = msg//NEW_LINE('a')
      msg = msg//TRIM(line)
   END SUBROUTINE add_line

   SUBROUTINE body_sizes(mesh, bsize)
      ! Largest bounding-box dimension of the nodes used by each body.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), ALLOCATABLE, INTENT(OUT) :: bsize(:)
      REAL(8), ALLOCATABLE :: lo(:,:), hi(:,:)
      INTEGER :: i, m, ib

      ALLOCATE (bsize(mesh%nbody), lo(3, mesh%nbody), hi(3, mesh%nbody))
      lo = HUGE(1D0)
      hi = -HUGE(1D0)
      DO i = 1, mesh%n
         ib = mesh%body(i)
         DO m = 1, mesh%ncn(i)
            lo(:, ib) = MIN(lo(:, ib), mesh%xnode(:, mesh%ncon(m, i)))
            hi(:, ib) = MAX(hi(:, ib), mesh%xnode(:, mesh%ncon(m, i)))
         END DO
      END DO
      DO ib = 1, mesh%nbody
         bsize(ib) = MAXVAL(hi(:, ib) - lo(:, ib))
         IF (.NOT. bsize(ib) > 0D0) bsize(ib) = 1D0   ! body without panels
      END DO
   END SUBROUTINE body_sizes

   SUBROUTINE merge_nodes(mesh, bsize, rep)
      ! rep(i) is the lowest node index among the nodes of the same body that lie within MERGE_TOL x body size of node i,
      ! chained transitively. Nodes are sorted by x and swept, so symmetry images join through their shared plane.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: bsize(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: rep(:)
      INTEGER, ALLOCATABLE :: lst(:), seen(:), ord(:)
      REAL(8), ALLOCATABLE :: key(:)
      REAL(8) :: tol
      INTEGER :: ib, i, m, k, a, b, ia, ib2, ra, rb, nl

      ALLOCATE (rep(mesh%nnode), seen(mesh%nnode), lst(mesh%nnode))
      rep = [(i, i = 1, mesh%nnode)]
      seen = 0
      DO ib = 1, mesh%nbody
         tol = MERGE_TOL*bsize(ib)
         nl = 0
         DO i = 1, mesh%n
            IF (mesh%body(i) /= ib) CYCLE
            DO m = 1, mesh%ncn(i)
               k = mesh%ncon(m, i)
               IF (seen(k) /= ib) THEN
                  seen(k) = ib
                  nl = nl + 1
                  lst(nl) = k
               END IF
            END DO
         END DO
         IF (nl == 0) CYCLE
         ALLOCATE (key(nl))
         key = mesh%xnode(1, lst(1:nl))
         CALL sort_index(key, ord)
         DO a = 1, nl
            ia = lst(ord(a))
            DO b = a + 1, nl
               ib2 = lst(ord(b))
               IF (mesh%xnode(1, ib2) - mesh%xnode(1, ia) > tol) EXIT
               IF (MAXVAL(ABS(mesh%xnode(:, ib2) - mesh%xnode(:, ia))) <= tol) THEN
                  ra = find_root(rep, ia)
                  rb = find_root(rep, ib2)
                  IF (ra < rb) THEN
                     rep(rb) = ra
                  ELSE IF (rb < ra) THEN
                     rep(ra) = rb
                  END IF
               END IF
            END DO
         END DO
         DEALLOCATE (key, ord)
      END DO
      DO i = 1, mesh%nnode
         rep(i) = find_root(rep, i)
      END DO
   END SUBROUTINE merge_nodes

   INTEGER FUNCTION find_root(rep, i)
      INTEGER, INTENT(INOUT) :: rep(:)
      INTEGER, INTENT(IN) :: i
      INTEGER :: r

      r = i
      DO WHILE (rep(r) /= r)
         rep(r) = rep(rep(r))   ! path halving
         r = rep(r)
      END DO
      find_root = r
   END FUNCTION find_root

   SUBROUTINE sort_index(key, idx)
      ! idx lists 1..SIZE(key) in ascending order of key (stable bottom-up merge sort).
      REAL(8), INTENT(IN) :: key(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: idx(:)
      INTEGER, ALLOCATABLE :: tmp(:)
      INTEGER :: m, w, lo, mid, hi, i, j, k

      m = SIZE(key)
      ALLOCATE (idx(m), tmp(m))
      idx = [(i, i = 1, m)]
      w = 1
      DO WHILE (w < m)
         DO lo = 1, m, 2*w
            mid = MIN(lo + w - 1, m)
            hi = MIN(lo + 2*w - 1, m)
            i = lo
            j = mid + 1
            k = lo
            DO WHILE (i <= mid .AND. j <= hi)
               IF (key(idx(i)) <= key(idx(j))) THEN
                  tmp(k) = idx(i)
                  i = i + 1
               ELSE
                  tmp(k) = idx(j)
                  j = j + 1
               END IF
               k = k + 1
            END DO
            DO WHILE (i <= mid)
               tmp(k) = idx(i)
               i = i + 1
               k = k + 1
            END DO
            DO WHILE (j <= hi)
               tmp(k) = idx(j)
               j = j + 1
               k = k + 1
            END DO
         END DO
         idx = tmp
         w = 2*w
      END DO
   END SUBROUTINE sort_index

   SUBROUTINE waterplane_properties(mesh, geom)
      ! Green's theorem over the counter-clockwise waterline segments, with coordinates relative to the body's xr.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(INOUT) :: geom
      REAL(8) :: p1(2), p2(2), c, sx(mesh%nbody), sy(mesh%nbody)
      INTEGER :: iw, ib

      sx = 0D0
      sy = 0D0
      DO iw = 1, geom%nwl
         ib = geom%wl_body(iw)
         p1 = geom%wl_x1(1:2, iw) - mesh%xr(1:2, ib)
         p2 = geom%wl_x2(1:2, iw) - mesh%xr(1:2, ib)
         c = p1(1)*p2(2) - p2(1)*p1(2)
         geom%awp(ib) = geom%awp(ib) + 0.5D0*c
         sx(ib) = sx(ib) + (p1(1) + p2(1))*c/6D0
         sy(ib) = sy(ib) + (p1(2) + p2(2))*c/6D0
         geom%l11(ib) = geom%l11(ib) + (p1(1)**2 + p1(1)*p2(1) + p2(1)**2)*c/12D0
         geom%l22(ib) = geom%l22(ib) + (p1(2)**2 + p1(2)*p2(2) + p2(2)**2)*c/12D0
         geom%l12(ib) = geom%l12(ib) + (p1(1)*p2(2) + 2D0*p1(1)*p1(2) + 2D0*p2(1)*p2(2) + p2(1)*p1(2))*c/24D0
      END DO
      DO ib = 1, mesh%nbody
         IF (geom%awp(ib) /= 0D0) THEN
            geom%xf(ib) = sx(ib)/geom%awp(ib)
            geom%yf(ib) = sy(ib)/geom%awp(ib)
         END IF
      END DO
   END SUBROUTINE waterplane_properties

   SUBROUTINE volume_properties(mesh, geom)
      ! Divergence theorem over the wetted panels, with the outward normal -nrm. The open waterplane at z = 0 adds
      ! nothing to int z n_z dS or to int x^2 n_x, int y^2 n_y, int z^2 n_z, so coordinates are absolute until the end.
      !   V = int z n_z dS,   int x dV = int x^2/2 n_x dS,   and the same for y and z.
      ! Panels are fanned into triangles; linear and squared linear integrands are integrated exactly on each.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(INOUT) :: geom
      REAL(8) :: mx(mesh%nbody), my(mesh%nbody), mz(mesh%nbody), v1(3), v2(3), v3(3), cr(3), at, no(3)
      INTEGER :: i, k, ib

      mx = 0D0
      my = 0D0
      mz = 0D0
      DO i = 1, mesh%n
         ib = mesh%body(i)
         no = -mesh%nrm(:, i)
         v1 = mesh%xnode(:, mesh%ncon(1, i))
         DO k = 2, mesh%ncn(i) - 1
            v2 = mesh%xnode(:, mesh%ncon(k, i))
            v3 = mesh%xnode(:, mesh%ncon(k + 1, i))
            cr = [(v2(2) - v1(2))*(v3(3) - v1(3)) - (v2(3) - v1(3))*(v3(2) - v1(2)), &
                  (v2(3) - v1(3))*(v3(1) - v1(1)) - (v2(1) - v1(1))*(v3(3) - v1(3)), &
                  (v2(1) - v1(1))*(v3(2) - v1(2)) - (v2(2) - v1(2))*(v3(1) - v1(1))]
            at = 0.5D0*SQRT(SUM(cr**2))
            geom%vol(ib) = geom%vol(ib) + at*(v1(3) + v2(3) + v3(3))/3D0*no(3)
            mx(ib) = mx(ib) + no(1)*at/12D0*quad_sum(v1(1), v2(1), v3(1))
            my(ib) = my(ib) + no(2)*at/12D0*quad_sum(v1(2), v2(2), v3(2))
            mz(ib) = mz(ib) + no(3)*at/12D0*quad_sum(v1(3), v2(3), v3(3))
         END DO
      END DO
      DO ib = 1, mesh%nbody
         IF (geom%vol(ib) > 0D0) THEN
            geom%xb(ib) = mx(ib)/geom%vol(ib) - mesh%xr(1, ib)
            geom%yb(ib) = my(ib)/geom%vol(ib) - mesh%xr(2, ib)
            geom%zb(ib) = mz(ib)/geom%vol(ib) - mesh%xr(3, ib)
         END IF
      END DO
   END SUBROUTINE volume_properties

   PURE REAL(8) FUNCTION quad_sum(a, b, c)
      ! 6 x the mean of f^2 over a triangle on which f is linear with vertex values a, b, c.
      REAL(8), INTENT(IN) :: a, b, c

      quad_sum = a*a + b*b + c*c + a*b + b*c + a*c
   END FUNCTION quad_sum

   SUBROUTINE build_neighbours(mesh, rep, nonconforming, sharp_deg, nbr_start, nbr)
      ! Neighbours of a panel: panels sharing a merged node. A panel with an open edge away from the waterline
      ! (nonconforming) also has the panels of the same body whose centroids are within 1.5 x the larger panel size, to
      ! reach across the gap. Other panels must not: this distance changes with the panel aspect ratio, and the velocity
      ! fit and the drift jump with it. Pairs with normals more than sharp_deg apart are dropped.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: rep(:)
      LOGICAL, INTENT(IN) :: nonconforming(:)
      REAL(8), INTENT(IN) :: sharp_deg
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nbr_start(:), nbr(:)
      INTEGER, ALLOCATABLE :: nstart(:), nlist(:), fillpos(:), mark(:), ord(:), rank(:)
      REAL(8), ALLOCATABLE :: key(:)
      REAL(8) :: cosmin, reach
      INTEGER :: n, pass, i, j, m, r, k, p, dirn

      n = mesh%n
      cosmin = COS(sharp_deg*PI/180D0) - 1D-12

      ! panels at each merged node
      ALLOCATE (nstart(mesh%nnode + 1), fillpos(mesh%nnode))
      nstart = 0
      DO i = 1, n
         DO m = 1, mesh%ncn(i)
            k = rep(mesh%ncon(m, i))
            nstart(k + 1) = nstart(k + 1) + 1
         END DO
      END DO
      nstart(1) = 1
      DO k = 1, mesh%nnode
         nstart(k + 1) = nstart(k) + nstart(k + 1)
      END DO
      fillpos = nstart(1:mesh%nnode)
      ALLOCATE (nlist(nstart(mesh%nnode + 1) - 1))
      DO i = 1, n
         DO m = 1, mesh%ncn(i)
            k = rep(mesh%ncon(m, i))
            nlist(fillpos(k)) = i
            fillpos(k) = fillpos(k) + 1
         END DO
      END DO

      ! panels in order of centroid x
      ALLOCATE (key(n), rank(n))
      key = mesh%xc(1, :)
      CALL sort_index(key, ord)
      DO r = 1, n
         rank(ord(r)) = r
      END DO
      reach = 1.5D0*MAXVAL(mesh%psize)

      ALLOCATE (nbr_start(n + 1), mark(n))
      ! pass 1 counts the neighbours, pass 2 stores them
      DO pass = 1, 2
         mark = 0
         nbr_start(1) = 1
         IF (pass == 2) ALLOCATE (nbr(nbr_start(n + 1) - 1))
         p = 0
         DO i = 1, n
            DO m = 1, mesh%ncn(i)   ! shared nodes
               k = rep(mesh%ncon(m, i))
               DO r = nstart(k), nstart(k + 1) - 1
                  CALL consider(nlist(r))
               END DO
            END DO
            DO dirn = -1, 1, 2   ! centroids close by
               IF (.NOT. nonconforming(i)) EXIT
               r = rank(i) + dirn
               DO WHILE (r >= 1 .AND. r <= n)
                  j = ord(r)
                  IF (ABS(mesh%xc(1, j) - mesh%xc(1, i)) > reach) EXIT
                  IF (mesh%body(j) == mesh%body(i) .AND. &
                      SQRT(SUM((mesh%xc(:, j) - mesh%xc(:, i))**2)) <= 1.5D0*MAX(mesh%psize(i), mesh%psize(j))) THEN
                     CALL consider(j)
                  END IF
                  r = r + dirn
               END DO
            END DO
            nbr_start(i + 1) = p + 1
         END DO
      END DO

   CONTAINS

      SUBROUTINE consider(j)
         INTEGER, INTENT(IN) :: j

         IF (j == i .OR. mark(j) == i) RETURN
         mark(j) = i
         IF (DOT_PRODUCT(mesh%nrm(:, i), mesh%nrm(:, j)) < cosmin) RETURN
         p = p + 1
         IF (pass == 2) nbr(p) = j
      END SUBROUTINE consider

   END SUBROUTINE build_neighbours

END MODULE SurfaceGeometry
