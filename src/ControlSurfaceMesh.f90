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
MODULE ControlSurfaceMesh
   ! Closed control surface C around each body for the control-surface mean drift, and the line Gamma_c where C meets
   ! the free surface z = 0. C is a vertical cylinder or box with a bottom; with a finite depth h it reaches the
   ! seabed and has no bottom. Its normals point towards the body, into the fluid domain between the body and C.
   ! Panel normals follow the HAMS rule (NormalProcess): for a quadrilateral n ~ (x4 - x2) x (x3 - x1), for a triangle
   ! n ~ (x1 - x2) x (x3 - x2). The wall rows are graded towards z = 0, where the waves and the waterline term need the
   ! finest sampling.
   ! The free-surface patch F is the part of z = 0 between the hull waterline and Gamma_c (Chen 2007, J. Eng. Math.
   ! 59:61-82, section 3.1). Its normal (0, 0, -1) points down, into the fluid domain.

   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, so_options, cs_points, CS_WALL, CS_LINE, CS_FREE
   USE SurfaceGeometry, ONLY: build_surface_geom, DEFAULT_SHARP_DEG
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cylinder_surface, box_surface, add_surface_points, build_control_surface, check_control_surface, &
             write_control_surface, read_control_mesh, free_surface_strips

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   REAL(8), PARAMETER :: WL_TOL = 1D-6        ! waterline tolerance for Gamma_c, relative to the size of C
   REAL(8), PARAMETER :: ROW_GROWTH = 1.25D0  ! thickness ratio of successive wall rows, until they reach the panel size
   REAL(8), PARAMETER :: DEFAULT_MARGIN_FACTOR = 4D0   ! default margin = 4 x the median hull panel size of the body
   REAL(8), PARAMETER :: FS_GRADING = 4D0     ! outer F row / inner F row thickness

CONTAINS

   ! ------------------------------------------------------------------------------------------------ surfaces

   SUBROUTINE cylinder_surface(centre, radius, depth, panel_size, h, cm)
      ! Cylinder of the given radius about the vertical axis through centre. depth is the wall depth; with h > 0 and
      ! depth >= h the wall reaches the seabed z = -h and there is no bottom.
      REAL(8), INTENT(IN) :: centre(2), radius, depth, panel_size, h
      TYPE(panel_mesh), INTENT(OUT) :: cm
      REAL(8), ALLOCATABLE :: zl(:)
      INTEGER :: np, nrow, nr, nw, nnode, npan, it, itn, iz, ir, ip
      REAL(8) :: th, r, dd
      LOGICAL :: bottom

      CALL wall_levels(depth, panel_size, h, zl, dd, bottom)
      nrow = SIZE(zl) - 1
      np = MAX(8, CEILING(2D0*PI*radius/panel_size))
      nr = MAX(1, CEILING(radius/panel_size))
      nw = np*(nrow + 1)
      nnode = nw
      npan = np*nrow
      IF (bottom) THEN
         nnode = nw + np*(nr - 1) + 1
         npan = npan + np*nr
      END IF
      CALL allocate_mesh(cm, npan, nnode, [centre(1), centre(2), 0D0])

      DO iz = 0, nrow
         DO it = 1, np
            th = 2D0*PI*REAL(it - 1, 8)/REAL(np, 8)
            cm%xnode(:, it + np*iz) = [centre(1) + radius*COS(th), centre(2) + radius*SIN(th), zl(iz)]
         END DO
      END DO
      IF (bottom) THEN
         DO ir = 1, nr - 1
            r = radius*REAL(ir, 8)/REAL(nr, 8)
            DO it = 1, np
               th = 2D0*PI*REAL(it - 1, 8)/REAL(np, 8)
               cm%xnode(:, nw + (ir - 1)*np + it) = [centre(1) + r*COS(th), centre(2) + r*SIN(th), -dd]
            END DO
         END DO
         cm%xnode(:, nnode) = [centre(1), centre(2), -dd]
      END IF

      ip = 0
      DO iz = 1, nrow
         DO it = 1, np
            itn = MOD(it, np) + 1
            CALL add_panel(cm, ip, [it + np*(iz - 1), it + np*iz, itn + np*iz, itn + np*(iz - 1)])
         END DO
      END DO
      IF (bottom) THEN
         DO ir = 1, nr
            DO it = 1, np
               itn = MOD(it, np) + 1
               IF (ir == 1) THEN
                  CALL add_panel(cm, ip, [nnode, ring_node(1, itn), ring_node(1, it)])
               ELSE
                  CALL add_panel(cm, ip, [ring_node(ir - 1, it), ring_node(ir - 1, itn), ring_node(ir, itn), ring_node(ir, it)])
               END IF
            END DO
         END DO
      END IF

   CONTAINS

      INTEGER FUNCTION ring_node(ring, j)
         INTEGER, INTENT(IN) :: ring, j
         IF (ring == nr) THEN
            ring_node = j + np*nrow
         ELSE
            ring_node = nw + (ring - 1)*np + j
         END IF
      END FUNCTION ring_node

   END SUBROUTINE cylinder_surface

   SUBROUTINE box_surface(centre, half, angle, depth, panel_size, h, cm)
      ! Rectangular box with half-lengths half(1), half(2) along its own axes, which are the global axes rotated by
      ! angle about z. Depth and seabed as in cylinder_surface.
      REAL(8), INTENT(IN) :: centre(2), half(2), angle, depth, panel_size, h
      TYPE(panel_mesh), INTENT(OUT) :: cm
      REAL(8), ALLOCATABLE :: zl(:), ls(:), lt(:)
      INTEGER, ALLOCATABLE :: gid(:,:)
      INTEGER :: nx, ny, np, nrow, nw, nnode, npan, it, itn, iz, i, j, ip, nint
      REAL(8) :: dd, cs, sn
      LOGICAL :: bottom

      CALL wall_levels(depth, panel_size, h, zl, dd, bottom)
      nrow = SIZE(zl) - 1
      nx = MAX(1, CEILING(2D0*half(1)/panel_size))
      ny = MAX(1, CEILING(2D0*half(2)/panel_size))
      np = 2*(nx + ny)
      nw = np*(nrow + 1)
      nnode = nw
      npan = np*nrow
      IF (bottom) THEN
         nnode = nw + (nx - 1)*(ny - 1)
         npan = npan + nx*ny
      END IF
      CALL allocate_mesh(cm, npan, nnode, [centre(1), centre(2), 0D0])
      cs = COS(angle)
      sn = SIN(angle)

      ! Loop of wall points, counter-clockwise from the corner (-hx, -hy), in the box's own axes.
      ALLOCATE (ls(np), lt(np))
      DO j = 0, nx - 1
         ls(j + 1) = -half(1) + 2D0*half(1)*REAL(j, 8)/REAL(nx, 8);   lt(j + 1) = -half(2)
         ls(nx + ny + j + 1) = half(1) - 2D0*half(1)*REAL(j, 8)/REAL(nx, 8);   lt(nx + ny + j + 1) = half(2)
      END DO
      DO j = 0, ny - 1
         ls(nx + j + 1) = half(1);   lt(nx + j + 1) = -half(2) + 2D0*half(2)*REAL(j, 8)/REAL(ny, 8)
         ls(2*nx + ny + j + 1) = -half(1);   lt(2*nx + ny + j + 1) = half(2) - 2D0*half(2)*REAL(j, 8)/REAL(ny, 8)
      END DO
      DO iz = 0, nrow
         DO it = 1, np
            cm%xnode(:, it + np*iz) = [centre(1) + cs*ls(it) - sn*lt(it), centre(2) + sn*ls(it) + cs*lt(it), zl(iz)]
         END DO
      END DO

      ip = 0
      DO iz = 1, nrow
         DO it = 1, np
            itn = MOD(it, np) + 1
            CALL add_panel(cm, ip, [it + np*(iz - 1), it + np*iz, itn + np*iz, itn + np*(iz - 1)])
         END DO
      END DO

      IF (bottom) THEN
         ! Grid node (i, j), i = 0..nx, j = 0..ny: the boundary ones are the bottom wall points, the rest are new.
         ALLOCATE (gid(0:nx, 0:ny))
         nint = nw
         DO j = 0, ny
            DO i = 0, nx
               IF (j == 0 .AND. i < nx) THEN
                  gid(i, j) = i + 1
               ELSE IF (i == nx .AND. j < ny) THEN
                  gid(i, j) = nx + j + 1
               ELSE IF (j == ny .AND. i > 0) THEN
                  gid(i, j) = nx + ny + (nx - i) + 1
               ELSE IF (i == 0 .AND. j > 0) THEN
                  gid(i, j) = 2*nx + ny + (ny - j) + 1
               ELSE IF (i == nx .AND. j == 0) THEN
                  gid(i, j) = nx + 1
               ELSE IF (i == 0 .AND. j == 0) THEN
                  gid(i, j) = 1
               END IF
               IF (i > 0 .AND. i < nx .AND. j > 0 .AND. j < ny) THEN
                  nint = nint + 1
                  cm%xnode(:, nint) = [centre(1) + cs*(-half(1) + 2D0*half(1)*REAL(i, 8)/REAL(nx, 8)) &
                                                  - sn*(-half(2) + 2D0*half(2)*REAL(j, 8)/REAL(ny, 8)), &
                                       centre(2) + sn*(-half(1) + 2D0*half(1)*REAL(i, 8)/REAL(nx, 8)) &
                                                  + cs*(-half(2) + 2D0*half(2)*REAL(j, 8)/REAL(ny, 8)), -dd]
                  gid(i, j) = nint
               ELSE
                  gid(i, j) = gid(i, j) + np*nrow   ! wall point at the bottom level
               END IF
            END DO
         END DO
         DO j = 0, ny - 1
            DO i = 0, nx - 1
               CALL add_panel(cm, ip, [gid(i, j), gid(i, j + 1), gid(i + 1, j + 1), gid(i + 1, j)])
            END DO
         END DO
      END IF
   END SUBROUTINE box_surface

   SUBROUTINE wall_levels(depth, panel_size, h, zl, dd, bottom)
      ! Depths zl(0:nrow) of the wall rows: zl(0) = 0, zl(nrow) = -dd. The top row is panel_size/2 thick and the rows
      ! grow by ROW_GROWTH to panel_size, then stay at panel_size. Rows below the first are scaled so that the last
      ! level is exactly -dd.
      REAL(8), INTENT(IN) :: depth, panel_size, h
      REAL(8), ALLOCATABLE, INTENT(OUT) :: zl(:)
      REAL(8), INTENT(OUT) :: dd
      LOGICAL, INTENT(OUT) :: bottom
      REAL(8), ALLOCATABLE :: rowt(:)
      REAL(8) :: cum, tk, t1
      INTEGER :: nrow, k

      bottom = .TRUE.
      dd = depth
      IF (h > 0D0 .AND. depth >= h) THEN
         dd = h
         bottom = .FALSE.
      END IF
      t1 = 0.5D0*panel_size
      ! Number of rows: add rows while at least half of the next one still fits.
      nrow = 1
      cum = t1
      tk = MIN(panel_size, t1*ROW_GROWTH)
      DO WHILE (cum + 0.5D0*tk < dd)
         nrow = nrow + 1
         cum = cum + tk
         tk = MIN(panel_size, tk*ROW_GROWTH)
      END DO
      ALLOCATE (rowt(nrow), zl(0:nrow))
      tk = t1
      DO k = 1, nrow
         rowt(k) = tk
         tk = MIN(panel_size, tk*ROW_GROWTH)
      END DO
      IF (nrow == 1) THEN
         rowt(1) = dd
      ELSE
         rowt(2:nrow) = rowt(2:nrow)*(dd - t1)/(SUM(rowt) - t1)
      END IF
      zl(0) = 0D0
      DO k = 1, nrow
         zl(k) = zl(k - 1) - rowt(k)
      END DO
      zl(nrow) = -dd
   END SUBROUTINE wall_levels

   SUBROUTINE allocate_mesh(cm, npan, nnode, xr)
      TYPE(panel_mesh), INTENT(INOUT) :: cm
      INTEGER, INTENT(IN) :: npan, nnode
      REAL(8), INTENT(IN) :: xr(3)
      INTEGER :: i

      cm%n = npan
      cm%nbase = npan
      cm%nnode = nnode
      cm%nbody = 1
      cm%nsys = 1
      ALLOCATE (cm%xc(3, npan), cm%nrm(3, npan), cm%area(npan), cm%psize(npan), cm%nj(6, npan))
      ALLOCATE (cm%xnode(3, nnode), cm%ncn(npan), cm%ncon(4, npan), cm%body(npan), cm%base(npan), cm%image(npan))
      ALLOCATE (cm%xr(3, 1))
      cm%xr(:, 1) = xr
      cm%ncon = 0
      cm%body = 1
      cm%image = 1
      cm%base = [(i, i = 1, npan)]
   END SUBROUTINE allocate_mesh

   SUBROUTINE add_panel(cm, ip, nodes)
      ! Panel ip+1 from 3 or 4 node numbers: centroid, area, size and the HAMS normal (see the module header).
      TYPE(panel_mesh), INTENT(INOUT) :: cm
      INTEGER, INTENT(INOUT) :: ip
      INTEGER, INTENT(IN) :: nodes(:)
      REAL(8) :: v(3, 4), d1(3), d2(3), c(3), un(3), x(3)
      INTEGER :: m, nv

      nv = SIZE(nodes)
      ip = ip + 1
      cm%ncn(ip) = nv
      cm%ncon(1:nv, ip) = nodes
      DO m = 1, nv
         v(:, m) = cm%xnode(:, nodes(m))
      END DO
      cm%xc(:, ip) = SUM(v(:, 1:nv), DIM=2)/REAL(nv, 8)
      IF (nv == 4) THEN
         d1 = v(:, 3) - v(:, 1)
         d2 = v(:, 4) - v(:, 2)
         c = cross(d2, d1)
      ELSE
         d1 = v(:, 2) - v(:, 1)
         d2 = v(:, 3) - v(:, 1)
         c = cross(v(:, 1) - v(:, 2), v(:, 3) - v(:, 2))
      END IF
      cm%area(ip) = 0.5D0*SQRT(SUM(cross(d1, d2)**2))
      cm%psize(ip) = SQRT(cm%area(ip))
      un = c/SQRT(SUM(c**2))
      cm%nrm(:, ip) = un
      x = cm%xc(:, ip) - cm%xr(:, 1)
      cm%nj(1:3, ip) = un
      cm%nj(4, ip) = x(2)*un(3) - x(3)*un(2)
      cm%nj(5, ip) = x(3)*un(1) - x(1)*un(3)
      cm%nj(6, ip) = x(1)*un(2) - x(2)*un(1)
   END SUBROUTINE add_panel

   PURE FUNCTION cross(a, b) RESULT(c)
      REAL(8), INTENT(IN) :: a(3), b(3)
      REAL(8) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross

   ! ------------------------------------------------------------------------------------------------ points

   SUBROUTINE add_surface_points(cm, b, pts, ok, msg)
      ! Appends the panels of C as CS_WALL points of body b and the segments of Gamma_c as CS_LINE points. Gamma_c comes
      ! from build_surface_geom: the segment midpoint wl_xm, the length wl_len and the horizontal unit normal wl_nh.
      ! The panel normal of C points towards the body, and so does wl_nh (the horizontal part of that normal), which
      ! is the direction n-bar into F that the control-surface method needs. No sign change is needed.
      ! ok and msg are optional: ok = .false. only for a stop of build_surface_geom, in which case no point is added.
      TYPE(panel_mesh), INTENT(IN) :: cm
      INTEGER, INTENT(IN) :: b
      TYPE(cs_points), INTENT(INOUT) :: pts
      LOGICAL, INTENT(OUT), OPTIONAL :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: msg
      TYPE(surface_geom) :: g
      LOGICAL :: ok0
      CHARACTER(LEN=:), ALLOCATABLE :: msg0
      INTEGER :: nw, nl, n0

      CALL build_surface_geom(cm, WL_TOL, DEFAULT_SHARP_DEG, g, ok0, msg0)
      IF (PRESENT(ok)) ok = ok0
      IF (PRESENT(msg)) msg = msg0
      IF (.NOT. ok0) RETURN
      nw = cm%n
      nl = g%nwl
      n0 = pts%n
      CALL grow_points(pts, nw + nl)
      pts%x(:, n0 + 1:n0 + nw) = cm%xc
      pts%nrm(:, n0 + 1:n0 + nw) = cm%nrm
      pts%w(n0 + 1:n0 + nw) = cm%area
      pts%psize(n0 + 1:n0 + nw) = cm%psize
      pts%kind(n0 + 1:n0 + nw) = CS_WALL
      pts%body(n0 + 1:n0 + nw) = b
      IF (nl > 0) THEN
         pts%x(:, n0 + nw + 1:n0 + nw + nl) = g%wl_xm
         pts%nrm(:, n0 + nw + 1:n0 + nw + nl) = g%wl_nh
         pts%w(n0 + nw + 1:n0 + nw + nl) = g%wl_len
         pts%psize(n0 + nw + 1:n0 + nw + nl) = g%wl_len
         pts%kind(n0 + nw + 1:n0 + nw + nl) = CS_LINE
         pts%body(n0 + nw + 1:n0 + nw + nl) = b
      END IF
   END SUBROUTINE add_surface_points

   SUBROUTINE grow_points(pts, nadd)
      ! Extends pts by nadd undefined entries, keeping what is there.
      TYPE(cs_points), INTENT(INOUT) :: pts
      INTEGER, INTENT(IN) :: nadd
      TYPE(cs_points) :: t
      INTEGER :: n0, n1

      n0 = pts%n
      n1 = n0 + nadd
      ALLOCATE (t%x(3, n1), t%nrm(3, n1), t%w(n1), t%psize(n1), t%kind(n1), t%body(n1))
      IF (n0 > 0) THEN
         t%x(:, 1:n0) = pts%x
         t%nrm(:, 1:n0) = pts%nrm
         t%w(1:n0) = pts%w
         t%psize(1:n0) = pts%psize
         t%kind(1:n0) = pts%kind
         t%body(1:n0) = pts%body
      END IF
      CALL MOVE_ALLOC(t%x, pts%x)
      CALL MOVE_ALLOC(t%nrm, pts%nrm)
      CALL MOVE_ALLOC(t%w, pts%w)
      CALL MOVE_ALLOC(t%psize, pts%psize)
      CALL MOVE_ALLOC(t%kind, pts%kind)
      CALL MOVE_ALLOC(t%body, pts%body)
      pts%n = n1
   END SUBROUTINE grow_points

   ! ------------------------------------------------------------------------------------------------ build

   SUBROUTINE build_control_surface(mesh, geom, opts, h, inputdir, pts, cms, ok, msg, with_free, free_built)
      ! Generates C for every body and turns it into evaluation points (default margins 4 x the median hull panel size).
      ! h < 0 is deep water.
      ! cs_shape = 3 reads C of body b from inputdir/CSMesh_<b>.pnl and checks that its normals point towards the body.
      ! For a file C, ok stays .true. when build_surface_geom only warns, and the warnings are returned in msg.
      ! with_free (default .false.) also appends F of each body, after the C points of that body: see add_free_surface.
      ! F is built also when C reaches the seabed. When a generated F of body b cannot be built (free_surface_strips), ok
      ! stays .true., msg says why and free_built(b) is .false.; the other bodies keep their F. An FSMesh_<b>.pnl that
      ! cannot be read or is not in z = 0 sets ok = .false.: the user asked for that F.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(so_options), INTENT(IN) :: opts
      REAL(8), INTENT(IN) :: h
      CHARACTER(LEN=*), INTENT(IN) :: inputdir
      TYPE(cs_points), INTENT(OUT) :: pts
      TYPE(panel_mesh), ALLOCATABLE, INTENT(OUT) :: cms(:)
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      LOGICAL, INTENT(IN), OPTIONAL :: with_free
      LOGICAL, ALLOCATABLE, INTENT(OUT), OPTIONAL :: free_built(:)   ! (nbody), F of body b was appended
      REAL(8), ALLOCATABLE :: psz(:)
      REAL(8) :: lo(2), hi(2), c(2), zmin, med, margin, bmargin, ps, radius, half(2), depth, xy(3)
      INTEGER :: b, i, m, nfound
      LOGICAL :: ok1, wf, from_file
      CHARACTER(LEN=:), ALLOCATABLE :: msg1
      CHARACTER(LEN=200) :: line
      CHARACTER(LEN=16) :: istr

      ok = .TRUE.
      msg = ''
      wf = .FALSE.
      IF (PRESENT(with_free)) wf = with_free
      IF (PRESENT(free_built)) THEN
         ALLOCATE (free_built(mesh%nbody))
         free_built = .FALSE.
      END IF
      IF (opts%cs_shape /= 1 .AND. opts%cs_shape /= 2 .AND. opts%cs_shape /= 3) THEN
         ok = .FALSE.
         WRITE (line, '(a,i0,a)') 'cs_shape = ', opts%cs_shape, ' is not known; use cylinder, box or file'
         msg = TRIM(line)
         RETURN
      END IF
      ALLOCATE (cms(mesh%nbody))
      DO b = 1, mesh%nbody
         CALL waterline_box(geom, b, lo, hi, nfound)
         IF (nfound == 0) THEN
            ok = .FALSE.
            WRITE (line, '(a,i0,a)') 'body ', b, ' has no waterline; a control surface needs one'
            msg = TRIM(line)
            RETURN
         END IF
         c = 0.5D0*(lo + hi)
         ! Median panel size of the body; the panel size of C and F.
         psz = PACK(mesh%psize, mesh%body == b)
         med = median(psz)
         ps = opts%cs_panel_size
         IF (.NOT. ps > 0D0) ps = med
         IF (opts%cs_shape == 3) THEN
            ! C from the file CSMesh_<b>.pnl. The warnings of build_surface_geom (gaps in the waterline of C) are
            ! appended to msg, with ok still .true.; the user's file can have such gaps.
            WRITE (istr, '(i0)') b
            CALL read_file_surface(TRIM(inputdir)//'/CSMesh_'//TRIM(istr)//'.pnl', 'CSMesh_'//TRIM(istr)//'.pnl', c, cms(b), &
                                   ok1, msg1)
            IF (.NOT. ok1) THEN
               ok = .FALSE.
               msg = msg1
               RETURN
            END IF
            CALL add_surface_points(cms(b), b, pts, ok1, msg1)
            IF (.NOT. ok1) THEN
               ok = .FALSE.
               WRITE (line, '(a,i0,a)') 'control surface of body ', b, ': '
               msg = TRIM(line)//msg1
               RETURN
            END IF
            IF (LEN(msg1) > 0) THEN
               WRITE (line, '(a,i0,a)') 'Warning, control surface of body ', b, ': '
               IF (LEN(msg) > 0) msg = msg//NEW_LINE('a')
               msg = msg//TRIM(line)//msg1
            END IF
            IF (wf) THEN
               CALL add_free_surface(mesh, geom, b, opts%cs_fs_rows, ps, inputdir, pts, ok1, msg1, from_file)
               CALL free_result(b, ok1, msg1, from_file)
               IF (.NOT. ok) RETURN
            END IF
            CYCLE
         END IF
         ! Deepest node and largest horizontal extent from c.
         zmin = 0D0
         radius = 0D0
         half = 0D0
         DO i = 1, mesh%n
            IF (mesh%body(i) /= b) CYCLE
            DO m = 1, mesh%ncn(i)
               xy = mesh%xnode(:, mesh%ncon(m, i))
               zmin = MIN(zmin, xy(3))
               radius = MAX(radius, HYPOT(xy(1) - c(1), xy(2) - c(2)))
               half = MAX(half, ABS(xy(1:2) - c))
            END DO
         END DO
         margin = opts%cs_margin
         IF (.NOT. margin > 0D0) margin = DEFAULT_MARGIN_FACTOR*med
         bmargin = opts%cs_bottom_margin
         IF (.NOT. bmargin > 0D0) bmargin = DEFAULT_MARGIN_FACTOR*med
         depth = -zmin + bmargin
         IF (opts%cs_shape == 1) THEN
            CALL cylinder_surface(c, radius + margin, depth, ps, h, cms(b))
         ELSE
            ! The box follows the global axes: the mesh does not carry the LCS angle of the body.
            CALL box_surface(c, half + margin, 0D0, depth, ps, h, cms(b))
         END IF
         CALL add_surface_points(cms(b), b, pts, ok1, msg1)
         IF (.NOT. ok1) THEN
            ok = .FALSE.
            WRITE (line, '(a,i0,a)') 'control surface of body ', b, ': '
            msg = TRIM(line)//msg1
            RETURN
         END IF
         IF (wf) THEN
            CALL add_free_surface(mesh, geom, b, opts%cs_fs_rows, ps, inputdir, pts, ok1, msg1, from_file)
            CALL free_result(b, ok1, msg1, from_file)
            IF (.NOT. ok) RETURN
         END IF
      END DO

   CONTAINS

      SUBROUTINE free_result(b, ok1, msg1, from_file)
         ! The outcome of add_free_surface for body b: a failing FSMesh file stops the build, a failing generated F is a
         ! warning in msg.
         INTEGER, INTENT(IN) :: b
         LOGICAL, INTENT(IN) :: ok1, from_file
         CHARACTER(LEN=*), INTENT(IN) :: msg1
         CHARACTER(LEN=16) :: bstr

         IF (ok1) THEN
            IF (PRESENT(free_built)) free_built(b) = .TRUE.
         ELSE IF (from_file) THEN
            ok = .FALSE.
            msg = msg1
         ELSE
            WRITE (bstr, '(i0)') b
            IF (LEN(msg) > 0) msg = msg//NEW_LINE('a')
            msg = msg//'Warning, '//msg1//' No free surface F, so F_z, M_x, M_y of body '//TRIM(bstr)//' are not computed.'
         END IF
      END SUBROUTINE free_result

   END SUBROUTINE build_control_surface

   SUBROUTINE waterline_box(geom, b, lo, hi, nfound)
      ! Horizontal bounding box lo, hi of the nfound hull waterline segments of body b.
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: b
      REAL(8), INTENT(OUT) :: lo(2), hi(2)
      INTEGER, INTENT(OUT) :: nfound
      INTEGER :: iw

      lo = HUGE(1D0)
      hi = -HUGE(1D0)
      nfound = 0
      DO iw = 1, geom%nwl
         IF (geom%wl_body(iw) /= b) CYCLE
         nfound = nfound + 1
         lo = MIN(lo, MIN(geom%wl_x1(1:2, iw), geom%wl_x2(1:2, iw)))
         hi = MAX(hi, MAX(geom%wl_x1(1:2, iw), geom%wl_x2(1:2, iw)))
      END DO
   END SUBROUTINE waterline_box

   ! ------------------------------------------------------------------------------------------------ free surface

   SUBROUTINE add_free_surface(mesh, geom, b, rows, ps, inputdir, pts, ok, msg, from_file)
      ! Appends F of body b to pts as CS_FREE points with nrm = (0, 0, -1), weight the panel area. F is read from
      ! inputdir/FSMesh_<b>.pnl when that file exists (its nodes must lie in z = 0; the node order does not matter),
      ! else it is built by free_surface_strips. The strips have rows rows; rows = 0 takes the default: the nearest
      ! integer to (d_c - d_w)/ps, at least 2, with d_c and d_w the mean distances from the centre of the Gamma_c
      ! segment midpoints and of the hull waterline segment midpoints, and ps the panel size of C.
      ! from_file is .true. when F came, or was to come, from FSMesh_<b>.pnl.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: b, rows
      REAL(8), INTENT(IN) :: ps
      CHARACTER(LEN=*), INTENT(IN) :: inputdir
      TYPE(cs_points), INTENT(INOUT) :: pts
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      LOGICAL, INTENT(OUT) :: from_file
      TYPE(panel_mesh) :: fm
      CHARACTER(LEN=16) :: istr
      CHARACTER(LEN=:), ALLOCATABLE :: name
      REAL(8) :: lo(2), hi(2), c(2), dc, dw
      INTEGER :: nrows, n0, iw, p, nw, nl
      LOGICAL :: exists

      WRITE (istr, '(i0)') b
      name = 'FSMesh_'//TRIM(istr)//'.pnl'
      INQUIRE (FILE=TRIM(inputdir)//'/'//name, EXIST=exists)
      from_file = exists
      IF (exists) THEN
         CALL read_control_mesh(TRIM(inputdir)//'/'//name, fm, ok, msg)
         IF (.NOT. ok) RETURN
         IF (MAXVAL(ABS(fm%xnode(3, :))) > WL_TOL*MAX(1D0, MAXVAL(ABS(fm%xnode(1:2, :))))) THEN
            ok = .FALSE.
            msg = name//': F must lie in the free surface z = 0.'
            RETURN
         END IF
      ELSE
         nrows = rows
         IF (nrows < 1) THEN
            CALL waterline_box(geom, b, lo, hi, nw)
            c = 0.5D0*(lo + hi)
            dw = 0D0
            DO iw = 1, geom%nwl
               IF (geom%wl_body(iw) == b) dw = dw + HYPOT(geom%wl_xm(1, iw) - c(1), geom%wl_xm(2, iw) - c(2))
            END DO
            dc = 0D0
            nl = 0
            DO p = 1, pts%n
               IF (pts%kind(p) /= CS_LINE .OR. pts%body(p) /= b) CYCLE
               nl = nl + 1
               dc = dc + HYPOT(pts%x(1, p) - c(1), pts%x(2, p) - c(2))
            END DO
            nrows = 2
            IF (nw > 0 .AND. nl > 0 .AND. ps > 0D0) nrows = MAX(2, NINT((dc/REAL(nl, 8) - dw/REAL(nw, 8))/ps))
         END IF
         CALL free_surface_strips(mesh, geom, b, pts, nrows, fm, ok, msg)
         IF (.NOT. ok) RETURN
      END IF
      n0 = pts%n
      CALL grow_points(pts, fm%n)
      pts%x(:, n0 + 1:) = fm%xc
      pts%nrm(1:2, n0 + 1:) = 0D0
      pts%nrm(3, n0 + 1:) = -1D0
      pts%w(n0 + 1:) = fm%area
      pts%psize(n0 + 1:) = fm%psize
      pts%kind(n0 + 1:) = CS_FREE
      pts%body(n0 + 1:) = b
   END SUBROUTINE add_free_surface

   SUBROUTINE free_surface_strips(mesh, geom, b, pts, nrows, fm, ok, msg)
      ! F of body b as quadrilateral strips in z = 0, from each hull waterline segment of body b out to Gamma_c of
      ! body b (its CS_LINE points in pts), along rays from the centre x_b of the waterplane bounding box (the centre
      ! of C). Each segment is split into columns about as wide as the Gamma_c segments where the strip meets Gamma_c,
      ! so that the outer edge of F follows Gamma_c. The outer nodes are the intersections of the rays with the Gamma_c
      ! segments. Along each ray there are nrows rows, graded towards the hull: the thicknesses grow geometrically and
      ! the inner row is 1/FS_GRADING of the outer one.
      ! The rays only cover F when the waterline is star-shaped about x_b: the polar angle about x_b must increase along
      ! every segment (they run counter-clockwise, see surface_geom) and turn once in total. Else ok = .false.
      ! Panel normals are (0, 0, -1), by the HAMS rule of add_panel. mesh is not used yet: geom and pts hold all the
      ! geometry the strips need.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: b, nrows
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(panel_mesh), INTENT(OUT) :: fm
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      REAL(8), ALLOCATABLE :: s1(:,:), s2(:,:), o1(:,:), o2(:,:), sk(:)
      INTEGER, ALLOCATABLE :: seg(:), ncol(:)
      REAL(8) :: lo(2), hi(2), c(2), v1(2), v2(2), turn, dth, lmean, q, xi(2), xo(2)
      INTEGER :: nr, nw, nl, iw, j, k, p, js, nnode, npan, ip, base
      CHARACTER(LEN=200) :: line
      LOGICAL :: hit1, hit2

      ok = .FALSE.
      msg = ''
      nr = MAX(1, nrows)
      CALL waterline_box(geom, b, lo, hi, nw)
      IF (nw == 0) THEN
         WRITE (line, '(a,i0,a)') 'body ', b, ' has no waterline; a control surface needs one'
         msg = TRIM(line)
         RETURN
      END IF
      c = 0.5D0*(lo + hi)
      seg = PACK([(iw, iw = 1, geom%nwl)], geom%wl_body(1:geom%nwl) == b)

      ! Star-shape test: every segment turns counter-clockwise about c, and all of them turn once.
      turn = 0D0
      DO js = 1, nw
         v1 = geom%wl_x1(1:2, seg(js)) - c
         v2 = geom%wl_x2(1:2, seg(js)) - c
         dth = ATAN2(v1(1)*v2(2) - v1(2)*v2(1), DOT_PRODUCT(v1, v2))
         IF (.NOT. dth > 0D0) EXIT
         turn = turn + dth
      END DO
      IF (js <= nw .OR. ABS(turn - 2D0*PI) > 1D-6*2D0*PI) THEN
         WRITE (line, '(a,i0,a,i0,a)') 'body ', b, ': the waterline is not star-shaped about its centre; give FSMesh_', &
            b, '.pnl.'
         msg = TRIM(line)
         RETURN
      END IF

      ! Gamma_c segments of body b, from the stored midpoint, normal and length (as in check_control_surface).
      nl = COUNT(pts%kind == CS_LINE .AND. pts%body == b)
      IF (nl == 0) THEN
         WRITE (line, '(a,i0,a)') 'body ', b, ': no Gamma_c to build the free surface F up to'
         msg = TRIM(line)
         RETURN
      END IF
      ALLOCATE (s1(2, nl), s2(2, nl))
      k = 0
      lmean = 0D0
      DO p = 1, pts%n
         IF (pts%kind(p) /= CS_LINE .OR. pts%body(p) /= b) CYCLE
         k = k + 1
         s1(:, k) = pts%x(1:2, p) - 0.5D0*pts%w(p)*[-pts%nrm(2, p), pts%nrm(1, p)]
         s2(:, k) = pts%x(1:2, p) + 0.5D0*pts%w(p)*[-pts%nrm(2, p), pts%nrm(1, p)]
         lmean = lmean + pts%w(p)
      END DO
      lmean = lmean/REAL(nl, 8)

      ! Outer end points of each segment's rays, and the number of columns.
      ALLOCATE (o1(2, nw), o2(2, nw), ncol(nw))
      DO js = 1, nw
         CALL ray_to_line(c, geom%wl_x1(1:2, seg(js)), s1, s2, o1(:, js), hit1)
         CALL ray_to_line(c, geom%wl_x2(1:2, seg(js)), s1, s2, o2(:, js), hit2)
         IF (.NOT. (hit1 .AND. hit2)) THEN
            WRITE (line, '(a,i0,a,i0,a)') 'body ', b, ': a ray from the centre through the waterline does not meet ' &
               //'Gamma_c outside the waterline; give FSMesh_', b, '.pnl.'
            msg = TRIM(line)
            RETURN
         END IF
         ncol(js) = MAX(1, CEILING(NORM2(o2(:, js) - o1(:, js))/lmean))
      END DO

      ! Radial node fractions sk(0:nr), 0 at the hull and 1 on Gamma_c; thicknesses in the ratio q per row.
      ALLOCATE (sk(0:nr))
      IF (nr == 1) THEN
         sk = [0D0, 1D0]
      ELSE
         q = FS_GRADING**(1D0/REAL(nr - 1, 8))
         sk = [((q**k - 1D0)/(q**nr - 1D0), k = 0, nr)]
         sk(nr) = 1D0
      END IF

      npan = SUM(ncol)*nr
      nnode = SUM(ncol + 1)*(nr + 1)
      CALL allocate_mesh(fm, npan, nnode, [c(1), c(2), 0D0])
      ip = 0
      base = 0
      DO js = 1, nw
         ! Node (j, k) of this strip: column j = 0..ncol along the segment, row level k = 0..nr out from the hull.
         DO j = 0, ncol(js)
            IF (j == 0) THEN
               xi = geom%wl_x1(1:2, seg(js))
               xo = o1(:, js)
            ELSE IF (j == ncol(js)) THEN
               xi = geom%wl_x2(1:2, seg(js))
               xo = o2(:, js)
            ELSE
               xi = geom%wl_x1(1:2, seg(js)) + REAL(j, 8)/REAL(ncol(js), 8)*(geom%wl_x2(1:2, seg(js)) - &
                                                                              geom%wl_x1(1:2, seg(js)))
               CALL ray_to_line(c, xi, s1, s2, xo, hit1)
               IF (.NOT. hit1) THEN
                  WRITE (line, '(a,i0,a,i0,a)') 'body ', b, ': a ray from the centre through the waterline does not ' &
                     //'meet Gamma_c outside the waterline; give FSMesh_', b, '.pnl.'
                  msg = TRIM(line)
                  RETURN
               END IF
            END IF
            DO k = 0, nr
               fm%xnode(:, node(j, k)) = [xi + sk(k)*(xo - xi), 0D0]
            END DO
         END DO
         ! The order inner j, outer j, outer j+1, inner j+1 gives the normal (0, 0, -1) with the HAMS rule.
         DO j = 0, ncol(js) - 1
            DO k = 0, nr - 1
               CALL add_panel(fm, ip, [node(j, k), node(j, k + 1), node(j + 1, k + 1), node(j + 1, k)])
            END DO
         END DO
         base = base + (ncol(js) + 1)*(nr + 1)
      END DO
      ok = .TRUE.

   CONTAINS

      INTEGER FUNCTION node(j, k)
         INTEGER, INTENT(IN) :: j, k
         node = base + j*(nr + 1) + k + 1
      END FUNCTION node

   END SUBROUTINE free_surface_strips

   SUBROUTINE ray_to_line(c, p, s1, s2, x, hit)
      ! First crossing x of the ray from c through p with the segments s1(:,k) to s2(:,k), beyond p. hit = .false.
      ! when the ray meets no segment beyond p.
      REAL(8), INTENT(IN) :: c(2), p(2), s1(:,:), s2(:,:)
      REAL(8), INTENT(OUT) :: x(2)
      LOGICAL, INTENT(OUT) :: hit
      REAL(8), PARAMETER :: UTOL = 1D-9
      REAL(8) :: d(2), e(2), w(2), den, t, u, tbest
      INTEGER :: k

      d = p - c
      tbest = HUGE(1D0)
      DO k = 1, SIZE(s1, 2)
         e = s2(:, k) - s1(:, k)
         den = d(1)*e(2) - d(2)*e(1)
         IF (ABS(den) <= 1D-14*NORM2(d)*NORM2(e)) CYCLE
         w = s1(:, k) - c
         t = (w(1)*e(2) - w(2)*e(1))/den
         u = (w(1)*d(2) - w(2)*d(1))/den
         IF (u < -UTOL .OR. u > 1D0 + UTOL .OR. .NOT. t > 1D0) CYCLE
         tbest = MIN(tbest, t)
      END DO
      hit = tbest < HUGE(1D0)
      x = p
      IF (hit) x = c + tbest*d
   END SUBROUTINE ray_to_line

   ! ------------------------------------------------------------------------------------------------ file

   SUBROUTINE read_control_mesh(path, cm, ok, msg)
      ! Reads C from a file in the HullMesh.pnl format (same reading sequence as HAMS_Prog). C is always full, so the
      ! symmetry flags ISX and ISY must be 0. Centroids, areas, sizes and normals come from add_panel, as for a
      ! generated C. xr is set to the origin: build_control_surface moves it to the body centre.
      CHARACTER(LEN=*), INTENT(IN) :: path
      TYPE(panel_mesh), INTENT(OUT) :: cm
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      INTEGER :: u, ios, k, np, nn, isx, isy, idx, nv, v(4), ip
      REAL(8) :: x(3), d1(3), d2(3)
      CHARACTER(LEN=512) :: rec
      CHARACTER(LEN=300) :: line

      ok = .FALSE.
      msg = ''
      OPEN (NEWUNIT=u, FILE=path, STATUS='old', ACTION='read', IOSTAT=ios)
      IF (ios /= 0) THEN
         msg = 'cannot open '//path
         RETURN
      END IF
      DO k = 1, 3
         READ (u, *, IOSTAT=ios)
         IF (ios /= 0) GO TO 900
      END DO
      READ (u, *, IOSTAT=ios) np, nn, isx, isy
      IF (ios /= 0) GO TO 900
      IF (isx /= 0 .OR. isy /= 0 .OR. np < 1 .OR. nn < 3) THEN
         WRITE (line, '(a,a,4(1x,i0))') TRIM(path), ': C must be a full surface (ISX = ISY = 0) with panels and nodes; found', &
            np, nn, isx, isy
         msg = TRIM(line)
         CLOSE (u)
         RETURN
      END IF
      CALL allocate_mesh(cm, np, nn, [0D0, 0D0, 0D0])
      DO k = 1, 2
         READ (u, *, IOSTAT=ios)
         IF (ios /= 0) GO TO 900
      END DO
      DO k = 1, nn
         READ (u, *, IOSTAT=ios) idx, x
         IF (ios /= 0) GO TO 900
         cm%xnode(:, k) = x
      END DO
      DO k = 1, 3
         READ (u, *, IOSTAT=ios)
         IF (ios /= 0) GO TO 900
      END DO
      DO ip = 1, np
         READ (u, '(a)', IOSTAT=ios) rec
         IF (ios /= 0) GO TO 900
         READ (rec, *, IOSTAT=ios) idx, nv
         IF (ios /= 0) GO TO 900
         IF (nv /= 3 .AND. nv /= 4) THEN
            WRITE (line, '(a,a,i0,a,i0,a)') TRIM(path), ': panel ', ip, ' has ', nv, ' vertices; use 3 or 4'
            msg = TRIM(line)
            CLOSE (u)
            RETURN
         END IF
         v = 0
         READ (rec, *, IOSTAT=ios) idx, nv, v(1:nv)
         IF (ios /= 0) GO TO 900
         IF (ANY(v(1:nv) < 1) .OR. ANY(v(1:nv) > nn)) THEN
            WRITE (line, '(a,a,i0,a)') TRIM(path), ': panel ', ip, ' refers to a node that does not exist'
            msg = TRIM(line)
            CLOSE (u)
            RETURN
         END IF
         d1 = cm%xnode(:, v(3)) - cm%xnode(:, v(1))
         d2 = cm%xnode(:, v(2)) - cm%xnode(:, v(1))
         IF (nv == 4) d2 = cm%xnode(:, v(4)) - cm%xnode(:, v(2))
         IF (SUM(cross(d1, d2)**2) <= 0D0) THEN
            WRITE (line, '(a,a,i0,a)') TRIM(path), ': panel ', ip, ' has zero area'
            msg = TRIM(line)
            CLOSE (u)
            RETURN
         END IF
         k = ip - 1
         CALL add_panel(cm, k, v(1:nv))
      END DO
      CLOSE (u)
      ok = .TRUE.
      RETURN

900   CONTINUE
      CLOSE (u)
      msg = 'cannot read '//path//' as a HullMesh.pnl file'
   END SUBROUTINE read_control_mesh

   SUBROUTINE read_file_surface(path, name, c, cm, ok, msg)
      ! C of one body from the file path (called name in messages). c is the centre of the waterplane bounding box of
      ! the body. The panel normals must point towards the body: S = sum dot(n_i, xc_i - c) A_i must be negative. This
      ! is a cheap sign test, not a containment test; containment is the user's job.
      CHARACTER(LEN=*), INTENT(IN) :: path, name
      REAL(8), INTENT(IN) :: c(2)
      TYPE(panel_mesh), INTENT(OUT) :: cm
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      LOGICAL :: exists
      REAL(8) :: s, x(3)
      INTEGER :: i

      ok = .FALSE.
      INQUIRE (FILE=path, EXIST=exists)
      IF (.NOT. exists) THEN
         msg = 'cs_shape = file needs CSMesh_'//name(8:LEN(name) - 4)//'.pnl in the input directory.'
         RETURN
      END IF
      CALL read_control_mesh(path, cm, ok, msg)
      IF (.NOT. ok) RETURN
      cm%xr(:, 1) = [c(1), c(2), 0D0]
      s = 0D0
      DO i = 1, cm%n
         x = cm%xc(:, i) - cm%xr(:, 1)
         s = s + DOT_PRODUCT(cm%nrm(:, i), x)*cm%area(i)
         cm%nj(4, i) = x(2)*cm%nrm(3, i) - x(3)*cm%nrm(2, i)
         cm%nj(5, i) = x(3)*cm%nrm(1, i) - x(1)*cm%nrm(3, i)
         cm%nj(6, i) = x(1)*cm%nrm(2, i) - x(2)*cm%nrm(1, i)
      END DO
      IF (.NOT. s < 0D0) THEN
         ok = .FALSE.
         msg = name//': the panel normals must point towards the body.'
      END IF
   END SUBROUTINE read_file_surface

   REAL(8) FUNCTION median(a)
      ! Median by sorting a copy with a shell sort.
      REAL(8), INTENT(IN) :: a(:)
      REAL(8), ALLOCATABLE :: s(:)
      REAL(8) :: t
      INTEGER :: n, gap, i, j

      n = SIZE(a)
      IF (n == 0) THEN
         median = 0D0
         RETURN
      END IF
      s = a
      gap = n/2
      DO WHILE (gap > 0)
         DO i = gap + 1, n
            t = s(i)
            j = i
            DO WHILE (j > gap)
               IF (s(j - gap) <= t) EXIT
               s(j) = s(j - gap)
               j = j - gap
            END DO
            s(j) = t
         END DO
         gap = gap/2
      END DO
      IF (MOD(n, 2) == 1) THEN
         median = s((n + 1)/2)
      ELSE
         median = 0.5D0*(s(n/2) + s(n/2 + 1))
      END IF
   END FUNCTION median

   ! ------------------------------------------------------------------------------------------------ check

   SUBROUTINE check_control_surface(mesh, pts, cms, generated, ok, msg)
      ! Two checks on C against the hulls.
      !   generated  .true. for a C from cylinder_surface or box_surface: every node of another body must lie outside
      !              body b's C. C counts as containing a node when the node is inside the horizontal polygon of the
      !              Gamma_c segments of body b and above the bottom of C (or anywhere, if C goes down to the seabed).
      !              The polygon is built from the exact segment ends, so the test is exact for a cylinder or box.
      !   always     no C point (wall panel or Gamma_c segment) may lie within one local panel size of a hull panel
      !              centroid. The local size is the larger of the sizes of the two panels.
      ! ok = .false. when a check fails, and msg names the bodies and the point.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(panel_mesh), INTENT(IN) :: cms(:)
      LOGICAL, INTENT(IN) :: generated
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      REAL(8), ALLOCATABLE :: sx1(:,:), sx2(:,:)
      REAL(8) :: zbot, d2, lim, t(2), x(3)
      INTEGER :: b, a, i, m, p, nseg, k
      LOGICAL, ALLOCATABLE :: bodyline(:)
      CHARACTER(LEN=300) :: line

      ok = .TRUE.
      msg = ''
      IF (generated) THEN
         DO b = 1, MIN(SIZE(cms), mesh%nbody)
            ! Gamma_c segments of body b, from the stored midpoint, normal and length.
            bodyline = pts%kind == CS_LINE .AND. pts%body == b
            nseg = COUNT(bodyline)
            IF (nseg == 0) CYCLE
            ALLOCATE (sx1(2, nseg), sx2(2, nseg))
            k = 0
            DO p = 1, pts%n
               IF (.NOT. bodyline(p)) CYCLE
               k = k + 1
               t = [-pts%nrm(2, p), pts%nrm(1, p)]
               sx1(:, k) = pts%x(1:2, p) - 0.5D0*pts%w(p)*t
               sx2(:, k) = pts%x(1:2, p) + 0.5D0*pts%w(p)*t
            END DO
            zbot = -HUGE(1D0)
            IF (ANY(cms(b)%nrm(3, :) > 0.5D0)) zbot = MINVAL(cms(b)%xnode(3, :))
            DO i = 1, mesh%n
               a = mesh%body(i)
               IF (a == b) CYCLE
               DO m = 1, mesh%ncn(i)
                  x = mesh%xnode(:, mesh%ncon(m, i))
                  IF (x(3) < zbot) CYCLE
                  IF (inside_polygon(x(1:2), sx1, sx2)) THEN
                     WRITE (line, '(a,i0,a,i0,a,3(1x,es11.3),a)') 'the control surface of body ', b, &
                        ' encloses a node of body ', a, ' at', x, '; use a smaller cs_margin or a different cs_shape'
                     msg = TRIM(line)
                     ok = .FALSE.
                     RETURN
                  END IF
               END DO
            END DO
            DEALLOCATE (sx1, sx2)
         END DO
      END IF

      ! No C point within one local panel size of a hull panel centroid.
      DO p = 1, pts%n
         IF (pts%kind(p) /= CS_WALL .AND. pts%kind(p) /= CS_LINE) CYCLE
         DO i = 1, mesh%n
            lim = MAX(pts%psize(p), mesh%psize(i))
            d2 = SUM((pts%x(:, p) - mesh%xc(:, i))**2)
            IF (d2 < lim**2) THEN
               WRITE (line, '(a,i0,a,i0,a,i0,a,es10.3,a)') 'a point of the control surface of body ', pts%body(p), &
                  ' lies within one panel size of hull panel ', i, ' of body ', mesh%body(i), ' (distance ', SQRT(d2), &
                  '); increase cs_margin or reduce cs_panel_size'
               msg = TRIM(line)
               ok = .FALSE.
               RETURN
            END IF
         END DO
      END DO
   END SUBROUTINE check_control_surface

   LOGICAL FUNCTION inside_polygon(p, s1, s2)
      ! Crossing-number test of the point p against the closed polygon given by its segments s1(:,k) to s2(:,k), in
      ! any order and orientation.
      REAL(8), INTENT(IN) :: p(2), s1(:,:), s2(:,:)
      INTEGER :: k, cnt
      REAL(8) :: xint

      cnt = 0
      DO k = 1, SIZE(s1, 2)
         IF ((s1(2, k) > p(2)) .NEQV. (s2(2, k) > p(2))) THEN
            xint = s1(1, k) + (p(2) - s1(2, k))*(s2(1, k) - s1(1, k))/(s2(2, k) - s1(2, k))
            IF (p(1) < xint) cnt = cnt + 1
         END IF
      END DO
      inside_polygon = MOD(cnt, 2) == 1
   END FUNCTION inside_polygon

   ! ------------------------------------------------------------------------------------------------ output

   SUBROUTINE write_control_surface(outputdir, cms)
      ! Writes CSMesh_{b}.pnl per body in the HullMesh.pnl format, which read_control_mesh reads back.
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      TYPE(panel_mesh), INTENT(IN) :: cms(:)
      CHARACTER(LEN=16) :: istr
      INTEGER :: b, u, i

      DO b = 1, SIZE(cms)
         WRITE (istr, '(i0)') b
         OPEN (NEWUNIT=u, FILE=TRIM(outputdir)//'/CSMesh_'//TRIM(istr)//'.pnl', STATUS='replace', ACTION='write')
         WRITE (u, '(a)') '    --------------Hull Mesh File---------------'
         WRITE (u, '(a)') ' '
         WRITE (u, '(a)') '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry'
         WRITE (u, '(4i12)') cms(b)%n, cms(b)%nnode, 0, 0
         WRITE (u, '(a)') ' '
         WRITE (u, '(a)') '    # Start Definition of Node Coordinates     ! node_number   x   y   z'
         DO i = 1, cms(b)%nnode
            WRITE (u, '(i5,3f20.12)') i, cms(b)%xnode(:, i)
         END DO
         WRITE (u, '(a)') '    # End Definition of Node Coordinates'
         WRITE (u, '(a)') ' '
         WRITE (u, '(a)') '  # Start Definition of Node Relations   ! panel_number  number_of_vertices   ' &
                          //'Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID)'
         DO i = 1, cms(b)%n
            WRITE (u, '(i5,i9,4i14)') i, cms(b)%ncn(i), cms(b)%ncon(1:cms(b)%ncn(i), i)
         END DO
         WRITE (u, '(a)') '    # End Definition of Node Relations'
         WRITE (u, '(a)') ' '
         WRITE (u, '(a)') '    --------------End Hull Mesh File---------------'
         CLOSE (u)
      END DO
   END SUBROUTINE write_control_surface

END MODULE ControlSurfaceMesh
