"""
Generate a HAMS hull mesh and a matching Hydrostatic.in for a surface-piercing vertical cylinder or box.

Usage:
    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --out HullMesh.pnl
    python make_mesh.py box --lx 20 --ly 10 --draft 5 --nx 10 --ny 5 --nz 4 --out HullMesh.pnl
    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --format nemoh --out hull.dat
    python make_mesh.py oc4 --multiple 1 [--out-dir folder]

The hull has no lid, so the top row of nodes lies on z = 0. --x0 and --y0 shift the body in the horizontal plane. Every
generated panel normal is checked to point out of the fluid (into the body), computed as HAMS does in CalTransNormals.
--half-y keeps the panels with y > 0 and sets ISY = 1, for a body symmetric about y = 0.
--waterplane-out also writes the interior free surface of a cylinder, for irregular frequency removal (IRSP = 1).

Hydrostatic.in is written next to the mesh (or to --hydro-out) for a freely floating, uniform solid prism, using the
discretized waterplane polygon: mass = rho*V, centre of gravity at half draft, all matrices about the point (x0, y0, 0),
which must be given as the rotation centre in ControlFile.in. External damping and stiffness are zero. These are the
values of the whole body, also with --half-y. write_body and write_cylinder also accept external restoring values, and
write_submerged_cylinder writes a closed body below the free surface, for use from make_drift_fixtures.py.
write_rounded_box writes a box with rounded vertical corners and, if wanted, a rounded bilge, for drift_convergence.py.

oc4 writes the OC4-DeepCwind semisubmersible (columns only, see OC4_* below): HullMesh.pnl, Hydrostatic.in with the mooring
as the external restoring, and the control surface C and free surface F as CSMesh_1.pnl and FSMesh_1.pnl for cs_shape = file,
then prints the checks against the exact columns and the report. Without --out-dir the files go to a temporary folder.

--format nemoh writes the same nodes and panels in the mesh format of Nemoh v3. The vertex order is kept: it gives a right-hand
normal into the fluid, as Nemoh documents, while HAMS computes the opposite normal from the same order. A triangle has its last
node repeated, since that format only has four node panels. Hydrostatic.in is still written in the HAMS format.
"""

import argparse
import math
import os

RHO = 1000.0     # as Const_mod in src/WavDynMods.f90
G = 9.80665
NODE_FORMAT = "%.10f"   # precision of the node coordinates in the mesh file
HYDROSTATIC_STEPS = 2000   # sections of a body that is not a prism


class Mesh:
    def __init__(self):
        self.nodes = []
        self.index = {}
        self.panels = []

    def node(self, x, y, z):
        """Return the 1-based index of a node, sharing nodes that coincide."""
        key = (round(x, 9), round(y, 9), round(z, 9))
        if key not in self.index:
            self.nodes.append((x, y, z))
            self.index[key] = len(self.nodes)
        return self.index[key]

    def add_panel(self, points):
        """Add a panel from 3 or 4 points, ordered counter-clockwise as seen from the fluid."""
        self.panels.append([self.node(*p) for p in points])

    def shifted(self, x0, y0):
        out = Mesh()
        out.nodes = [(x + x0, y + y0, z) for x, y, z in self.nodes]
        out.panels = self.panels
        return out


def half_in_y(mesh):
    """The panels with y > 0 of a mesh symmetric about y = 0, which must not have panels crossing that plane."""
    out = Mesh()
    for panel in mesh.panels:
        points = [mesh.nodes[i - 1] for i in panel]
        if min(p[1] for p in points) < -1e-9 and max(p[1] for p in points) > 1e-9:
            raise ValueError("a panel crosses y = 0, so the mesh cannot be halved")
        if max(p[1] for p in points) > 1e-9:
            out.add_panel(points)
    return out


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def hams_normal(points):
    """Unit normal with the vertex order used by CalTransNormals in src/NormalProcess.f90."""
    if len(points) == 3:
        un = cross(sub(points[0], points[1]), sub(points[2], points[1]))
    else:
        un = cross(sub(points[3], points[1]), sub(points[2], points[0]))
    length = math.sqrt(sum(c * c for c in un))
    if length < 1e-12:
        raise ValueError("degenerate panel: " + str(points))
    return tuple(c / length for c in un)


def check_normals(mesh, inside):
    """Stop unless every panel normal points from the panel towards the point inside the body."""
    for number, panel in enumerate(mesh.panels, start=1):
        points = [mesh.nodes[i - 1] for i in panel]
        centre = tuple(sum(p[k] for p in points) / len(points) for k in range(3))
        normal = hams_normal(points)
        if sum(n * d for n, d in zip(normal, sub(inside, centre))) <= 0.0:
            raise ValueError("panel %d has a normal pointing into the fluid" % number)


def cylinder(radius, draft, ntheta, nz, nr, graded=False, edge_power=None):
    """Side wall and flat bottom, with triangles at the centre of the bottom. Nodes lie on y = 0, so the mesh is
    symmetric about y = 0. With graded, rows are cosine spaced, finer at the waterline and the bilge. With edge_power q
    (instead of graded), the rows are power-law graded: on the side wall s^q / (s^q + (1 - s)^q) of the draft, finer at
    the waterline and the bilge; on the bottom the distance from the bilge edge is s^q of the radius. The first row is then
    about (1/n)^q of the draft or radius, against about 2.5/n^2 (side) and 1.2/n^2 (bottom) for cosine spacing.
    Returns the mesh and its waterline polygon, counter-clockwise from above."""
    if edge_power is not None and not edge_power >= 1.0:
        raise ValueError("edge_power must be at least 1")
    if ntheta % 2 != 0:
        raise ValueError("ntheta must be even, to keep the mesh symmetric about both axes")
    mesh = Mesh()
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta + 1)]
    angles[-1] = 0.0

    def ring(r, j, z):
        return (r * math.cos(angles[j]), r * math.sin(angles[j]), z)

    def depth(i):
        s = i / nz
        if edge_power is not None:
            return -draft * s ** edge_power / (s ** edge_power + (1.0 - s) ** edge_power)
        return -draft * (0.5 * (1.0 - math.cos(math.pi * s)) if graded else s)

    def bottom_radius(i):
        s = i / nr
        if edge_power is not None:
            return radius * (1.0 - s ** edge_power)
        return radius * (math.cos(0.5 * math.pi * s) if graded else 1.0 - s)

    for j in range(ntheta):
        for i in range(nz):
            upper = depth(i)
            lower = depth(i + 1)
            mesh.add_panel([ring(radius, j, lower), ring(radius, j + 1, lower),
                            ring(radius, j + 1, upper), ring(radius, j, upper)])
        for i in range(nr):
            outer = bottom_radius(i)
            inner = bottom_radius(i + 1)
            if i < nr - 1:
                mesh.add_panel([ring(outer, j + 1, -draft), ring(outer, j, -draft),
                                ring(inner, j, -draft), ring(inner, j + 1, -draft)])
            else:
                mesh.add_panel([ring(outer, j + 1, -draft), ring(outer, j, -draft), (0.0, 0.0, -draft)])
    polygon = [(radius * math.cos(angles[j]), radius * math.sin(angles[j])) for j in range(ntheta)]
    return mesh, polygon


def waterplane_disc(radius, ntheta, nr, graded=False):
    """Interior free surface of a cylinder: rings at z = 0 with triangles at the centre. Normals point down, as in the
    waterplane meshes of examples/. HAMS needs z = 0 exactly here."""
    mesh = Mesh()
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta)] + [0.0]

    def ring(r, j):
        return (r * math.cos(angles[j]), r * math.sin(angles[j]), 0.0)

    def ring_radius(i):
        return radius * (math.cos(0.5 * math.pi * i / nr) if graded else (nr - i) / nr)

    for j in range(ntheta):
        for i in range(nr):
            outer = ring_radius(i)
            inner = ring_radius(i + 1)
            if i < nr - 1:
                mesh.add_panel([ring(inner, j + 1), ring(inner, j), ring(outer, j), ring(outer, j + 1)])
            else:
                mesh.add_panel([(0.0, 0.0, 0.0), ring(outer, j), ring(outer, j + 1)])
    return mesh


def revolved(profile, ntheta):
    """Surface of revolution of the profile (r, z), from the waterline down to the axis (r = 0). Returns the mesh and its
    waterline polygon, counter-clockwise from above."""
    if ntheta % 2 != 0:
        raise ValueError("ntheta must be even, to keep the mesh symmetric about both axes")
    mesh = Mesh()
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta)] + [0.0]

    def point(p, j):
        return (p[0] * math.cos(angles[j]), p[0] * math.sin(angles[j]), p[1])

    for j in range(ntheta):
        for i in range(len(profile) - 1):
            a, b = profile[i], profile[i + 1]
            if b[0] > 1e-12:
                mesh.add_panel([point(b, j), point(b, j + 1), point(a, j + 1), point(a, j)])
            else:
                mesh.add_panel([(0.0, 0.0, b[1]), point(a, j + 1), point(a, j)])
    polygon = [(profile[0][0] * math.cos(angles[j]), profile[0][0] * math.sin(angles[j])) for j in range(ntheta)]
    return mesh, polygon


def rounded_profile(radius, draft, bilge_radius, nz, narc, nr, uniform_rows=False):
    """Wall with rows cosine spaced towards the waterline, or uniform in z with uniform_rows, a circular bilge, and a flat
    bottom, as (r, z)."""
    wall = draft - bilge_radius
    out = [(radius, -wall * (i / nz if uniform_rows else 1.0 - math.cos(0.5 * math.pi * i / nz))) for i in range(nz + 1)]
    for i in range(1, narc + 1):
        angle = 0.5 * math.pi * i / narc
        out.append((radius - bilge_radius + bilge_radius * math.cos(angle), -wall - bilge_radius * math.sin(angle)))
    out += [((radius - bilge_radius) * (1.0 - i / nr), -draft) for i in range(1, nr + 1)]
    return out


def rounded_radius(radius, draft, bilge_radius, z):
    """Hull radius at depth z of the rounded cylinder."""
    below = -z - (draft - bilge_radius)
    if below <= 0.0:
        return radius
    return radius - bilge_radius + math.sqrt(max(bilge_radius ** 2 - below ** 2, 0.0))


def half_ellipsoid(ax, ay, draft, ntheta, nphi):
    """Lower half of an ellipsoid with semi-axes ax, ay and draft, and rows of equal polar angle. It has no edges, and
    at an oblique heading it has a yaw moment of its own. Returns the mesh and its waterline polygon."""
    mesh, polygon = revolved([(math.cos(0.5 * math.pi * i / nphi), -math.sin(0.5 * math.pi * i / nphi) if i else 0.0)
                              for i in range(nphi)] + [(0.0, -1.0)], ntheta)
    mesh.nodes = [(ax * x, ay * y, draft * z) for x, y, z in mesh.nodes]
    return mesh, [(ax * x, ay * y) for x, y in polygon]


def closed_cylinder(radius, top, height, ntheta, nz, nr):
    """Side wall, top and bottom of a cylinder entirely below the free surface, with its top at z = -top."""
    mesh = Mesh()
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta)] + [0.0]

    def ring(r, j, z):
        return (r * math.cos(angles[j]), r * math.sin(angles[j]), z)

    for j in range(ntheta):
        for i in range(nz):
            upper = -top - height * i / nz
            lower = -top - height * (i + 1) / nz
            mesh.add_panel([ring(radius, j, lower), ring(radius, j + 1, lower),
                            ring(radius, j + 1, upper), ring(radius, j, upper)])
        for i in range(nr):
            outer = radius * (nr - i) / nr
            inner = radius * (nr - i - 1) / nr
            if i < nr - 1:
                mesh.add_panel([ring(outer, j + 1, -top - height), ring(outer, j, -top - height),
                                ring(inner, j, -top - height), ring(inner, j + 1, -top - height)])
                mesh.add_panel([ring(outer, j, -top), ring(outer, j + 1, -top),
                                ring(inner, j + 1, -top), ring(inner, j, -top)])
            else:
                mesh.add_panel([ring(outer, j + 1, -top - height), ring(outer, j, -top - height),
                                (0.0, 0.0, -top - height)])
                mesh.add_panel([ring(outer, j, -top), ring(outer, j + 1, -top), (0.0, 0.0, -top)])
    return mesh


def box(lx, ly, draft, nx, ny, nz):
    """Four walls and a flat bottom, centred on the origin. Returns the mesh and its waterline polygon."""
    mesh = Mesh()
    corners = [(-lx / 2, -ly / 2), (lx / 2, -ly / 2), (lx / 2, ly / 2), (-lx / 2, ly / 2)]
    for c in range(4):
        a, b = corners[c], corners[(c + 1) % 4]
        count = nx if c % 2 == 0 else ny

        def wall(s, z):
            return (a[0] + (b[0] - a[0]) * s / count, a[1] + (b[1] - a[1]) * s / count, z)

        for s in range(count):
            for i in range(nz):
                upper = -draft * i / nz
                lower = -draft * (i + 1) / nz
                mesh.add_panel([wall(s, lower), wall(s + 1, lower), wall(s + 1, upper), wall(s, upper)])

    def bottom(ix, iy):
        return (-lx / 2 + lx * ix / nx, -ly / 2 + ly * iy / ny, -draft)

    for ix in range(nx):
        for iy in range(ny):
            mesh.add_panel([bottom(ix, iy), bottom(ix, iy + 1), bottom(ix + 1, iy + 1), bottom(ix + 1, iy)])
    return mesh, corners


def rounded_rectangle(lx, ly, corner_radius, nx, ny, ncorner, inset=0.0):
    """Closed polygon of a rectangle with rounded corners, moved inwards by inset, counter-clockwise from above. It starts in
    the middle of the corner at (-x, -y). It has nx and ny segments on the straight sides and ncorner on each corner, so
    the sides between corner middles have nx + ncorner and ny + ncorner segments."""
    half_x = lx / 2.0 - corner_radius
    half_y = ly / 2.0 - corner_radius
    radius = corner_radius - inset
    centres = [(-half_x, -half_y), (half_x, -half_y), (half_x, half_y), (-half_x, half_y)]

    def arc(c, i):
        angle = math.pi + 0.5 * math.pi * (c + i / ncorner)
        return (centres[c][0] + radius * math.cos(angle), centres[c][1] + radius * math.sin(angle))

    points = []
    for c in range(4):
        points += [arc(c, i) for i in range(ncorner // 2 if c == 0 else 0, ncorner + 1)]
        ax, ay = arc(c, ncorner)
        bx, by = arc((c + 1) % 4, 0)
        count = nx if c % 2 == 0 else ny
        points += [(ax + (bx - ax) * i / count, ay + (by - ay) * i / count) for i in range(1, count)]
    return points + [arc(0, i) for i in range(ncorner // 2)]


def rounded_box(lx, ly, draft, corner_radius, bilge_radius, nx, ny, ncorner, nz, narc):
    """Walls with rounded vertical corners, rows cosine spaced towards the waterline, a rounded bilge if bilge_radius is not
    zero, and a flat bottom of quadrilaterals. ncorner must be even. Returns the mesh and its waterline polygon."""
    if ncorner % 2 != 0:
        raise ValueError("ncorner must be even, since the sides of the bottom meet in the middle of a corner")
    if bilge_radius >= corner_radius:
        raise ValueError("the bilge radius must be less than the corner radius, to keep the corners of the bottom rounded")
    profile = rounded_profile(bilge_radius, draft, bilge_radius, nz, narc, 0)
    rings = [rounded_rectangle(lx, ly, corner_radius, nx, ny, ncorner, bilge_radius - r) for r, _ in profile]
    mesh = Mesh()
    count = len(rings[0])
    for level in range(len(profile) - 1):
        upper, lower = rings[level], rings[level + 1]
        z_upper, z_lower = profile[level][1], profile[level + 1][1]
        for j in range(count):
            k = (j + 1) % count
            mesh.add_panel([lower[j] + (z_lower,), lower[k] + (z_lower,), upper[k] + (z_upper,), upper[j] + (z_upper,)])

    # The bottom is a transfinite grid between the four sides of the last ring
    edge = rings[-1]
    side_x = nx + ncorner
    side_y = ny + ncorner
    south = edge[:side_x + 1]
    east = edge[side_x:side_x + side_y + 1]
    north = edge[side_x + side_y:2 * side_x + side_y + 1][::-1]
    west = (edge[2 * side_x + side_y:] + edge[:1])[::-1]

    def bottom(i, j):
        if j == 0:
            return south[i] + (-draft,)
        if j == side_y:
            return north[i] + (-draft,)
        if i == 0:
            return west[j] + (-draft,)
        if i == side_x:
            return east[j] + (-draft,)
        s, t = i / side_x, j / side_y
        return tuple((1 - t) * south[i][k] + t * north[i][k] + (1 - s) * west[j][k] + s * east[j][k]
                     - (1 - s) * (1 - t) * south[0][k] - s * (1 - t) * south[-1][k]
                     - s * t * north[-1][k] - (1 - s) * t * north[0][k] for k in range(2)) + (-draft,)

    for i in range(side_x):
        for j in range(side_y):
            mesh.add_panel([bottom(i, j), bottom(i, j + 1), bottom(i + 1, j + 1), bottom(i + 1, j)])
    return mesh, rings[0]


def polygon_moments(polygon):
    """Area and second moments int y^2 dA, int x^2 dA of a counter-clockwise polygon centred on the origin."""
    area = ixx = iyy = 0.0
    for k, (x1, y1) in enumerate(polygon):
        x2, y2 = polygon[(k + 1) % len(polygon)]
        c = x1 * y2 - x2 * y1
        area += 0.5 * c
        ixx += c * (y1 * y1 + y1 * y2 + y2 * y2) / 12.0
        iyy += c * (x1 * x1 + x1 * x2 + x2 * x2) / 12.0
    return area, ixx, iyy


def hydrostatics(polygon, draft, scale=lambda z: 1.0, steps=1, section=None):
    """Mass and restoring matrices about the waterplane centre of a uniform solid floating at its draft, whose sections
    are the waterline polygon scaled by scale(z), or the polygon section(z) if that is given. The default is a prism, for
    which one step of the rule is exact."""
    area, ixx, iyy = polygon_moments(polygon)
    vol = zmoment = i44 = i55 = i66 = 0.0
    dz = draft / steps
    for step in range(steps):
        # Exact sections are integrated by the midpoint rule, or for a prism by its exact second moment of z
        z = -(step + 0.5) * dz
        f = scale(z)
        z2 = z * z + dz * dz / 12.0
        a, xx, yy = polygon_moments(section(z)) if section else (area, ixx, iyy)
        vol += f ** 2 * a * dz
        zmoment += f ** 2 * a * z * dz
        i44 += RHO * (f ** 4 * xx + f ** 2 * a * z2) * dz
        i55 += RHO * (f ** 4 * yy + f ** 2 * a * z2) * dz
        i66 += RHO * f ** 4 * (xx + yy) * dz
    mass = RHO * vol
    zg = zb = zmoment / vol

    m = [[0.0] * 6 for _ in range(6)]
    for k in range(3):
        m[k][k] = mass
    m[0][4] = m[4][0] = mass * zg
    m[1][3] = m[3][1] = -mass * zg
    m[3][3], m[4][4], m[5][5] = i44, i55, i66

    c = [[0.0] * 6 for _ in range(6)]
    c[2][2] = RHO * G * area
    c[3][3] = RHO * G * (ixx + vol * zb) - mass * G * zg
    c[4][4] = RHO * G * (iyy + vol * zb) - mass * G * zg
    if c[3][3] <= 0.0 or c[4][4] <= 0.0:
        raise ValueError("the body is unstable: a metacentric height is not positive")
    return {"vol": vol, "area": area, "zg": zg, "zb": zb, "mass": m, "restoring": c}


def write_mesh(path, mesh, isy=0, title="Hull Mesh File"):
    with open(path, "w") as f:
        f.write("    --------------%s---------------\n\n" % title)
        f.write("    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n")
        f.write("    %8d  %8d  %8d  %8d\n\n" % (len(mesh.panels), len(mesh.nodes), 0, isy))
        f.write("    # Start Definition of Node Coordinates     ! node_number   x   y   z\n")
        for number, (x, y, z) in enumerate(mesh.nodes, start=1):
            f.write("    %6d  " % number + "  ".join((NODE_FORMAT % c).rjust(18) for c in (x, y, z)) + "\n")
        f.write("    # End Definition of Node Coordinates\n\n")
        f.write("    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ...\n")
        for number, panel in enumerate(mesh.panels, start=1):
            f.write("    %6d  %3d  " % (number, len(panel)) + "  ".join("%6d" % i for i in panel) + "\n")
        f.write("    # End Definition of Node Relations\n\n")
        f.write("    --------------End %s---------------\n" % title)


def write_nemoh_mesh(path, mesh, isym=0):
    """Write a mesh in the format of Nemoh v3. A triangle has its last node repeated."""
    with open(path, "w") as f:
        f.write("2 %d\n" % isym)
        for number, (x, y, z) in enumerate(mesh.nodes, start=1):
            f.write("%d " % number + " ".join(NODE_FORMAT % c for c in (x, y, z)) + "\n")
        f.write("0 0. 0. 0.\n")
        for panel in mesh.panels:
            f.write(" ".join("%d" % i for i in (panel + panel[-1:])[:4]) + "\n")
        f.write("0 0 0 0\n")


def write_matrix(f, matrix):
    for row in matrix:
        f.write("".join("  %12.5E" % value for value in row) + "\n")


def write_hydrostatic(path, hydro, x0, y0, external_restoring=None, external_damping=None):
    """external_restoring and external_damping map a 0-based mode number to a value on the diagonal of the external
    restoring and the external linear damping matrix. external_restoring may also be a full 6 x 6 matrix (rows of lists)."""
    zero = [[0.0] * 6 for _ in range(6)]
    if isinstance(external_restoring, (list, tuple)):
        external = [list(row) for row in external_restoring]
        external_restoring = None
    else:
        external = [[0.0] * 6 for _ in range(6)]
    for mode, value in (external_restoring or {}).items():
        external[mode][mode] = value
    damping = [[0.0] * 6 for _ in range(6)]
    for mode, value in (external_damping or {}).items():
        damping[mode][mode] = value
    with open(path, "w") as f:
        f.write(" Center of Gravity:\n")
        f.write("  %.12E  %.12E  %.12E\n" % (x0, y0, hydro["zg"]))
        f.write(" Body Mass Matrix:\n")
        write_matrix(f, hydro["mass"])
        f.write(" External Linear Damping Matrix:\n")
        write_matrix(f, damping)
        f.write(" External Quadratic Damping Matrix:\n")
        write_matrix(f, zero)
        f.write(" Hydrostatic Restoring Matrix:\n")
        write_matrix(f, hydro["restoring"])
        f.write(" External Restoring Matrix:\n")
        write_matrix(f, external)


def write_body(mesh, polygon, draft, x0, y0, half_y, mesh_path, hydro_path, external_restoring=None, mesh_format="hams",
              external_damping=None):
    """Check and write a mesh centred on the origin, and its hydrostatics. Returns the written mesh and hydrostatics."""
    if half_y:
        if y0 != 0.0:
            raise ValueError("a half mesh must be centred on y = 0")
        mesh = half_in_y(mesh)
    mesh = mesh.shifted(x0, y0)
    check_normals(mesh, (x0, y0 + (draft / 2.0 if half_y else 0.0), -draft / 2.0))
    hydro = hydrostatics(polygon, draft)
    if mesh_format == "nemoh":
        write_nemoh_mesh(mesh_path, mesh, 1 if half_y else 0)
    else:
        write_mesh(mesh_path, mesh, 1 if half_y else 0)
    write_hydrostatic(hydro_path, hydro, x0, y0, external_restoring, external_damping)
    return mesh, hydro


def write_cylinder(mesh_path, hydro_path, radius, draft, ntheta, nz, nr, graded=False, x0=0.0, y0=0.0, half_y=False,
                   external_restoring=None, edge_power=None, external_damping=None):
    """Write the mesh and hydrostatics of a cylinder. Returns the written mesh."""
    mesh, polygon = cylinder(radius, draft, ntheta, nz, nr, graded, edge_power)
    return write_body(mesh, polygon, draft, x0, y0, half_y, mesh_path, hydro_path, external_restoring,
                      external_damping=external_damping)[0]


def write_rounded_cylinder(mesh_path, hydro_path, radius, draft, bilge_radius, ntheta, nz, narc, nr, external_restoring=None,
                           uniform_rows=False, half_y=False):
    """Write the mesh and hydrostatics of a cylinder with a rounded bilge, which has no edge. Returns the written mesh. With
    half_y, the mesh is the half y > 0 with a y symmetry plane; the hydrostatics stay those of the whole body."""
    mesh, polygon = revolved(rounded_profile(radius, draft, bilge_radius, nz, narc, nr, uniform_rows), ntheta)
    scale = lambda z: rounded_radius(radius, draft, bilge_radius, z) / radius
    if half_y:
        mesh = half_in_y(mesh)
    check_normals(mesh, (0.0, radius / 2.0 if half_y else 0.0, -draft / 2.0))
    write_mesh(mesh_path, mesh, 1 if half_y else 0)
    write_hydrostatic(hydro_path, hydrostatics(polygon, draft, scale, HYDROSTATIC_STEPS), 0.0, 0.0, external_restoring)
    return mesh


def write_half_ellipsoid(mesh_path, hydro_path, ax, ay, draft, ntheta, nphi):
    """Write the mesh and hydrostatics of a half ellipsoid. Returns the written mesh."""
    mesh, polygon = half_ellipsoid(ax, ay, draft, ntheta, nphi)
    scale = lambda z: math.sqrt(max(1.0 - (z / draft) ** 2, 0.0))
    check_normals(mesh, (0.0, 0.0, -draft / 2.0))
    write_mesh(mesh_path, mesh)
    write_hydrostatic(hydro_path, hydrostatics(polygon, draft, scale, HYDROSTATIC_STEPS), 0.0, 0.0)
    return mesh


def write_box(mesh_path, hydro_path, lx, ly, draft, nx, ny, nz, x0=0.0, y0=0.0):
    """Write the mesh and hydrostatics of a box. Returns the written mesh."""
    mesh, polygon = box(lx, ly, draft, nx, ny, nz)
    return write_body(mesh, polygon, draft, x0, y0, False, mesh_path, hydro_path)[0]


def write_rounded_box(mesh_path, hydro_path, lx, ly, draft, corner_radius, bilge_radius, nx, ny, ncorner, nz, narc):
    """Write the mesh and hydrostatics of a box with rounded corners. Returns the written mesh."""
    mesh, polygon = rounded_box(lx, ly, draft, corner_radius, bilge_radius, nx, ny, ncorner, nz, narc)

    def section(z):
        inset = bilge_radius - rounded_radius(bilge_radius, draft, bilge_radius, z) if bilge_radius > 0.0 else 0.0
        return rounded_rectangle(lx, ly, corner_radius, nx, ny, ncorner, inset)

    check_normals(mesh, (0.0, 0.0, -draft / 2.0))
    write_mesh(mesh_path, mesh)
    write_hydrostatic(hydro_path, hydrostatics(polygon, draft, steps=HYDROSTATIC_STEPS, section=section), 0.0, 0.0)
    return mesh


def write_waterplane(path, radius, ntheta, nr, graded=False, x0=0.0, y0=0.0):
    """Write the interior free surface mesh of a cylinder, checking that every normal points down."""
    mesh = waterplane_disc(radius, ntheta, nr, graded).shifted(x0, y0)
    for number, panel in enumerate(mesh.panels, start=1):
        if hams_normal([mesh.nodes[i - 1] for i in panel])[2] >= 0.0:
            raise ValueError("waterplane panel %d has a normal pointing up" % number)
    write_mesh(path, mesh, title="Waterplane Mesh File")
    return mesh


def write_submerged_cylinder(mesh_path, hydro_path, radius, top, height, ntheta, nz, nr, x0=0.0, y0=0.0):
    """Write a closed cylinder below the free surface and the hydrostatics of a neutrally buoyant body with its centre
    of gravity below the centre of buoyancy, about (x0, y0, 0). Only the roll and pitch restoring is not zero."""
    mesh = closed_cylinder(radius, top, height, ntheta, nz, nr).shifted(x0, y0)
    check_normals(mesh, (x0, y0, -top - height / 2.0))
    polygon = [(radius * math.cos(2.0 * math.pi * j / ntheta), radius * math.sin(2.0 * math.pi * j / ntheta))
               for j in range(ntheta)]
    area, ixx, iyy = polygon_moments(polygon)
    vol = area * height
    mass = RHO * vol
    zb = -top - height / 2.0
    zg = zb - height / 4.0
    gyration = 0.25 * radius * radius + height * height / 12.0   # uniform solid about a horizontal axis through its centre
    m = [[0.0] * 6 for _ in range(6)]
    for k in range(3):
        m[k][k] = mass
    m[0][4] = m[4][0] = mass * zg
    m[1][3] = m[3][1] = -mass * zg
    m[3][3] = m[4][4] = mass * (gyration + zg * zg)
    m[5][5] = RHO * height * (ixx + iyy)
    c = [[0.0] * 6 for _ in range(6)]
    c[3][3] = c[4][4] = RHO * G * vol * zb - mass * G * zg
    write_mesh(mesh_path, mesh)
    write_hydrostatic(hydro_path, {"zg": zg, "mass": m, "restoring": c}, x0, y0)
    return mesh


# ------------------------------------------------------------------------------------------------ OC4 semisubmersible
# The OC4-DeepCwind semisubmersible of Robertson et al. (2014), Definition of the Semisubmersible Floating System for Phase
# II of OC4, NREL/TP-5000-60601: geometry from Tables 3-1 and 3-2, platform mass, CM and inertias from Table 3-3 (the
# platform alone, without tower and turbine: the comparison model), check values from Eqs. 4-2 and 4-3, linearized mooring
# from Eq. 5-15. Only the columns are panelled: the pontoons and cross braces (D = 1.6 m) are left out, as in the mesh of the
# Nemoh 3.0.3 OC4 test case (no node there is farther than 12.5 m from a column axis). The volume and the restoring of the
# report include the braces, so the mesh values differ from them by the braces (oc4_check prints both).
OC4_RHO = 1025.0          # water density of the report; HAMS has RHO = 1000 built in (oc4_matrices)
OC4_DEPTH = 200.0
OC4_DRAFT = 20.0
OC4_BASE_TOP = 14.0       # depth of the top of the base columns, where the upper columns end
OC4_MAIN_RADIUS = 3.25    # main column, D = 6.5 m, on the axis
OC4_UPPER_RADIUS = 6.0    # upper (offset) columns, D = 12 m
OC4_BASE_RADIUS = 12.0    # base columns, D = 24 m, below the upper columns
OC4_OFFSET_AXES = ((14.43, 25.0), (-28.87, 0.0), (14.43, -25.0))
OC4_MASS = 1.3473e7       # kg
OC4_ZG = -13.46           # m, CM on the centreline
OC4_I_ROLL = 6.827e9      # kg m^2 about the CM, also pitch
OC4_I_YAW = 1.226e10
OC4_REPORT = dict(rho_g_v0=1.3989e8, c33=3.836e6, c44=-3.776e8)   # Eqs. 4-2, 4-3; C44 = C55 without gravity
# Eq. 5-15, about the origin on the SWL, 0-based (surge, sway, heave, roll, pitch, yaw). The extracted text has the
# magnitudes of the couplings by their units: 1.08e5 N/rad in rows 1 and 2 (C15, C24), 1.07e5 N m/m in rows 4 and 5 (C42,
# C51). The signs are [derived]: a fairlead 14 m below the origin gives C15 = C51 < 0 and C24 = C42 = -C15 > 0; see
# oc4_mooring_signs.
OC4_MOORING = [[7.08e4, 0.0, 0.0, 0.0, -1.08e5, 0.0],
               [0.0, 7.08e4, 0.0, 1.08e5, 0.0, 0.0],
               [0.0, 0.0, 1.91e4, 0.0, 0.0, 0.0],
               [0.0, 1.07e5, 0.0, 8.73e7, 0.0, 0.0],
               [-1.07e5, 0.0, 0.0, 0.0, 8.73e7, 0.0],
               [0.0, 0.0, 0.0, 0.0, 0.0, 1.17e8]]
# Panels at multiple 1, all counts times the multiple: node counts around the main, upper and base columns; rows down the main
# column, the upper columns and the base columns (cosine spaced, finer at the waterline and the edges); rows across the
# annulus on top of a base column, of the base bottom and of the main bottom (finer at the outer edge)
OC4_MESH = dict(n_main=12, n_upper=16, n_base=32, nz_main=8, nz_upper=6, nz_base=3, nr_annulus=2, nr_base=3, nr_main=2)
# Control surface C around the whole platform (a vertical cylinder with a bottom) and its panel size; the same at every
# hull mesh. F has panels of about this size away from the columns
OC4_CS = dict(margin=8.0, bottom_margin=8.0, panel_size=4.0)
RING_HALVING = 0.6        # a ring of a disc takes half the nodes of the ring outside it when its arc is below 0.6 of the edge arc
FS_GROWTH = 1.3           # growth of the point spacing of F per ring away from a waterline


def _cosine(a, b, n):
    """n + 1 levels from a to b, cosine spaced: finer at both ends."""
    return [a + (b - a) * 0.5 * (1.0 - math.cos(math.pi * i / n)) for i in range(n + 1)]


def _ring_point(cx, cy, r, n, k, z, offset=0.0):
    angle = offset + 2.0 * math.pi * (k % n) / n
    return (cx + r * math.cos(angle), cy + r * math.sin(angle), z)


def _add_oriented(mesh, points, want):
    """Add a panel with its HAMS normal along want (reversing the vertex order if needed)."""
    if sum(a * b for a, b in zip(hams_normal(points), want)) < 0.0:
        points = points[::-1]
    mesh.add_panel(points)


def ring_counts(radii, n_outer):
    """Node counts of the rings of a disc at the radii (outer first, the last may be 0, the centre): a ring takes half the count
    of the ring outside it when its arc would fall below RING_HALVING of the arc at the edge, while that stays even and at
    least 6."""
    edge = 2.0 * math.pi * radii[0] / n_outer
    counts = [n_outer]
    for r in radii[1:]:
        count = counts[-1]
        if r > 0.0 and count % 2 == 0 and count // 2 >= 6 and 2.0 * math.pi * r / count < RING_HALVING * edge:
            count //= 2
        counts.append(count)
    return counts


def _disc(mesh, cx, cy, radii, counts, z, up):
    """Horizontal disc or annulus at z between rings at radii (outer first) with the node counts given: quadrilaterals between
    rings of the same count, three triangles per sector where the count halves (symmetric about each sector's middle), and a
    fan to the centre if the last radius is 0. HAMS normals +z if up, else -z."""
    want = (0.0, 0.0, 1.0 if up else -1.0)
    for i in range(len(radii) - 1):
        ra, rb, na, nb = radii[i], radii[i + 1], counts[i], counts[i + 1]

        def a(k):
            return _ring_point(cx, cy, ra, na, k, z)

        def b(k):
            return _ring_point(cx, cy, rb, nb, k, z)

        if rb <= 0.0:
            for k in range(na):
                _add_oriented(mesh, [a(k), a(k + 1), (cx, cy, z)], want)
        elif na == nb:
            for k in range(na):
                _add_oriented(mesh, [a(k), a(k + 1), b(k + 1), b(k)], want)
        elif na == 2 * nb:
            for k in range(nb):
                _add_oriented(mesh, [b(k), a(2 * k), a(2 * k + 1)], want)
                _add_oriented(mesh, [b(k), a(2 * k + 1), b(k + 1)], want)
                _add_oriented(mesh, [b(k + 1), a(2 * k + 1), a(2 * k + 2)], want)
        else:
            raise ValueError("ring counts %d and %d do not match" % (na, nb))


def _wall(mesh, cx, cy, r, n, levels, towards_axis=True):
    """Vertical wall of radius r about (cx, cy) between the z levels (top first). HAMS normals towards the axis, or away."""
    for k in range(n):
        for i in range(len(levels) - 1):
            points = [_ring_point(cx, cy, r, n, k, levels[i + 1]), _ring_point(cx, cy, r, n, k + 1, levels[i + 1]),
                      _ring_point(cx, cy, r, n, k + 1, levels[i]), _ring_point(cx, cy, r, n, k, levels[i])]
            centre = tuple(sum(p[c] for p in points) / 4.0 for c in range(3))
            sign = 1.0 if towards_axis else -1.0
            _add_oriented(mesh, points, (sign * (cx - centre[0]), sign * (cy - centre[1]), 0.0))


def oc4_hull(multiple=1, mesh=None):
    """The wetted columns of OC4, full (no symmetry). Returns the mesh, the waterline polygons (one per column,
    counter-clockwise from above, main column first) and, per panel, a point inside the body next to it (for the normal check).
    """
    p = dict(OC4_MESH if mesh is None else mesh)
    p = {key: int(round(value * multiple)) for key, value in p.items()}
    out = Mesh()
    inside = []

    def track(start, point_of):
        for number in range(start, len(out.panels)):
            points = [out.nodes[i - 1] for i in out.panels[number]]
            inside.append(point_of(tuple(sum(q[c] for q in points) / len(points) for c in range(3))))

    # Main column: wall from the waterline to the draft and the bottom disc
    start = len(out.panels)
    _wall(out, 0.0, 0.0, OC4_MAIN_RADIUS, p["n_main"], _cosine(0.0, -OC4_DRAFT, p["nz_main"]))
    radii = [OC4_MAIN_RADIUS * math.cos(0.5 * math.pi * i / p["nr_main"]) for i in range(p["nr_main"] + 1)]
    radii[-1] = 0.0
    _disc(out, 0.0, 0.0, radii, ring_counts(radii, p["n_main"]), -OC4_DRAFT, up=True)
    track(start, lambda c: (0.0, 0.0, -0.5 * OC4_DRAFT))
    polygons = [[_ring_point(0.0, 0.0, OC4_MAIN_RADIUS, p["n_main"], k, 0.0)[:2] for k in range(p["n_main"])]]

    if p["n_base"] != 2 * p["n_upper"]:
        raise ValueError("the base columns need twice the nodes of the upper columns around")
    middle = -0.5 * (OC4_DRAFT + OC4_BASE_TOP)
    for cx, cy in OC4_OFFSET_AXES:
        # Upper column wall, annulus on top of the base column, base column wall, base bottom
        start = len(out.panels)
        _wall(out, cx, cy, OC4_UPPER_RADIUS, p["n_upper"], _cosine(0.0, -OC4_BASE_TOP, p["nz_upper"]))
        _wall(out, cx, cy, OC4_BASE_RADIUS, p["n_base"], _cosine(-OC4_BASE_TOP, -OC4_DRAFT, p["nz_base"]))
        track(start, lambda c, cx=cx, cy=cy: (cx, cy, c[2]))
        start = len(out.panels)
        radii = _cosine(OC4_BASE_RADIUS, OC4_UPPER_RADIUS, p["nr_annulus"])
        counts = ring_counts(radii, p["n_base"])
        counts[-1] = p["n_upper"]
        _disc(out, cx, cy, radii, counts, -OC4_BASE_TOP, up=False)
        radii = [OC4_BASE_RADIUS * math.cos(0.5 * math.pi * i / p["nr_base"]) for i in range(p["nr_base"] + 1)]
        radii[-1] = 0.0
        _disc(out, cx, cy, radii, ring_counts(radii, p["n_base"]), -OC4_DRAFT, up=True)
        track(start, lambda c: (c[0], c[1], middle))
        polygons.append([_ring_point(cx, cy, OC4_UPPER_RADIUS, p["n_upper"], k, 0.0)[:2] for k in range(p["n_upper"])])
    return out, polygons, inside


def check_normals_each(mesh, inside):
    """Stop unless the normal of every panel points from the panel towards its own inside point."""
    for number, panel in enumerate(mesh.panels, start=1):
        points = [mesh.nodes[i - 1] for i in panel]
        centre = tuple(sum(q[k] for q in points) / len(points) for k in range(3))
        if sum(n * d for n, d in zip(hams_normal(points), sub(inside[number - 1], centre))) <= 0.0:
            raise ValueError("panel %d has a normal pointing into the fluid" % number)


def open_edges(mesh, tol=1e-9):
    """Check that a mesh is closed and consistently oriented except at z = 0: every edge not in z = 0 is used by two panels in
    opposite directions, every edge in z = 0 by one. Returns the edges in z = 0 as (start node, end node), 1-based."""
    use = {}
    for panel in mesh.panels:
        for k in range(len(panel)):
            a, b = panel[k], panel[(k + 1) % len(panel)]
            use.setdefault((min(a, b), max(a, b)), []).append((a, b))
    edges = []
    for (a, b), directed in use.items():
        on_surface = abs(mesh.nodes[a - 1][2]) < tol and abs(mesh.nodes[b - 1][2]) < tol
        if on_surface:
            if len(directed) != 1:
                raise ValueError("an edge in z = 0 is used by %d panels" % len(directed))
            edges.append(directed[0])
        elif len(directed) != 2 or directed[0] == directed[1]:
            raise ValueError("the mesh is open or not consistently oriented at the edge of nodes %d and %d" % (a, b))
    return edges


def mesh_properties(mesh):
    """Volume, centre of buoyancy and the waterplane integrals of a hull with HAMS normals (into the body), open at z = 0 only,
    from the panels by the divergence theorem: V = -sum n_z int z dA, V zb = -1/2 sum n_z int z^2 dA, and so on, and int f dA
    over the waterplane = sum n_z int f dA over the hull for f = 1, x, y, x^2, y^2, xy. Each panel is cut into triangles
    from its first vertex; the mid-edge rule is exact for these quadratic integrands on a plane panel."""
    sums = dict(vol=0.0, vzb=0.0, vxb=0.0, vyb=0.0, area=0.0, sx=0.0, sy=0.0, sxx=0.0, syy=0.0, sxy=0.0)
    for panel in mesh.panels:
        points = [mesh.nodes[i - 1] for i in panel]
        nz = hams_normal(points)[2]
        for k in range(1, len(points) - 1):
            t = (points[0], points[k], points[k + 1])
            area = 0.5 * math.sqrt(sum(c * c for c in cross(sub(t[1], t[0]), sub(t[2], t[0]))))
            mids = [tuple(0.5 * (t[i][c] + t[(i + 1) % 3][c]) for c in range(3)) for i in range(3)]

            def integral(f):
                return area * sum(f(*m) for m in mids) / 3.0
            sums["vol"] -= nz * integral(lambda x, y, z: z)
            sums["vzb"] -= 0.5 * nz * integral(lambda x, y, z: z * z)
            sums["vxb"] -= nz * integral(lambda x, y, z: x * z)
            sums["vyb"] -= nz * integral(lambda x, y, z: y * z)
            sums["area"] += nz * area
            sums["sx"] += nz * integral(lambda x, y, z: x)
            sums["sy"] += nz * integral(lambda x, y, z: y)
            sums["sxx"] += nz * integral(lambda x, y, z: x * x)
            sums["syy"] += nz * integral(lambda x, y, z: y * y)
            sums["sxy"] += nz * integral(lambda x, y, z: x * y)
    out = dict(sums)
    out.update(xb=sums["vxb"] / sums["vol"], yb=sums["vyb"] / sums["vol"], zb=sums["vzb"] / sums["vol"])
    return out


def oc4_exact():
    """Volume, zb and waterplane integrals of the columns (no braces) as exact cylinders."""
    pi = math.pi
    main_area = pi * OC4_MAIN_RADIUS ** 2
    upper_area = pi * OC4_UPPER_RADIUS ** 2
    base_area = pi * OC4_BASE_RADIUS ** 2
    vol = main_area * OC4_DRAFT + 3.0 * (upper_area * OC4_BASE_TOP + base_area * (OC4_DRAFT - OC4_BASE_TOP))
    zmom = -main_area * OC4_DRAFT ** 2 / 2.0 + 3.0 * (-upper_area * OC4_BASE_TOP ** 2 / 2.0
                                                     - base_area * (OC4_DRAFT ** 2 - OC4_BASE_TOP ** 2) / 2.0)
    area = main_area + 3.0 * upper_area
    syy = pi * OC4_MAIN_RADIUS ** 4 / 4.0 + sum(pi * OC4_UPPER_RADIUS ** 4 / 4.0 + upper_area * y * y for _, y in OC4_OFFSET_AXES)
    sxx = pi * OC4_MAIN_RADIUS ** 4 / 4.0 + sum(pi * OC4_UPPER_RADIUS ** 4 / 4.0 + upper_area * x * x for x, _ in OC4_OFFSET_AXES)
    return dict(vol=vol, zb=zmom / vol, area=area, sxx=sxx, syy=syy)


def oc4_matrices(props, rho=RHO, g=G):
    """Mass, hydrostatic restoring and mooring matrices about the origin on the SWL, for the water density rho of the solver.
    The report's mass, inertias and mooring are for OC4_RHO; they are scaled by rho / OC4_RHO, which leaves the motions and the
    loads over rho g of the report's platform unchanged (every term of the equation of motion then scales with rho) [derived].
    Restoring: C33 = rho g A, C34 = rho g Sy, C35 = -rho g Sx, C44 = rho g (Syy + V zb) - M g zg, C45 = -rho g Sxy,
    C46 = -rho g V xb, C55 = rho g (Sxx + V zb) - M g zg, C56 = -rho g V yb (CM on the centreline)."""
    scale = rho / OC4_RHO
    mass = scale * OC4_MASS
    zg = OC4_ZG
    m = [[0.0] * 6 for _ in range(6)]
    for k in range(3):
        m[k][k] = mass
    m[0][4] = m[4][0] = mass * zg
    m[1][3] = m[3][1] = -mass * zg
    m[3][3] = m[4][4] = scale * OC4_I_ROLL + mass * zg * zg
    m[5][5] = scale * OC4_I_YAW
    rg = rho * g
    c = [[0.0] * 6 for _ in range(6)]
    c[2][2] = rg * props["area"]
    c[2][3] = c[3][2] = rg * props["sy"]
    c[2][4] = c[4][2] = -rg * props["sx"]
    c[3][3] = rg * (props["syy"] + props["vol"] * props["zb"]) - mass * g * zg
    c[3][4] = c[4][3] = -rg * props["sxy"]
    c[3][5] = -rg * props["vol"] * props["xb"]
    c[4][4] = rg * (props["sxx"] + props["vol"] * props["zb"]) - mass * g * zg
    c[4][5] = -rg * props["vol"] * props["yb"]
    mooring = [[scale * value for value in row] for row in OC4_MOORING]
    return dict(mass=m, restoring=c, mooring=mooring, zg=zg)


def oc4_mooring_signs():
    """Estimate of the surge-pitch coupling of the mooring, for its sign [derived]. Three catenary lines at 120 degrees, fairleads
    at z_f = -14 m and radius R_f = 40.868 m (Table 5-1), h = 186 m above the anchors. A surge x1 moves each fairlead
    along its line by -x1 e_ix; the line's horizontal tension H and vertical pull V change with the horizontal distance X
    to the anchor at the rates H' and V' = (h / L_s) H' (catenary, suspended length L_s). About the origin
    C51 = -dMy/dx1 = -|z_f| k + 1.5 R_f V' with k = C11 = 1.5 H' (the horizontal part: the restoring force acts below the
    origin; the vertical part: the pull grows on the line behind and drops on the two ahead). L_s = V / w from the vertical
    pretension 1.839e6 N / 3 (Eq. 5-14) and the wet weight w = 108.63 x 9.80665 N/m (Table 5-1). The two parts nearly cancel,
    leaving C51 < 0; a pitch rotation lifts and lowers the fairleads, which gives C15 the same value (a catenary line is
    conservative, so dH/dz_f = dV/dX), and sway and roll give C24 = C42 = -C51 by the symmetry of the layout. Leaving out
    the transverse stiffness H / X of each line (about 1.1e3 N/m, with H about 8.5e5 N and X about 797 m) makes H' about 2%
    high; with it C51 is about -7.7e4, so the sign is firmer than the near cancellation suggests. The extracted Eq. 5-15 has
    two minus signs, which fits two negative entries (C15, C51).
    Returns (C51 estimate, |z_f| k, the vertical part)."""
    k = OC4_MOORING[0][0]
    h_prime = k / 1.5
    suspended = (1.839e6 / 3.0) / (108.63 * G)
    vertical = 1.5 * 40.868 * (186.0 / suspended) * h_prime
    return -14.0 * k + vertical, 14.0 * k, vertical


def _wall_levels(depth, panel_size):
    """Depths of the wall rows of a control surface, as wall_levels in src/ControlSurfaceMesh.f90: a top row of half the panel
    size, rows growing by 1.25 to the panel size, scaled so that the last level is -depth."""
    t1 = 0.5 * panel_size
    rows = [t1]
    tk = min(panel_size, t1 * 1.25)
    cum = t1
    while cum + 0.5 * tk < depth:
        rows.append(tk)
        cum += tk
        tk = min(panel_size, tk * 1.25)
    if len(rows) == 1:
        rows = [depth]
    else:
        factor = (depth - t1) / (sum(rows) - t1)
        rows = [t1] + [r * factor for r in rows[1:]]
    levels = [0.0]
    for r in rows:
        levels.append(levels[-1] - r)
    levels[-1] = -depth
    return levels


def oc4_control_surface(margin=8.0, bottom_margin=8.0, panel_size=4.0):
    """C around the whole platform: a vertical cylinder about the origin, margin beyond the base columns, with a bottom
    bottom_margin below the draft, about panel_size panels, HAMS normals towards the body (inwards, and up on the bottom).
    Returns the mesh and its radius."""
    radius = max(math.hypot(x, y) for x, y in OC4_OFFSET_AXES) + OC4_BASE_RADIUS + margin
    depth = OC4_DRAFT + bottom_margin
    n = 16 * int(math.ceil(2.0 * math.pi * radius / panel_size / 16.0))   # halves four times down the bottom
    out = Mesh()
    _wall(out, 0.0, 0.0, radius, n, _wall_levels(depth, panel_size))
    rows = int(math.ceil(radius / panel_size))
    radii = [radius * (1.0 - i / rows) for i in range(rows + 1)]
    _disc(out, 0.0, 0.0, radii, ring_counts(radii, n), -depth, up=True)
    return out, radius, n


class _Triangulation:
    """Delaunay triangulation of points in the plane by Bowyer-Watson, in pure Python (no SciPy on the cluster nodes)."""

    def __init__(self, points, order=None):
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        size = 100.0 * max(max(xs) - min(xs), max(ys) - min(ys))
        cx, cy = 0.5 * (max(xs) + min(xs)), 0.5 * (max(ys) + min(ys))
        self.points = list(points) + [(cx - size, cy - size), (cx + size, cy - size), (cx, cy + size)]
        n = len(points)
        self.triangles = {}
        self._add((n, n + 1, n + 2))
        for i in (range(n) if order is None else order):
            self._insert(i)
        self.triangles = [t for t in self.triangles if max(t) < n]

    def _add(self, t):
        (ax, ay), (bx, by), (cx, cy) = (self.points[i] for i in t)
        d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
        uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
        self.triangles[t] = (ux, uy, (ax - ux) ** 2 + (ay - uy) ** 2)

    def _insert(self, i):
        px, py = self.points[i]
        bad = [t for t, (ux, uy, r2) in self.triangles.items() if (px - ux) ** 2 + (py - uy) ** 2 < r2]
        count = {}
        for t in bad:
            for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                key = (min(e), max(e))
                count[key] = count.get(key, 0) + 1
            del self.triangles[t]
        for t in bad:
            for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                if count[(min(a, b), max(a, b))] == 1:
                    self._add((a, b, i))


def _signed_area(a, b, c):
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))


def oc4_free_surface(polygons, c_radius, c_count, panel_size=4.0):
    """F, the free surface z = 0 between the column waterlines (the polygons of oc4_hull) and Gamma_c (the top ring of C, of
    c_count nodes at c_radius), as triangles with HAMS normals -z. Its boundary nodes are those of the hull waterlines and of
    Gamma_c. Points: the boundary nodes, staggered rings around each waterline and inside Gamma_c whose spacing grows by
    FS_GROWTH per ring to that of triangles of area panel_size^2, and a hexagonal lattice of that spacing elsewhere, all
    symmetric about y = 0; then a Delaunay triangulation, with each column centre added (and its triangles removed) so that the
    holes stay holes. Checked: the triangles tile the region exactly (area), and the boundary edges are the waterline and
    Gamma_c edges. Returns the mesh."""
    far = panel_size / math.sqrt(math.sqrt(3.0) / 4.0)   # side of an equilateral triangle of area panel_size^2
    columns = [(0.0, 0.0, OC4_MAIN_RADIUS, len(polygons[0]))]
    columns += [(x, y, OC4_UPPER_RADIUS, len(poly)) for (x, y), poly in zip(OC4_OFFSET_AXES, polygons[1:])]
    fixed = [p for poly in polygons for p in poly]
    centres = [(x, y) for x, y, _, _ in columns]
    boundary = [_ring_point(0.0, 0.0, c_radius, c_count, k, 0.0)[:2] for k in range(c_count)]
    accepted = []

    def clear(p, spacing, own):
        """p keeps clear of Gamma_c, of the columns other than own, and of the points of other rings."""
        if math.hypot(p[0], p[1]) > c_radius - 0.5 * spacing:
            return False
        for index, (x, y, r, _) in enumerate(columns):
            if math.hypot(p[0] - x, p[1] - y) < r + (0.3 if index == own else 0.6) * spacing:
                return False
        return all(math.hypot(p[0] - q[0], p[1] - q[1]) >= 0.7 * spacing for q, owner in accepted if owner != own)

    # Rings out from each waterline (owner = column index) and in from Gamma_c (owner = -1), one level at a time
    sources = [(x, y, r, 2.0 * math.pi * r / n, 1.0, index) for index, (x, y, r, n) in enumerate(columns)]
    sources.append((0.0, 0.0, c_radius, 2.0 * math.pi * c_radius / c_count, -1.0, -1))
    state = [dict(radius=r, spacing=a) for _, _, r, a, _, _ in sources]
    level = 0
    while any(s["spacing"] < far for s in state):
        level += 1
        candidates = []
        for (x, y, _, a0, direction, owner), s in zip(sources, state):
            if s["spacing"] >= far:
                continue
            spacing = min(far, a0 * FS_GROWTH ** level)
            s["radius"] += direction * 0.5 * math.sqrt(3.0) * 0.5 * (s["spacing"] + spacing)
            s["spacing"] = spacing
            count = max(6, int(round(2.0 * math.pi * s["radius"] / spacing)))
            offset = math.pi / count if level % 2 else 0.0
            ring = [_ring_point(x, y, s["radius"], count, k, 0.0, offset)[:2] for k in range(count)]
            candidates += [(p, spacing, owner) for p in ring]
        kept = [(p, owner) for p, spacing, owner in candidates if clear(p, spacing, owner)
                and all(math.hypot(p[0] - q[0], p[1] - q[1]) >= 0.7 * spacing for q, _, o in candidates if o != owner)]
        accepted += kept
    rows = int(c_radius / (0.5 * math.sqrt(3.0) * far)) + 1
    for row in range(-rows, rows + 1):
        y = row * 0.5 * math.sqrt(3.0) * far
        shift = 0.5 * far if abs(row) % 2 else 0.0
        for col in range(-rows - 1, rows + 2):
            p = (col * far + shift, y)
            if clear(p, far, None) and all(math.hypot(p[0] - q[0], p[1] - q[1]) >= 0.75 * far for q, _ in accepted):
                accepted.append((p, -2))
    free = [p for p, _ in accepted]
    # The free points get a tiny deterministic shift, so that no four points are on a circle (symmetric pairs about y = 0 are)
    jitter = [((math.sin(12.9898 * i + 78.233) * 43758.5453) % 1.0 - 0.5, (math.sin(39.346 * i + 11.135) * 24634.6345) % 1.0
               - 0.5) for i in range(len(free))]
    free = [(p[0] + 1e-6 * j[0], p[1] + 1e-6 * j[1]) for p, j in zip(free, jitter)]
    points = fixed + boundary + centres + free
    first_centre = len(fixed) + len(boundary)
    # The centres and the free points go in first, so that no point is inserted on the circle of three others of a waterline
    order = list(range(first_centre, len(points))) + list(range(first_centre))
    triangles = [t for t in _Triangulation(points, order).triangles
                 if not any(first_centre <= i < first_centre + len(centres) for i in t)]

    out = Mesh()
    for t in triangles:
        corners = [points[i] + (0.0,) for i in t]
        _add_oriented(out, corners, (0.0, 0.0, -1.0))
    # Checks: positive areas that sum to the area of the region, and boundary edges on the waterlines and Gamma_c only
    total = sum(abs(_signed_area(*[points[i] for i in t])) for t in triangles)
    region = polygon_moments(boundary)[0] - sum(polygon_moments(poly)[0] for poly in polygons)
    if abs(total - region) > 1e-9 * region:
        raise ValueError("F does not tile the free surface: %.9f against %.9f m^2" % (total, region))
    smallest = min(abs(_signed_area(*[points[i] for i in t])) for t in triangles)
    if smallest < 1e-6 * panel_size ** 2:
        raise ValueError("F has a degenerate triangle")
    use = {}
    for t in triangles:
        for k in range(3):
            key = (min(t[k], t[(k + 1) % 3]), max(t[k], t[(k + 1) % 3]))
            use[key] = use.get(key, 0) + 1
    if max(use.values()) > 2:
        raise ValueError("an edge of F is used by more than two triangles")
    expected = set()
    start = 0
    for poly in polygons + [boundary]:
        for k in range(len(poly)):
            a, b = start + k, start + (k + 1) % len(poly)
            expected.add((min(a, b), max(a, b)))
        start += len(poly)
    found = set(key for key, count in use.items() if count == 1)
    if found != expected:
        raise ValueError("the boundary of F is not the waterlines and Gamma_c: %d edges against %d" % (len(found), len(expected)))
    return out


def triangle_quality(mesh):
    """Smallest angle in degrees of the triangles of a mesh."""
    worst = 180.0
    for panel in mesh.panels:
        if len(panel) != 3:
            continue
        p = [mesh.nodes[i - 1] for i in panel]
        for k in range(3):
            u, v = sub(p[(k + 1) % 3], p[k]), sub(p[(k + 2) % 3], p[k])
            cosine = sum(a * b for a, b in zip(u, v)) / math.sqrt(sum(a * a for a in u) * sum(b * b for b in v))
            worst = min(worst, math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return worst


def check_control_clearance(hull, control):
    """As check_control_surface of src/ControlSurfaceMesh.f90: no point of C (panel centroids, size sqrt(area), and the
    midpoints of the Gamma_c segments, size their length) within one size (the larger of the two) of a hull panel centroid.
    Returns the smallest distance over that size, or 2 if it is at least 2."""
    def centroids(mesh):
        out = []
        for panel in mesh.panels:
            points = [mesh.nodes[i - 1] for i in panel]
            if len(points) == 4:
                area = 0.5 * math.sqrt(sum(c * c for c in cross(sub(points[2], points[0]), sub(points[3], points[1]))))
            else:
                area = 0.5 * math.sqrt(sum(c * c for c in cross(sub(points[1], points[0]), sub(points[2], points[0]))))
            out.append((tuple(sum(q[k] for q in points) / len(points) for k in range(3)), math.sqrt(area)))
        return out
    hull_points = centroids(hull)
    control_points = centroids(control)
    for a, b in open_edges(control):
        pa, pb = control.nodes[a - 1], control.nodes[b - 1]
        control_points.append((tuple(0.5 * (u + v) for u, v in zip(pa, pb)), math.sqrt(sum((u - v) ** 2 for u, v in zip(pa, pb)))))
    # Pairs farther apart than twice the largest size cannot fail and are not looked at (bins of that size)
    cell = 2.0 * max(size for _, size in hull_points + control_points)
    bins = {}
    for h, hsize in hull_points:
        bins.setdefault(tuple(int(math.floor(v / cell)) for v in h), []).append((h, hsize))
    worst = 2.0
    for c, size in control_points:
        key = tuple(int(math.floor(v / cell)) for v in c)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for h, hsize in bins.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        worst = min(worst, math.sqrt(sum((a - b) ** 2 for a, b in zip(c, h))) / max(size, hsize))
    if worst < 1.0:
        raise ValueError("a panel of C is within one panel size of the hull")
    return min(worst, 2.0)


def read_mesh(path):
    """A mesh file in the HullMesh.pnl format, as written by write_mesh."""
    lines = open(path).read().splitlines()
    npan, nnode = (int(v) for v in lines[3].split()[:2])
    mesh = Mesh()
    mesh.nodes = [tuple(float(v) for v in line.split()[1:4]) for line in lines[6:6 + nnode]]
    start = 6 + nnode + 3
    mesh.panels = [[int(v) for v in line.split()[2:]] for line in lines[start:start + npan]]
    return mesh


def check_written_surfaces(folder):
    """The files as HAMS reads them: F lies in z = 0 and its boundary edges are the hull waterline edges and the edges of
    Gamma_c, node for node (to 1e-8 m); C passes the normal test of read_file_surface (sum n.(x - c) A < 0 about the centre c
    of the waterline box) and has a bottom (a normal with n_z > 0.5, as b2_control_check asks)."""
    hull, control, free = (read_mesh(os.path.join(folder, name)) for name in ("HullMesh.pnl", "CSMesh_1.pnl", "FSMesh_1.pnl"))

    def key(mesh, a, b):
        pa, pb = (tuple(round(c, 8) for c in mesh.nodes[i - 1][:2]) for i in (a, b))
        return (min(pa, pb), max(pa, pb))
    if any(abs(z) > 0.0 for _, _, z in free.nodes):
        raise ValueError("FSMesh_1.pnl is not in z = 0")
    use = {}
    for panel in free.panels:
        for k in range(len(panel)):
            edge = key(free, panel[k], panel[(k + 1) % len(panel)])
            use[edge] = use.get(edge, 0) + 1
    free_boundary = set(edge for edge, count in use.items() if count == 1)
    expected = set(key(hull, a, b) for a, b in open_edges(hull)) | set(key(control, a, b) for a, b in open_edges(control))
    if free_boundary != expected:
        raise ValueError("FSMesh_1.pnl does not meet the hull waterlines and Gamma_c: %d boundary edges against %d"
                         % (len(free_boundary), len(expected)))
    xs = [hull.nodes[i - 1][0] for edge in open_edges(hull) for i in edge]
    ys = [hull.nodes[i - 1][1] for edge in open_edges(hull) for i in edge]
    centre = (0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)), 0.0)
    total, bottom = 0.0, False
    for panel in control.panels:
        points = [control.nodes[i - 1] for i in panel]
        n = hams_normal(points)
        c = tuple(sum(q[k] for q in points) / len(points) for k in range(3))
        area = 0.5 * math.sqrt(sum(v * v for v in (cross(sub(points[2], points[0]), sub(points[3], points[1])) if len(points) == 4
                                                    else cross(sub(points[1], points[0]), sub(points[2], points[0])))))
        total += sum(a * b for a, b in zip(n, sub(c, centre))) * area
        bottom = bottom or n[2] > 0.5
    if not total < 0.0 or not bottom:
        raise ValueError("CSMesh_1.pnl fails the normal test of HAMS or has no bottom")


def write_oc4(folder, multiple=1, external_damping=None, check_clearance=True):
    """Write HullMesh.pnl, Hydrostatic.in (mooring as the external restoring, origin on the SWL as the rotation centre),
    CSMesh_1.pnl and FSMesh_1.pnl of OC4 into folder, after the checks: the hull is closed below z = 0 and consistently
    oriented, its normals point into the body, C has its normals towards the body and a bottom, F tiles the free surface
    between the waterlines and Gamma_c. Returns a dict with the meshes, the mesh properties and the matrices."""
    hull, polygons, inside = oc4_hull(multiple)
    check_normals_each(hull, inside)
    waterline = open_edges(hull)
    if len(waterline) != sum(len(poly) for poly in polygons):
        raise ValueError("the hull waterline has %d edges, not those of the four columns" % len(waterline))
    props = mesh_properties(hull)
    matrices = oc4_matrices(props)
    control, c_radius, c_count = oc4_control_surface(**OC4_CS)
    c_edges = open_edges(control)
    if len(c_edges) != c_count:
        raise ValueError("C is not closed below z = 0")
    # C normals towards the body: the sign test of read_file_surface, and every panel separately
    for panel in control.panels:
        points = [control.nodes[i - 1] for i in panel]
        centre = tuple(sum(q[k] for q in points) / len(points) for k in range(3))
        n = hams_normal(points)
        if abs(n[2]) > 0.5:
            if n[2] < 0.0:
                raise ValueError("a bottom panel of C has its normal pointing down")
        elif n[0] * centre[0] + n[1] * centre[1] >= 0.0:
            raise ValueError("a wall panel of C has its normal pointing outwards")
    free = oc4_free_surface(polygons, c_radius, c_count, OC4_CS["panel_size"])
    for panel in free.panels:
        if hams_normal([free.nodes[i - 1] for i in panel])[2] >= 0.0:
            raise ValueError("a panel of F has its normal pointing up")
    clearance = check_control_clearance(hull, control) if check_clearance else None
    os.makedirs(folder, exist_ok=True)
    write_mesh(os.path.join(folder, "HullMesh.pnl"), hull)
    write_hydrostatic(os.path.join(folder, "Hydrostatic.in"),
                      dict(zg=matrices["zg"], mass=matrices["mass"], restoring=matrices["restoring"]), 0.0, 0.0,
                      matrices["mooring"], external_damping)
    write_mesh(os.path.join(folder, "CSMesh_1.pnl"), control)
    write_mesh(os.path.join(folder, "FSMesh_1.pnl"), free, title="Free Surface Mesh File")
    check_written_surfaces(folder)
    return dict(hull=hull, control=control, free=free, props=props, matrices=matrices, c_radius=c_radius,
                clearance=clearance)


def oc4_check(multiple=1, folder=None):
    """Write the OC4 files at a multiple (into a temporary folder by default) and print the checks: panel counts, the mesh
    volume and waterplane against the exact columns and against the report (at the report's rho), and the triangle quality
    of F."""
    import tempfile
    with tempfile.TemporaryDirectory() as scratch:
        out = write_oc4(folder or scratch, multiple)
    props, exact = out["props"], oc4_exact()
    rg = OC4_RHO * G
    c44 = rg * (props["syy"] + props["vol"] * props["zb"])
    c44_exact = rg * (exact["syy"] + exact["vol"] * exact["zb"])
    lines = ["OC4 m%g: hull %d panels, %d nodes; C %d panels (radius %.2f m); F %d triangles, smallest angle %.1f deg; "
             "C clearance %s panel sizes" % (multiple, len(out["hull"].panels), len(out["hull"].nodes),
                                               len(out["control"].panels), out["c_radius"], len(out["free"].panels),
                                               triangle_quality(out["free"]),
                                               "%.2f" % out["clearance"] if out["clearance"] < 2.0 else ">= 2")]
    for label, mesh_value, exact_value, report in (
            ("rho g V0 (N)", rg * props["vol"], rg * exact["vol"], OC4_REPORT["rho_g_v0"]),
            ("C33 (N/m)", rg * props["area"], rg * exact["area"], OC4_REPORT["c33"]),
            ("C44 w/o gravity (N m/rad)", c44, c44_exact, OC4_REPORT["c44"]),
            ("C55 w/o gravity (N m/rad)", rg * (props["sxx"] + props["vol"] * props["zb"]), c44_exact - rg * (
                exact["syy"] - exact["sxx"]), OC4_REPORT["c44"])):
        lines.append("  %-26s mesh %.5e  exact columns %.5e (%+.2f%%)  report %.5e (%+.2f%%)" % (
            label, mesh_value, exact_value, 100.0 * (mesh_value / exact_value - 1.0), report,
            100.0 * (mesh_value / report - 1.0)))
    lines.append("  V %.3f m^3 (exact columns %.3f), zb %.4f m (%.4f), xb %.2e, yb %.2e, waterplane Sx %.3f, Sy %.2e m^3" % (
        props["vol"], exact["vol"], props["zb"], exact["zb"], props["xb"], props["yb"], props["sx"], props["sy"]))
    return "\n".join(lines), out


def main():
    parser = argparse.ArgumentParser(description="Generate a HAMS hull mesh and Hydrostatic.in.")
    sub_parsers = parser.add_subparsers(dest="shape")
    sub_parsers.required = True

    cyl = sub_parsers.add_parser("cylinder")
    cyl.add_argument("--radius", type=float, required=True)
    cyl.add_argument("--ntheta", type=int, required=True)
    cyl.add_argument("--nr", type=int, required=True)
    cyl.add_argument("--graded", action="store_true", help="cosine spacing, finer at the waterline and the bilge")
    cyl.add_argument("--waterplane-out", help="also write the interior free surface mesh to this path, for IRSP = 1")

    bx = sub_parsers.add_parser("box")
    bx.add_argument("--lx", type=float, required=True)
    bx.add_argument("--ly", type=float, required=True)
    bx.add_argument("--nx", type=int, required=True)
    bx.add_argument("--ny", type=int, required=True)

    for p in (cyl, bx):
        p.add_argument("--draft", type=float, required=True)
        p.add_argument("--nz", type=int, required=True)
        p.add_argument("--x0", type=float, default=0.0)
        p.add_argument("--y0", type=float, default=0.0)
        p.add_argument("--half-y", action="store_true", help="keep y > 0 and set the y symmetry flag")
        p.add_argument("--out", required=True, help="path of the HullMesh file")
        p.add_argument("--format", choices=("hams", "nemoh"), default="hams", help="format of the mesh file")
        p.add_argument("--hydro-out", help="path of the hydrostatic file (default: Hydrostatic.in next to --out)")

    o = sub_parsers.add_parser("oc4")
    o.add_argument("--multiple", type=int, default=1, help="mesh multiple (1: about 1200 panels, 2: about 4600)")
    o.add_argument("--out-dir", help="folder of the files (default: a temporary folder)")

    args = parser.parse_args()
    if args.shape == "oc4":
        print(oc4_check(args.multiple, args.out_dir)[0])
        estimate, horizontal, vertical = oc4_mooring_signs()
        print("  mooring C15 = C51 estimate %.3g N (horizontal part %.3g, vertical part %.3g); OC4_MOORING has %.3g, %.3g"
              % (estimate, -horizontal, vertical, OC4_MOORING[0][4], OC4_MOORING[4][0]))
        return
    if args.shape == "cylinder":
        mesh, polygon = cylinder(args.radius, args.draft, args.ntheta, args.nz, args.nr, args.graded)
    else:
        mesh, polygon = box(args.lx, args.ly, args.draft, args.nx, args.ny, args.nz)
    hydro_path = args.hydro_out or os.path.join(os.path.dirname(args.out), "Hydrostatic.in")
    mesh, hydro = write_body(mesh, polygon, args.draft, args.x0, args.y0, args.half_y, args.out, hydro_path,
                             mesh_format=args.format)

    print("%d panels, %d nodes, normals checked" % (len(mesh.panels), len(mesh.nodes)))
    if args.shape == "cylinder" and args.waterplane_out:
        # Uniform even for a graded hull
        plane = write_waterplane(args.waterplane_out, args.radius, args.ntheta, args.nr, x0=args.x0, y0=args.y0)
        print("%d waterplane panels, %d nodes" % (len(plane.panels), len(plane.nodes)))
    print("V = %.6f, Awp = %.6f, zg = zb = %.6f" % (hydro["vol"], hydro["area"], hydro["zg"]))
    print("M = %.6E, I44 = %.6E, I55 = %.6E, I66 = %.6E" % (hydro["mass"][0][0], hydro["mass"][3][3],
                                                             hydro["mass"][4][4], hydro["mass"][5][5]))
    print("C33 = %.6E, C44 = %.6E, C55 = %.6E" % (hydro["restoring"][2][2], hydro["restoring"][3][3],
                                                   hydro["restoring"][4][4]))


if __name__ == "__main__":
    main()
