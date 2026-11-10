"""
Generate a HAMS hull mesh and a matching Hydrostatic.in for a surface-piercing vertical cylinder or box.

Usage:
    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --out HullMesh.pnl
    python make_mesh.py box --lx 20 --ly 10 --draft 5 --nx 10 --ny 5 --nz 4 --out HullMesh.pnl
    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --format nemoh --out hull.dat

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


def write_hydrostatic(path, hydro, x0, y0, external_restoring=None):
    """external_restoring maps a 0-based mode number to a value on the diagonal of the external restoring matrix."""
    zero = [[0.0] * 6 for _ in range(6)]
    external = [[0.0] * 6 for _ in range(6)]
    for mode, value in (external_restoring or {}).items():
        external[mode][mode] = value
    with open(path, "w") as f:
        f.write(" Center of Gravity:\n")
        f.write("  %.12E  %.12E  %.12E\n" % (x0, y0, hydro["zg"]))
        f.write(" Body Mass Matrix:\n")
        write_matrix(f, hydro["mass"])
        f.write(" External Linear Damping Matrix:\n")
        write_matrix(f, zero)
        f.write(" External Quadratic Damping Matrix:\n")
        write_matrix(f, zero)
        f.write(" Hydrostatic Restoring Matrix:\n")
        write_matrix(f, hydro["restoring"])
        f.write(" External Restoring Matrix:\n")
        write_matrix(f, external)


def write_body(mesh, polygon, draft, x0, y0, half_y, mesh_path, hydro_path, external_restoring=None, mesh_format="hams"):
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
    write_hydrostatic(hydro_path, hydro, x0, y0, external_restoring)
    return mesh, hydro


def write_cylinder(mesh_path, hydro_path, radius, draft, ntheta, nz, nr, graded=False, x0=0.0, y0=0.0, half_y=False,
                   external_restoring=None, edge_power=None):
    """Write the mesh and hydrostatics of a cylinder. Returns the written mesh."""
    mesh, polygon = cylinder(radius, draft, ntheta, nz, nr, graded, edge_power)
    return write_body(mesh, polygon, draft, x0, y0, half_y, mesh_path, hydro_path, external_restoring)[0]


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

    args = parser.parse_args()
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
