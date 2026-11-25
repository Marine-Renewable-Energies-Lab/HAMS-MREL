"""
Simple geometric meshes for HAMS-MREL: single bodies, arrays of bodies, and dipole (zero-thickness) sheets.

Hull shapes (each is a function returning a Body, and a subcommand of the command line):
    cylinder      vertical surface-piercing cylinder (wall and flat bottom)
    box           rectangular surface-piercing box
    hemisphere    lower half of a sphere, cut by the free surface
    sphere        sphere entirely below the free surface
    spar          surface of revolution of a polyline profile: a truncated cone, a stepped or tapered spar
    cone          shortcut of spar: a truncated cone from a top radius to a bottom radius
    heave-plate   cylinder column with a dipole sheet: the annulus of a plate at the column bottom
Dipole sheets, added to a hull body (add_sheets, or --disc and --plate on any subcommand, repeatable):
    disc          circular dipole sheet
    plate         rectangular dipole sheet
A sheet cannot be written alone: HAMS stops at input when HullMesh.pnl has no hull panels (NELEM = 0). The sheets add no
volume, so Hydrostatic.in comes from the hull, or from --mass, --cog and --gyration.

Fineness: set --panel-size (a target panel edge, from which the panel counts follow) or the counts of each shape
(--ntheta, --nz, --nr, ...). Counts override the panel size. Two gradings move the nodes towards the edges, between
0 (uniform) and 1 (sine spacing):
    --waterline-grading G   rows finer towards the waterline z = 0
    --edge-grading G        rows finer towards sharp edges (the bilge corner of a cylinder or box, the corners of a
                            spar profile) and towards the free rim of a dipole sheet
A segment fine at one end has the nodes (1 - G) u + G f(u), f(u) = 1 - cos(pi u / 2) (fine at its start) or sin(pi u / 2)
(fine at its end); fine at both ends, f(u) = (1 - cos(pi u)) / 2 with the larger G. The counts are those of uniform
spacing at the panel size, so grading makes the end panels smaller than the panel size.

Outputs, in --out-dir (default the current folder), for one body:
    HullMesh.pnl         hull panels and, for a sheet, the dipole panels after "# Start Definition of Dipole Panels",
                         with one shared node list and line 4 "NELEM NTND ISX ISY NDIPOLE" (NDIPOLE only when > 0)
    WaterPlaneMesh.pnl   with --lid: the interior free surface for irregular frequency removal, IRSP = 1
    Hydrostatic.in       mass matrix of a uniform solid filling the submerged volume (or --mass, --cog, --gyration),
                         hydrostatic restoring from the panels, zero external damping and restoring
and the ControlFile.in lines to paste (rotation centre, If_remove_irr_freq). With more than one --at X Y [ANGLE], the
same body is written at each position as the multi-body set HullMesh_i.pnl, WaterPlaneMesh_i.pnl, Hydrostatic_i.in
(meshes in the body's own frame, placed by the LCS lines; centre of gravity, rotation centre and matrices global), and the
Number_of_bodies, LCS_i and Reference_body_center_i lines are printed. Dipole sheets support one body only.

Every mesh is checked before it is written:
    - normals as HAMS computes them (CalTransNormals, src/NormalProcess.f90) point out of the fluid, into the body;
    - the hull is closed, apart from its waterline on z = 0, and every edge is shared by two panels in opposite order;
    - no node above z = 0, and no hull panel centroid on z = 0;
    - every dipole node at least one panel size (the longest panel edge, CalPanelChartLength) below z = 0, no dipole
      centroid on a hull centroid, and no dipole panel overlapping a coplanar hull panel, by the rules of DIPOLE_SETUP
      (src/ReadPanelMesh.f90);
    - lid nodes exactly on z = 0, and a warning if a lid centroid is closer to the hull wall than 0.3 waterline rows
      (docs/input_files.md, "Lid mesh requirements").

Examples:
    python hams_mesh.py cylinder --radius 5 --draft 10 --panel-size 0.8 --waterline-grading 0.5 --edge-grading 0.5 --lid
    python hams_mesh.py box --length 20 --width 10 --draft 5 --nx 20 --ny 10 --nz 6 --out-dir Input
    python hams_mesh.py spar --profile 3.25 0 3.25 -4 4.7 -12 4.7 -120 --panel-size 1.5 --edge-grading 0.5
    python hams_mesh.py heave-plate --column-radius 10 --draft 9.9 --plate-radius 15 --panel-size 0.75 --edge-grading 1
    python hams_mesh.py cylinder --radius 2 --draft 4 --panel-size 0.4 --at 0 0 --at 20 0 --at 10 15 30 --lid
    python hams_mesh.py cylinder --radius 3 --draft 6 --panel-size 0.5 --disc 2 12 --plate 4 2 10 8 0 --edge-grading 1

Constants as Const_mod in src/WavDynMods.f90: rho = 1000 kg/m^3, g = 9.80665 m/s^2.
"""

import argparse
import math
import os

import numpy as np

RHO = 1000.0
G = 9.80665
NODE_FORMAT = "%.10f"
DIPOLE_HEADER = "# Start Definition of Dipole Panels"
TOL = 1e-9


# ----------------------------------------------------------------------------------------------------------------------
# Mesh and body


class Mesh:
    """Nodes shared by the hull panels, the dipole panels and (separately) a lid. Panels are lists of 1-based node
    indices, 3 or 4 of them."""

    def __init__(self):
        self.nodes = []
        self.index = {}
        self.panels = []
        self.dipoles = []

    def node(self, x, y, z):
        key = (round(x, 9) + 0.0, round(y, 9) + 0.0, round(z, 9) + 0.0)
        if key not in self.index:
            self.nodes.append((float(x), float(y), float(z)))
            self.index[key] = len(self.nodes)
        return self.index[key]

    def add(self, points, dipole=False):
        ids = [self.node(*p) for p in points]
        if len(set(ids)) != len(ids):
            raise ValueError("degenerate panel with repeated nodes: %s" % (points,))
        (self.dipoles if dipole else self.panels).append(ids)

    def points(self, panel):
        return [self.nodes[i - 1] for i in panel]


class Body:
    """A mesh in the body's own frame and what is needed to cover its waterplane with a lid: waterline = ("circle", R) or
    ("rect", lx, ly), with the height of its waterline row, or None for a body that does not cut the free surface."""

    def __init__(self, mesh, name, waterline=None, waterline_row=None, ntheta=None, exact=None):
        self.mesh = mesh
        self.name = name
        self.waterline = waterline
        self.waterline_row = waterline_row
        self.ntheta = ntheta
        self.exact = exact or {}


# ----------------------------------------------------------------------------------------------------------------------
# Spacing


def spacing(n, fine_start=0.0, fine_end=0.0):
    """n + 1 values from 0 to 1, graded towards a fine end with strengths between 0 (uniform) and 1 (sine)."""
    if n < 1:
        raise ValueError("a segment needs at least one panel")
    for g in (fine_start, fine_end):
        if not 0.0 <= g <= 1.0:
            raise ValueError("a grading must be between 0 and 1")
    g = max(fine_start, fine_end)
    out = []
    for i in range(n + 1):
        u = i / n
        if fine_start > 0.0 and fine_end > 0.0:
            f = 0.5 * (1.0 - math.cos(math.pi * u))
        elif fine_start > 0.0:
            f = 1.0 - math.cos(0.5 * math.pi * u)
        elif fine_end > 0.0:
            f = math.sin(0.5 * math.pi * u)
        else:
            f = u
        out.append((1.0 - g) * u + g * f)
    out[0], out[-1] = 0.0, 1.0
    return out


def count(length, panel_size, given=None, minimum=1):
    """Panel count of a segment: given, or the length over the panel size, rounded up."""
    if given is not None:
        if given < minimum:
            raise ValueError("a panel count must be at least %d" % minimum)
        return int(given)
    if panel_size is None or panel_size <= 0.0:
        raise ValueError("give --panel-size or the panel counts")
    return max(minimum, int(math.ceil(length / panel_size - 1e-9)))


def count_around(radius, panel_size, given=None):
    """Panels around a circle: given (even), or the circumference over the panel size rounded up to a multiple of 4,
    at least 8, which keeps the mesh symmetric about both axes."""
    if given is not None:
        if given < 4 or given % 2:
            raise ValueError("ntheta must be even and at least 4")
        return int(given)
    n = count(2.0 * math.pi * radius, panel_size, minimum=8)
    return n + (-n) % 4


# ----------------------------------------------------------------------------------------------------------------------
# Shapes


def revolve(mesh, profile, ntheta, dipole=False):
    """Surface of revolution of the profile [(r, z), ...], from the waterline (or the top pole) down to the axis or the
    sheet's end, with ntheta panels around. Points on the axis give triangles."""
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta)] + [0.0]

    def point(p, j):
        return (p[0] * math.cos(angles[j]), p[0] * math.sin(angles[j]), p[1])

    for a, b in zip(profile[:-1], profile[1:]):
        if a[0] < TOL and b[0] < TOL:
            continue
        for j in range(ntheta):
            if b[0] < TOL:
                mesh.add([(0.0, 0.0, b[1]), point(a, j + 1), point(a, j)], dipole)
            elif a[0] < TOL:
                mesh.add([point(b, j), point(b, j + 1), (0.0, 0.0, a[1])], dipole)
            else:
                mesh.add([point(b, j), point(b, j + 1), point(a, j + 1), point(a, j)], dipole)


def polyline(vertices, panel_size, counts, waterline_grading, edge_grading):
    """Subdivide the profile polyline [(r, z), ...]. A segment is fine at its ends by waterline_grading at z = 0 and by
    edge_grading at a corner of the polyline; not at the axis. Returns the points and the length of the first row."""
    if counts is not None and len(counts) != len(vertices) - 1:
        raise ValueError("give one count per profile segment (%d)" % (len(vertices) - 1))

    def kind(k):
        r, z = vertices[k]
        if r < TOL:
            return 0.0
        if k == 0:
            return waterline_grading if abs(z) < TOL else 0.0
        if k == len(vertices) - 1:
            return edge_grading
        d1 = np.subtract(vertices[k], vertices[k - 1])
        d2 = np.subtract(vertices[k + 1], vertices[k])
        cosine = np.dot(d1, d2) / (np.linalg.norm(d1) * np.linalg.norm(d2))
        return edge_grading if cosine < math.cos(math.radians(1.0)) else 0.0

    points = [tuple(vertices[0])]
    first_row = None
    for k in range(len(vertices) - 1):
        a, b = np.array(vertices[k], float), np.array(vertices[k + 1], float)
        length = float(np.linalg.norm(b - a))
        if length < TOL:
            raise ValueError("repeated profile vertex %s" % (vertices[k],))
        n = count(length, panel_size, None if counts is None else counts[k])
        for u in spacing(n, kind(k), kind(k + 1))[1:]:
            points.append(tuple(a + (b - a) * u))
        if first_row is None:
            first_row = float(np.linalg.norm(np.array(points[1]) - np.array(points[0])))
    return points, first_row


def spar(profile, panel_size=None, ntheta=None, counts=None, waterline_grading=0.0, edge_grading=0.0):
    """Body of revolution of the profile [(r0, 0), (r1, z1), ...] from the waterline downwards; the bottom disc to the axis
    is added. counts: one per segment, the bottom disc included."""
    vertices = [tuple(map(float, p)) for p in profile]
    if abs(vertices[0][1]) > TOL or vertices[0][0] <= 0.0:
        raise ValueError("the profile must start on the waterline, at a radius > 0 and z = 0")
    if any(z > TOL for _, z in vertices) or any(r < -TOL for r, _ in vertices):
        raise ValueError("the profile must lie at r >= 0 and z <= 0")
    if vertices[-1][0] > TOL:
        vertices.append((0.0, vertices[-1][1]))
    nt = count_around(max(r for r, _ in vertices), panel_size, ntheta)
    points, first_row = polyline(vertices, panel_size, counts, waterline_grading, edge_grading)
    mesh = Mesh()
    revolve(mesh, points, nt)
    return Body(mesh, "spar", ("circle", vertices[0][0]), first_row, nt)


def cylinder(radius, draft, panel_size=None, ntheta=None, nz=None, nr=None, waterline_grading=0.0, edge_grading=0.0):
    counts = None if nz is None and nr is None else [count(draft, panel_size, nz), count(radius, panel_size, nr)]
    body = spar([(radius, 0.0), (radius, -draft)], panel_size, ntheta, counts, waterline_grading, edge_grading)
    body.name = "cylinder"
    n = body.ntheta
    area = 0.5 * n * radius ** 2 * math.sin(2.0 * math.pi / n)
    body.exact = {"volume": math.pi * radius ** 2 * draft, "waterplane": math.pi * radius ** 2,
                  "volume_polygon": area * draft, "waterplane_polygon": area}
    return body


def cone(top_radius, bottom_radius, draft, panel_size=None, ntheta=None, nz=None, nr=None, waterline_grading=0.0,
         edge_grading=0.0):
    counts = None if nz is None and nr is None else [
        count(math.hypot(top_radius - bottom_radius, draft), panel_size, nz), count(bottom_radius, panel_size, nr)]
    body = spar([(top_radius, 0.0), (bottom_radius, -draft)], panel_size, ntheta, counts, waterline_grading, edge_grading)
    body.name = "cone"
    n = body.ntheta
    k = 0.5 * n * math.sin(2.0 * math.pi / n)
    a1, a2 = top_radius ** 2, bottom_radius ** 2
    body.exact = {"volume": math.pi * draft / 3.0 * (a1 + math.sqrt(a1 * a2) + a2), "waterplane": math.pi * a1,
                  "volume_polygon": k * draft / 3.0 * (a1 + math.sqrt(a1 * a2) + a2), "waterplane_polygon": k * a1}
    return body


def box(length, width, draft, panel_size=None, nx=None, ny=None, nz=None, waterline_grading=0.0, edge_grading=0.0):
    """Box centred on the axis, length along x. The x and y rows are graded towards the vertical and bottom edges, the
    z rows towards the waterline and the bottom edges."""
    xs = [-length / 2 + length * u for u in spacing(count(length, panel_size, nx), edge_grading, edge_grading)]
    ys = [-width / 2 + width * u for u in spacing(count(width, panel_size, ny), edge_grading, edge_grading)]
    zs = [-draft * u for u in spacing(count(draft, panel_size, nz), waterline_grading, edge_grading)]
    perimeter = ([(x, ys[0]) for x in xs] + [(xs[-1], y) for y in ys[1:]] + [(x, ys[-1]) for x in xs[-2::-1]]
                 + [(xs[0], y) for y in ys[-2:0:-1]])
    mesh = Mesh()
    for k in range(len(perimeter)):
        a, b = perimeter[k], perimeter[(k + 1) % len(perimeter)]
        for upper, lower in zip(zs[:-1], zs[1:]):
            mesh.add([a + (lower,), b + (lower,), b + (upper,), a + (upper,)])
    for ix in range(len(xs) - 1):
        for iy in range(len(ys) - 1):
            mesh.add([(xs[ix], ys[iy], -draft), (xs[ix], ys[iy + 1], -draft), (xs[ix + 1], ys[iy + 1], -draft),
                      (xs[ix + 1], ys[iy], -draft)])
    body = Body(mesh, "box", ("rect", length, width), abs(zs[1]))
    body.exact = {"volume": length * width * draft, "waterplane": length * width,
                  "volume_polygon": length * width * draft, "waterplane_polygon": length * width}
    return body


def hemisphere(radius, panel_size=None, ntheta=None, nphi=None, waterline_grading=0.0):
    """Lower half of a sphere centred on the waterline, rows of polar angle graded towards the waterline."""
    n = count(0.5 * math.pi * radius, panel_size, nphi)
    profile = [(radius * math.cos(0.5 * math.pi * u), -radius * math.sin(0.5 * math.pi * u))
               for u in spacing(n, waterline_grading, 0.0)]
    profile[-1] = (0.0, -radius)
    nt = count_around(radius, panel_size, ntheta)
    mesh = Mesh()
    revolve(mesh, profile, nt)
    row = math.hypot(profile[1][0] - profile[0][0], profile[1][1] - profile[0][1])
    body = Body(mesh, "hemisphere", ("circle", radius), row, nt)
    body.exact = {"volume": 2.0 / 3.0 * math.pi * radius ** 3, "waterplane": math.pi * radius ** 2,
                  "waterplane_polygon": 0.5 * nt * radius ** 2 * math.sin(2.0 * math.pi / nt)}
    return body


def sphere(radius, centre_depth, panel_size=None, ntheta=None, nphi=None):
    """Sphere entirely below the free surface, centre at z = -centre_depth, uniform rows of polar angle."""
    if centre_depth <= radius:
        raise ValueError("the sphere must lie below z = 0: centre_depth > radius")
    n = count(math.pi * radius, panel_size, nphi, minimum=2)
    profile = [(radius * math.sin(math.pi * u), -centre_depth + radius * math.cos(math.pi * u)) for u in spacing(n)]
    profile[0], profile[-1] = (0.0, -centre_depth + radius), (0.0, -centre_depth - radius)
    nt = count_around(radius, panel_size, ntheta)
    mesh = Mesh()
    revolve(mesh, profile, nt)
    body = Body(mesh, "sphere", None, None, nt)
    body.exact = {"volume": 4.0 / 3.0 * math.pi * radius ** 3, "waterplane": 0.0}
    return body


def annulus_profile(r_in, r_out, z, n, rim_grading):
    """Rings from r_in to r_out at depth z, graded towards the free outer rim."""
    return [(r_in + (r_out - r_in) * u, z) for u in spacing(n, 0.0, rim_grading)]


def heave_plate(column_radius, draft, plate_radius, panel_size=None, ntheta=None, nz=None, nr=None, nplate=None,
                waterline_grading=0.0, edge_grading=0.0, rim_grading=None):
    """Cylinder column (wall and bottom disc, the hull) with a dipole sheet: the annulus column_radius <= r <= plate_radius
    at the column bottom, graded towards its free outer rim (rim_grading, default edge_grading), as the heave plate of
    docs/input_files.md. The annulus shares the nodes of the column's bottom ring."""
    if plate_radius <= column_radius:
        raise ValueError("the plate radius must be larger than the column radius")
    nt = count_around(plate_radius, panel_size, ntheta)
    body = cylinder(column_radius, draft, panel_size, nt, nz, nr, waterline_grading, edge_grading)
    n = count(plate_radius - column_radius, panel_size, nplate)
    revolve(body.mesh, annulus_profile(column_radius, plate_radius, -draft, n,
                                       edge_grading if rim_grading is None else rim_grading)[::-1], nt, dipole=True)
    body.name = "heave-plate"
    body.exact["sheet_area_polygon"] = 0.5 * nt * math.sin(2.0 * math.pi / nt) * (plate_radius ** 2 - column_radius ** 2)
    body.exact["sheet_area"] = math.pi * (plate_radius ** 2 - column_radius ** 2)
    return body


def disc(radius, depth, panel_size=None, ntheta=None, nrings=None, rim_grading=0.0, centre=(0.0, 0.0)):
    """Circular dipole sheet of the given radius at z = -depth, centred on (x, y) = centre, rings graded towards the rim,
    triangles at the centre. A sheet alone cannot be written: add it to a hull body with add_sheets."""
    nt = count_around(radius, panel_size, ntheta)
    n = count(radius, panel_size, nrings)
    mesh = Mesh()
    revolve(mesh, [(radius * u, -depth) for u in spacing(n, 0.0, rim_grading)][::-1], nt, dipole=True)
    mesh.nodes = [(x + centre[0], y + centre[1], z) for x, y, z in mesh.nodes]
    body = Body(mesh, "disc", None, None, nt)
    body.exact = {"sheet_area": math.pi * radius ** 2,
                  "sheet_area_polygon": 0.5 * nt * radius ** 2 * math.sin(2.0 * math.pi / nt)}
    return body


def plate(length, width, depth, panel_size=None, nx=None, ny=None, rim_grading=0.0, centre=(0.0, 0.0)):
    """Rectangular dipole sheet at z = -depth, centred on (x, y) = centre, graded towards its four edges. A sheet alone
    cannot be written: add it to a hull body with add_sheets."""
    xs = [centre[0] - length / 2 + length * u for u in spacing(count(length, panel_size, nx), rim_grading, rim_grading)]
    ys = [centre[1] - width / 2 + width * u for u in spacing(count(width, panel_size, ny), rim_grading, rim_grading)]
    mesh = Mesh()
    for ix in range(len(xs) - 1):
        for iy in range(len(ys) - 1):
            mesh.add([(xs[ix], ys[iy], -depth), (xs[ix + 1], ys[iy], -depth), (xs[ix + 1], ys[iy + 1], -depth),
                      (xs[ix], ys[iy + 1], -depth)], dipole=True)
    body = Body(mesh, "plate", None, None)
    body.exact = {"sheet_area": length * width, "sheet_area_polygon": length * width}
    return body


def add_sheets(body, *sheets):
    """Add the dipole panels of the sheets (disc, plate, or another body's dipoles) to the body, in its frame. Nodes that
    coincide with the body's nodes are shared. The sheets add no volume, so the hydrostatics stay those of the hull.
    Returns the body."""
    for sheet in sheets:
        for panel in sheet.mesh.dipoles:
            body.mesh.add(sheet.mesh.points(panel), dipole=True)
        body.name += " + " + sheet.name
    return body


def lid(body, panel_size=None):
    """Interior free surface of a body (WaterPlaneMesh.pnl, IRSP = 1), nodes on z = 0. The outer ring or frame is one
    waterline row wide, so its centroids sit about half a row from the wall (docs/input_files.md, "Lid mesh requirements");
    the inside is uniform at the panel size (or at the width of the outer ring). Normals point up."""
    if body.waterline is None:
        raise ValueError("the %s does not cut the free surface, so it has no lid" % body.name)
    w = body.waterline_row
    size = panel_size or w
    mesh = Mesh()
    if body.waterline[0] == "circle":
        radius = body.waterline[1]
        w = min(w, 0.5 * radius)
        n = count(radius - w, size)
        radii = [radius] + [(radius - w) * (1.0 - u) for u in spacing(n)]
        nt = body.ntheta
        angles = [2.0 * math.pi * j / nt for j in range(nt)] + [0.0]
        for a, b in zip(radii[:-1], radii[1:]):
            for j in range(nt):
                pa = [(a * math.cos(angles[k]), a * math.sin(angles[k]), 0.0) for k in (j, j + 1)]
                if b < TOL:
                    mesh.add([pa[0], pa[1], (0.0, 0.0, 0.0)])
                else:
                    mesh.add([pa[0], pa[1], (b * math.cos(angles[j + 1]), b * math.sin(angles[j + 1]), 0.0),
                              (b * math.cos(angles[j]), b * math.sin(angles[j]), 0.0)])
    else:
        lx, ly = body.waterline[1:]
        w = min(w, 0.25 * lx, 0.25 * ly)

        def nodes(length):
            n = count(length - 2.0 * w, size)
            return [-length / 2] + [-length / 2 + w + (length - 2.0 * w) * u for u in spacing(n)] + [length / 2]
        xs, ys = nodes(lx), nodes(ly)
        for ix in range(len(xs) - 1):
            for iy in range(len(ys) - 1):
                mesh.add([(xs[ix], ys[iy], 0.0), (xs[ix + 1], ys[iy], 0.0), (xs[ix + 1], ys[iy + 1], 0.0),
                          (xs[ix], ys[iy + 1], 0.0)])
    for panel in mesh.panels:
        if hams_normal(mesh.points(panel))[2] < 0.0:
            panel.reverse()
    return mesh


# ----------------------------------------------------------------------------------------------------------------------
# Geometry as HAMS computes it


def hams_normal(points):
    """Unit normal as CalTransNormals (src/NormalProcess.f90): (x4 - x2) x (x3 - x1) for a quadrilateral,
    (x1 - x2) x (x3 - x2) for a triangle. For a hull panel it points out of the fluid."""
    p = np.array(points, float)
    un = np.cross(p[0] - p[1], p[2] - p[1]) if len(p) == 3 else np.cross(p[3] - p[1], p[2] - p[0])
    length = np.linalg.norm(un)
    if length < 1e-14:
        raise ValueError("degenerate panel: %s" % (points,))
    return un / length


def panel_size(points):
    """Longest edge, as CalPanelChartLength."""
    p = np.array(points, float)
    return max(np.linalg.norm(p[(k + 1) % len(p)] - p[k]) for k in range(len(p)))


def area_centroid(points):
    """Area-weighted centroid of the panel split into triangles (1, 2, 3) and (1, 3, 4)."""
    p = np.array(points, float)
    total, moment = 0.0, np.zeros(3)
    for tri in ([0, 1, 2], [0, 2, 3])[:len(p) - 2]:
        a = 0.5 * np.linalg.norm(np.cross(p[tri[1]] - p[tri[0]], p[tri[2]] - p[tri[0]]))
        total += a
        moment += a * p[tri].mean(axis=0)
    return moment / total


def triangles(mesh, panels):
    """Triangles of the panels, as an array (n, 3, 3), ordered so that their right-hand normal points into the fluid
    (the opposite of the HAMS normal)."""
    out = []
    for panel in panels:
        p = mesh.points(panel)
        out.append([p[0], p[1], p[2]])
        if len(p) == 4:
            out.append([p[0], p[2], p[3]])
    return np.array(out, float).reshape(-1, 3, 3)


# Degree-5 rule on the triangle, 7 points (barycentric, weights summing to 1)
_A1, _B1, _A2, _B2 = 0.059715871789770, 0.470142064105115, 0.797426985353087, 0.101286507323456
_BARY = np.array([[1 / 3, 1 / 3, 1 / 3], [_A1, _B1, _B1], [_B1, _A1, _B1], [_B1, _B1, _A1],
                  [_A2, _B2, _B2], [_B2, _A2, _B2], [_B2, _B2, _A2]])
_W = np.array([0.225] + [0.132394152788506] * 3 + [0.125939180544827] * 3)


def surface_integral(tris, f, component):
    """Sum over the triangles of the integral of f(x, y, z) n_k dS, n out of the body (into the fluid)."""
    pts = np.einsum("qv,tvk->tqk", _BARY, tris)
    s = 0.5 * np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    values = f(pts[..., 0], pts[..., 1], pts[..., 2]) @ _W
    return float(np.sum(values * s[:, component]))


def volume_moment(tris, a, b, c):
    """Integral of x^a y^b z^c over the volume closed by the hull and the waterplane z = 0 (which adds nothing)."""
    if c >= 1:
        return surface_integral(tris, lambda x, y, z: x ** a * y ** b * z ** (c + 1) / (c + 1), 2)
    return surface_integral(tris, lambda x, y, z: x ** (a + 1) / (a + 1) * y ** b, 0)


def waterplane_moment(tris, a, b):
    """Integral of x^a y^b over the waterplane that closes the hull, from the hull alone (divergence theorem)."""
    return -surface_integral(tris, lambda x, y, z: x ** a * y ** b + 0.0 * z, 2)


# ----------------------------------------------------------------------------------------------------------------------
# Checks


def edges(panels):
    use = {}
    for number, panel in enumerate(panels, start=1):
        for k in range(len(panel)):
            a, b = panel[k], panel[(k + 1) % len(panel)]
            use.setdefault((min(a, b), max(a, b)), []).append((a, b, number))
    return use


def check_hull(mesh):
    """Raise ValueError unless the hull is closed (open only on z = 0), consistently ordered, with HAMS normals into the
    body, and below the free surface. Returns the boundary edges on z = 0."""
    if not mesh.panels:
        return []
    for number, panel in enumerate(mesh.panels, start=1):
        p = mesh.points(panel)
        if max(q[2] for q in p) > TOL:
            raise ValueError("hull panel %d has a node above z = 0" % number)
        if abs(area_centroid(p)[2]) < 1e-10:
            raise ValueError("hull panel %d has its centroid on z = 0, where HAMS stops (CalNormals)" % number)
        tri_normal = np.cross(np.subtract(p[1], p[0]), np.subtract(p[2], p[0]))
        if len(p) == 4:
            tri_normal = tri_normal + np.cross(np.subtract(p[2], p[0]), np.subtract(p[3], p[0]))
        if np.dot(hams_normal(p), tri_normal) >= 0.0:
            raise ValueError("hull panel %d: the HAMS normal does not oppose the node order (warped panel)" % number)
    boundary = []
    for (i, j), uses in edges(mesh.panels).items():
        if len(uses) > 2:
            raise ValueError("edge %d-%d is shared by %d hull panels" % (i, j, len(uses)))
        if len(uses) == 2 and uses[0][:2] == uses[1][:2]:
            raise ValueError("hull panels %d and %d run edge %d-%d in the same direction: normals flip there"
                             % (uses[0][2], uses[1][2], i, j))
        if len(uses) == 1:
            if abs(mesh.nodes[i - 1][2]) > TOL or abs(mesh.nodes[j - 1][2]) > TOL:
                raise ValueError("the hull is not closed: edge %d-%d of panel %d is free below z = 0"
                                 % (i, j, uses[0][2]))
            boundary.append(uses[0][:2])
    tris = triangles(mesh, mesh.panels)
    volume = volume_moment(tris, 0, 0, 0)
    if volume <= 0.0:
        raise ValueError("the hull normals point into the fluid (volume %.6g)" % volume)
    return boundary


def inside_panel(x, points, n, size):
    """INSIDE_PANEL of src/ReadPanelMesh.f90: x in the plane of the panel lies inside its outline or on it."""
    pos = neg = False
    p = np.array(points, float)
    for k in range(len(p)):
        e = p[(k + 1) % len(p)] - p[k]
        s = np.dot(np.cross(e, np.asarray(x) - p[k]), n) / max(np.linalg.norm(e), 1e-300)
        pos = pos or s > 1e-9 * size
        neg = neg or s < -1e-9 * size
    return not (pos and neg)


def check_dipoles(mesh):
    """Raise ValueError on the input checks of DIPOLE_SETUP: nodes below z = 0 (DIPOLE_CHECK_GEOMETRY) and by at least one
    panel size (DIPOLE_CHECK_CLEARANCE), no dipole centroid within 1e-6 panel sizes of a hull centroid, and no overlap
    with a coplanar hull panel (DIPOLE_CHECK_COPLANAR)."""
    if not mesh.dipoles:
        return
    hull = [mesh.points(p) for p in mesh.panels]
    hc = np.array([area_centroid(p) for p in hull]).reshape(-1, 3)
    hn = np.array([hams_normal(p) for p in hull]).reshape(-1, 3)
    hs = np.array([panel_size(p) for p in hull])
    for number, panel in enumerate(mesh.dipoles, start=1):
        p = mesh.points(panel)
        size = panel_size(p)
        if max(q[2] for q in p) >= 0.0:
            raise ValueError("dipole panel %d has a node not below z = 0" % number)
        if max(q[2] for q in p) > -size:
            raise ValueError("dipole panel %d has a node less than one panel size (%.4g) below z = 0" % (number, size))
        if not hull:
            continue
        c, n = area_centroid(p), hams_normal(p)
        near = np.nonzero(np.linalg.norm(hc - c, axis=1) < 1e-6 * size)[0]
        if near.size:
            raise ValueError("dipole panel %d coincides with hull panel %d" % (number, near[0] + 1))
        s = np.maximum(size, hs)
        # A centroid can lie inside the other panel only within that panel's diameter (at most twice its longest edge)
        # of the panel's own centroid: a filter that changes nothing in the result of DIPOLE_CHECK_COPLANAR
        candidates = np.nonzero((np.linalg.norm(np.cross(n, hn), axis=1) < 1e-6)
                                & (np.abs((hc - c) @ n) < 1e-6 * s)
                                & (np.linalg.norm(hc - c, axis=1) <= 2.0 * np.maximum(size, hs) * (1 + 1e-9)))[0]
        for j in candidates:
            if inside_panel(hc[j], p, n, s[j]) or inside_panel(c, hull[j], n, s[j]):
                raise ValueError("dipole panel %d overlaps hull panel %d in the same plane" % (number, j + 1))
    for (i, j), uses in edges(mesh.dipoles).items():
        if len(uses) > 2:
            raise ValueError("edge %d-%d is shared by %d dipole panels" % (i, j, len(uses)))


def check_lid(lid_mesh, hull, boundary, row):
    """Raise unless every lid node is on z = 0; return the warnings on lid centroids close to the hull waterline."""
    if any(abs(z) > 1e-10 for _, _, z in lid_mesh.nodes):
        raise ValueError("a lid node is not on z = 0")
    segments = [(np.array(hull.nodes[a - 1][:2]), np.array(hull.nodes[b - 1][:2])) for a, b in boundary]
    closest = np.inf
    for panel in lid_mesh.panels:
        c = area_centroid(lid_mesh.points(panel))[:2]
        for a, b in segments:
            t = np.clip(np.dot(c - a, b - a) / np.dot(b - a, b - a), 0.0, 1.0)
            closest = min(closest, np.linalg.norm(c - a - t * (b - a)))
    if closest < 0.3 * row:
        return ["a lid centroid is %.3g from the hull wall, less than 0.3 of the waterline row (%.3g)" % (closest, row)]
    return []


# ----------------------------------------------------------------------------------------------------------------------
# Hydrostatics


def _comb(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def transform(nodes, x0, y0, angle):
    """Local to global, as ReadBodyMeshMulti: rotation by angle (rad) about z, then the shift (x0, y0)."""
    c, s = math.cos(angle), math.sin(angle)
    return [(x0 + c * x - s * y, y0 + s * x + c * y, z) for x, y, z in nodes]


def placed(mesh, x0, y0, angle):
    out = Mesh()
    out.nodes = transform(mesh.nodes, x0, y0, angle)
    out.panels, out.dipoles = mesh.panels, mesh.dipoles
    return out


def hydrostatics(mesh, centre, mass=None, cog=None, gyration=None):
    """Volume, buoyancy centre, waterplane, mass matrix and hydrostatic restoring about centre = (xr, yr, zr), in the
    frame of the mesh. The mass is that of a uniform solid filling the submerged volume (rho V, centre of gravity at the
    centre of buoyancy, its inertia from the panels), unless mass, cog and gyration (radii about the centre of gravity,
    along the axes of the mesh) are given. The restoring is that of WAMIT (Lee and Newman, WAMIT User Manual, Sec. 3.1),
    from the waterplane moments and the volume, with the gravity term of the centre of gravity."""
    xr, yr, zr = centre
    # The moments are taken in the mesh frame, where the waterplane is z = 0, then moved to the centre
    tris = triangles(mesh, mesh.panels) if mesh.panels else None
    raw_v, raw_w = {}, {}

    def vmc(a, b, c):
        """Volume integral of (x - xr)^a (y - yr)^b (z - zr)^c."""
        total = 0.0
        for i in range(a + 1):
            for j in range(b + 1):
                for k in range(c + 1):
                    if (i, j, k) not in raw_v:
                        raw_v[i, j, k] = volume_moment(tris, i, j, k)
                    total += (_comb(a, i) * _comb(b, j) * _comb(c, k) * (-xr) ** (a - i) * (-yr) ** (b - j)
                              * (-zr) ** (c - k) * raw_v[i, j, k])
        return total

    def wmc(a, b):
        """Waterplane integral of (x - xr)^a (y - yr)^b."""
        total = 0.0
        for i in range(a + 1):
            for j in range(b + 1):
                if (i, j) not in raw_w:
                    raw_w[i, j] = waterplane_moment(tris, i, j)
                total += _comb(a, i) * _comb(b, j) * (-xr) ** (a - i) * (-yr) ** (b - j) * raw_w[i, j]
        return total

    if tris is not None:
        vol = vmc(0, 0, 0)
        xb = np.array([vmc(1, 0, 0), vmc(0, 1, 0), vmc(0, 0, 1)]) / vol
        awp, sx, sy = wmc(0, 0), wmc(1, 0), wmc(0, 1)
        sxx, syy, sxy = wmc(2, 0), wmc(0, 2), wmc(1, 1)
    else:
        vol, xb, awp, sx, sy, sxx, syy, sxy = 0.0, np.zeros(3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    if mass is None:
        m = RHO * vol
        xg = xb.copy()
        if vol > 0.0:
            vm = lambda a, b, c: RHO * vmc(a, b, c)
            ixx, iyy, izz = vm(2, 0, 0), vm(0, 2, 0), vm(0, 0, 2)
            inertia = np.array([[iyy + izz, -vm(1, 1, 0), -vm(1, 0, 1)],
                                [-vm(1, 1, 0), ixx + izz, -vm(0, 1, 1)],
                                [-vm(1, 0, 1), -vm(0, 1, 1), ixx + iyy]])
        else:
            inertia = np.zeros((3, 3))
    else:
        m = float(mass)
        xg = np.array(cog, float) - np.array(centre, float)
        r = np.array(gyration if gyration is not None else (0.0, 0.0, 0.0), float)
        inertia = np.diag(m * r ** 2) + m * (np.dot(xg, xg) * np.eye(3) - np.outer(xg, xg))
    mm = np.zeros((6, 6))
    mm[:3, :3] = m * np.eye(3)
    skew = np.array([[0.0, -xg[2], xg[1]], [xg[2], 0.0, -xg[0]], [-xg[1], xg[0], 0.0]])
    mm[:3, 3:] = -m * skew
    mm[3:, :3] = m * skew
    mm[3:, 3:] = inertia
    c = np.zeros((6, 6))
    c[2, 2] = RHO * G * awp
    c[2, 3] = c[3, 2] = RHO * G * sy
    c[2, 4] = c[4, 2] = -RHO * G * sx
    c[3, 3] = RHO * G * (syy + vol * xb[2]) - m * G * xg[2]
    c[3, 4] = c[4, 3] = -RHO * G * sxy
    c[3, 5] = -RHO * G * vol * xb[0] + m * G * xg[0]
    c[4, 4] = RHO * G * (sxx + vol * xb[2]) - m * G * xg[2]
    c[4, 5] = -RHO * G * vol * xb[1] + m * G * xg[1]
    return {"volume": vol, "buoyancy_centre": xb + np.array(centre), "waterplane": awp, "mass": m,
            "cog": xg + np.array(centre), "mass_matrix": mm, "restoring": c}


# ----------------------------------------------------------------------------------------------------------------------
# Files


def write_pnl(path, nodes, panels, dipoles=(), title="Hull Mesh File"):
    """HullMesh.pnl / WaterPlaneMesh.pnl in the layout of examples/ and make_mesh.py; dipole panels, if any, after the
    hull panels and the line DIPOLE_HEADER, with NDIPOLE as a fifth number on line 4."""
    with open(path, "w") as f:
        f.write("    --------------%s---------------\n\n" % title)
        if dipoles:
            f.write("    # Number of Panels, Nodes, X-Symmetry, Y-Symmetry and Dipole Panels\n")
            f.write("    %8d  %8d  %8d  %8d  %8d\n\n" % (len(panels), len(nodes), 0, 0, len(dipoles)))
        else:
            f.write("    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n")
            f.write("    %8d  %8d  %8d  %8d\n\n" % (len(panels), len(nodes), 0, 0))
        f.write("    # Start Definition of Node Coordinates     ! node_number   x   y   z\n")
        for number, xyz in enumerate(nodes, start=1):
            f.write("    %6d  " % number + "  ".join((NODE_FORMAT % (c + 0.0)).rjust(18) for c in xyz) + "\n")
        f.write("    # End Definition of Node Coordinates\n\n")
        f.write("    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ...\n")
        for number, panel in enumerate(panels, start=1):
            f.write("    %6d  %3d  " % (number, len(panel)) + "  ".join("%6d" % i for i in panel) + "\n")
        if dipoles:
            f.write("    %s\n" % DIPOLE_HEADER)
            for number, panel in enumerate(dipoles, start=1):
                f.write("    %6d  %3d  " % (number, len(panel)) + "  ".join("%6d" % i for i in panel) + "\n")
        f.write("    # End Definition of Node Relations\n\n")
        f.write("    --------------End %s---------------\n" % title)


def read_pnl(path):
    """Read a .pnl file as HAMS does (list-directed reads): 3 lines, line 4 (NELEM NTND ISX ISY [NDIPOLE]), 2 lines,
    NTND nodes, 3 lines, NELEM panels, and with NDIPOLE > 0 one line and NDIPOLE panels. Returns nodes, panels, dipoles,
    (isx, isy)."""
    with open(path) as f:
        lines = f.read().splitlines()
    head = [int(v) for v in lines[3].split()]
    nelem, ntnd, isx, isy = head[:4]
    ndip = head[4] if len(head) > 4 else 0
    k = 6
    nodes = []
    for i in range(ntnd):
        v = lines[k + i].split()
        if int(v[0]) != i + 1:
            raise ValueError("%s: node %d is numbered %s" % (path, i + 1, v[0]))
        nodes.append(tuple(float(c) for c in v[1:4]))
    k += ntnd + 3

    def read_panels(start, n):
        out = []
        for i in range(n):
            v = [int(c) for c in lines[start + i].split()[:6]]
            out.append(v[2:2 + v[1]])
            if v[1] not in (3, 4) or any(not 1 <= j <= ntnd for j in out[-1]):
                raise ValueError("%s: bad panel line %d" % (path, start + i + 1))
        return out
    panels = read_panels(k, nelem)
    dipoles = []
    if ndip:
        if lines[k + nelem].strip() != DIPOLE_HEADER:
            raise ValueError("%s: expected '%s' after the hull panels" % (path, DIPOLE_HEADER))
        dipoles = read_panels(k + nelem + 1, ndip)
    return nodes, panels, dipoles, (isx, isy)


def write_matrix(f, matrix):
    for row in matrix:
        f.write("".join("  %12.5E" % (v + 0.0) for v in row) + "\n")   # 6(2x,E12.5), ReadHydroStaticMulti


def write_hydrostatic(path, hydro):
    zero = np.zeros((6, 6))
    with open(path, "w") as f:
        f.write(" Center of Gravity:\n")
        f.write("  %.12E  %.12E  %.12E\n" % tuple(v + 0.0 for v in hydro["cog"]))
        f.write(" Body Mass Matrix:\n")
        write_matrix(f, hydro["mass_matrix"])
        f.write(" External Linear Damping Matrix:\n")
        write_matrix(f, zero)
        f.write(" External Quadratic Damping Matrix:\n")
        write_matrix(f, zero)
        f.write(" Hydrostatic Restoring Matrix:\n")
        write_matrix(f, hydro["restoring"])
        f.write(" External Restoring Matrix:\n")
        write_matrix(f, zero)


def read_hydrostatic(path):
    """Read Hydrostatic.in with the formats of ReadHydroStaticMulti: CG free, matrices 6(2x,E12.5)."""
    with open(path) as f:
        lines = f.read().splitlines()
    cog = [float(v) for v in lines[1].split()]
    blocks = []
    for start in (3, 10, 17, 24, 31):
        blocks.append(np.array([[float(lines[start + i][2 + 14 * j:14 * (j + 1)]) for j in range(6)] for i in range(6)]))
    return cog, blocks


def control_lines(positions, centres, irsp):
    """The ControlFile.in lines of the bodies (formats of ReadControlFile, src/InputFiles.f90: LCS 26x,4f12.3 with the
    angle in radians, rotation centre 28x,3f12.3)."""
    out = []
    if len(positions) > 1:
        out.append("   # Start of definition of bodies for multi-body interaction")
        out.append("    Number_of_bodies          %d" % len(positions))
        for i, (x, y, a) in enumerate(positions, start=1):
            out.append(("    LCS_%d" % i).ljust(26) + "%12.3f%12.3f%12.3f%12.3f" % (x, y, 0.0, a))
        out.append("   # End of definition of bodies for multi-body interaction")
        for i, c in enumerate(centres, start=1):
            out.append(("    Reference_body_center_%d" % i).ljust(28) + "%12.3f%12.3f%12.3f" % tuple(c))
    else:
        out.append("    Reference_body_center".ljust(28) + "%12.3f%12.3f%12.3f" % tuple(centres[0]))
    out.append("    If_remove_irr_freq       %d" % (1 if irsp else 0))
    return out


def write_input(body, out_dir=".", positions=((0.0, 0.0, 0.0),), irsp=False, centre_z=0.0, mass=None, cog=None,
                gyration=None, lid_panel_size=None):
    """Check and write the input files of the body at each position (x, y, angle in rad, rounded as HAMS reads them),
    one body without index, more with _i. mass, cog (in the body's own frame) and gyration replace the uniform solid.
    Returns the hydrostatics of each body, the ControlFile.in lines and the warnings."""
    mesh = body.mesh
    positions = [(round(x, 3), round(y, 3), round(a, 3)) for x, y, a in positions]
    if not mesh.panels:
        raise ValueError("the %s has no hull panels: HAMS stops at input when HullMesh.pnl has NELEM = 0, so a dipole "
                         "sheet cannot be written alone. Add it to a hull body (add_sheets, or --disc / --plate on a "
                         "hull subcommand)" % body.name)
    if mesh.dipoles and len(positions) > 1:
        raise ValueError("dipole panels support one body only")
    boundary = check_hull(mesh)
    check_dipoles(mesh)
    warnings = []
    lid_mesh = None
    if irsp:
        lid_mesh = lid(body, lid_panel_size)
        warnings += check_lid(lid_mesh, mesh, boundary, body.waterline_row)
    os.makedirs(out_dir, exist_ok=True)
    multi = len(positions) > 1
    results, centres = [], []
    for i, (x, y, a) in enumerate(positions, start=1):
        suffix = "_%d" % i if multi else ""
        if not multi and (x, y, a) != (0.0, 0.0, 0.0):
            nodes = transform(mesh.nodes, x, y, a)
        else:
            nodes = mesh.nodes
        write_pnl(os.path.join(out_dir, "HullMesh%s.pnl" % suffix), nodes, mesh.panels, mesh.dipoles)
        if lid_mesh is not None:
            lnodes = transform(lid_mesh.nodes, x, y, a) if not multi else lid_mesh.nodes
            write_pnl(os.path.join(out_dir, "WaterPlaneMesh%s.pnl" % suffix), lnodes, lid_mesh.panels,
                      title="Waterplane Mesh File")
        centre = (x, y, centre_z)
        gcog = None if cog is None else transform([tuple(cog)], x, y, a)[0]
        hydro = hydrostatics(placed(mesh, x, y, a), centre, mass, gcog, gyration)
        if hydro["volume"] > 0.0 and body.waterline is not None and min(hydro["restoring"][3, 3],
                                                                         hydro["restoring"][4, 4]) <= 0.0:
            warnings.append("body %d: a metacentric height is not positive, the body is unstable" % i)
        if mass is not None and hydro["volume"] > 0.0 and abs(mass - RHO * hydro["volume"]) > 1e-3 * mass:
            warnings.append("body %d: the mass %.6g differs from rho V = %.6g: not in equilibrium"
                            % (i, mass, RHO * hydro["volume"]))
        write_hydrostatic(os.path.join(out_dir, "Hydrostatic%s.in" % suffix), hydro)
        results.append(hydro)
        centres.append(centre)
    return results, control_lines(positions, centres, irsp), warnings


# ----------------------------------------------------------------------------------------------------------------------
# Plot


def plot(bodies, path, titles=None):
    """PNG of the meshes: hull panels grey, dipole panels orange, lid panels blue, with a legend."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    n = len(bodies)
    cols = min(n, 3)
    rows = (n + cols - 1) // cols
    fig = plt.figure(figsize=(6.0 * cols, 5.2 * rows))
    for k, (body, lid_mesh) in enumerate(bodies):
        ax = fig.add_subplot(rows, cols, k + 1, projection="3d")
        mesh = body.mesh
        groups = [(mesh, mesh.panels, "0.75", "hull panels (HullMesh.pnl)"),
                  (mesh, mesh.dipoles, "tab:orange", "dipole panels (after %s)" % DIPOLE_HEADER)]
        if lid_mesh is not None:
            groups.append((lid_mesh, lid_mesh.panels, "tab:blue", "lid panels (WaterPlaneMesh.pnl)"))
        handles, allp = [], []
        for m, panels, colour, label in groups:
            if not panels:
                continue
            polys = [m.points(p) for p in panels]
            allp += [q for p in polys for q in p]
            ax.add_collection3d(Poly3DCollection(polys, facecolor=colour, edgecolor="k", linewidths=0.2, alpha=0.9))
            handles.append(Patch(facecolor=colour, edgecolor="k", label="%s: %d" % (label, len(panels))))
        a = np.array(allp)
        lo, hi = a.min(axis=0), a.max(axis=0)
        mid, half = (lo + hi) / 2, max(hi - lo) / 2
        ax.set_xlim(mid[0] - half, mid[0] + half)
        ax.set_ylim(mid[1] - half, mid[1] + half)
        ax.set_zlim(mid[2] - half, mid[2] + half)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.view_init(elev=-25 if body.mesh.dipoles else 20, azim=-60)
        ax.set_title(titles[k] if titles else body.name, fontsize=9)
        ax.legend(handles=handles, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------
# Command line


def build(args):
    """The hull of the subcommand, with the sheets of --disc and --plate added."""
    body = build_hull(args)
    sheets = []
    for values in args.disc or []:
        if len(values) not in (2, 4):
            raise SystemExit("--disc takes RADIUS DEPTH or RADIUS DEPTH X Y")
        sheets.append(disc(values[0], values[1], args.panel_size, rim_grading=args.edge_grading,
                           centre=tuple(values[2:4]) or (0.0, 0.0)))
    for values in args.plate or []:
        if len(values) not in (3, 5):
            raise SystemExit("--plate takes LENGTH WIDTH DEPTH or LENGTH WIDTH DEPTH X Y")
        sheets.append(plate(values[0], values[1], values[2], args.panel_size, rim_grading=args.edge_grading,
                            centre=tuple(values[3:5]) or (0.0, 0.0)))
    return add_sheets(body, *sheets)


def build_hull(args):
    common = dict(panel_size=args.panel_size)
    wl, edge = args.waterline_grading, args.edge_grading
    if args.shape == "cylinder":
        return cylinder(args.radius, args.draft, ntheta=args.ntheta, nz=args.nz, nr=args.nr, waterline_grading=wl,
                        edge_grading=edge, **common)
    if args.shape == "box":
        return box(args.length, args.width, args.draft, nx=args.nx, ny=args.ny, nz=args.nz, waterline_grading=wl,
                   edge_grading=edge, **common)
    if args.shape == "hemisphere":
        return hemisphere(args.radius, ntheta=args.ntheta, nphi=args.nphi, waterline_grading=wl, **common)
    if args.shape == "sphere":
        return sphere(args.radius, args.centre_depth, ntheta=args.ntheta, nphi=args.nphi, **common)
    if args.shape == "spar":
        if len(args.profile) % 2 or len(args.profile) < 4:
            raise ValueError("--profile takes pairs r z, at least two")
        prof = list(zip(args.profile[::2], args.profile[1::2]))
        return spar(prof, ntheta=args.ntheta, counts=args.counts, waterline_grading=wl, edge_grading=edge, **common)
    if args.shape == "cone":
        return cone(args.top_radius, args.bottom_radius, args.draft, ntheta=args.ntheta, nz=args.nz, nr=args.nr,
                    waterline_grading=wl, edge_grading=edge, **common)
    return heave_plate(args.column_radius, args.draft, args.plate_radius, ntheta=args.ntheta, nz=args.nz,
                       nr=args.nr, nplate=args.nplate, waterline_grading=wl, edge_grading=edge, **common)


def parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog="Run a subcommand with -h for its options. See the module docstring for details.")
    sub = p.add_subparsers(dest="shape")
    sub.required = True

    def shape(name, help_text, waterline=True):
        s = sub.add_parser(name, help=help_text, description=help_text)
        s.add_argument("--panel-size", type=float, help="target panel edge [m]; panel counts override it")
        s.add_argument("--edge-grading", type=float, default=0.0,
                       help="0 to 1: finer towards sharp edges and free dipole rims (default 0, uniform)")
        if waterline:
            s.add_argument("--waterline-grading", type=float, default=0.0,
                           help="0 to 1: finer towards the waterline (default 0, uniform)")
        else:
            s.set_defaults(waterline_grading=0.0)
        s.add_argument("--out-dir", default=".", help="folder of the input files (default .)")
        s.add_argument("--at", nargs="+", type=float, action="append", metavar="X Y [ANGLE]",
                       help="position of a body and its rotation about z [deg]; repeat for an array (multi-body files)")
        s.add_argument("--centre-z", type=float, default=0.0, help="z of the rotation centre (default 0)")
        s.add_argument("--lid", action="store_true", help="also write WaterPlaneMesh.pnl, for IRSP = 1")
        s.add_argument("--lid-panel-size", type=float, help="panel size inside the lid (default the panel size)")
        s.add_argument("--mass", type=float, help="body mass [kg], instead of a uniform solid filling the volume")
        s.add_argument("--cog", type=float, nargs=3, metavar=("X", "Y", "Z"), help="centre of gravity, body frame")
        s.add_argument("--gyration", type=float, nargs=3, metavar=("RX", "RY", "RZ"),
                       help="radii of gyration about the centre of gravity [m]")
        s.add_argument("--plot", metavar="PNG", help="write a picture of the mesh")
        s.add_argument("--disc", nargs="+", type=float, action="append", metavar="R DEPTH [X Y]",
                       help="add a circular dipole sheet of radius R at z = -DEPTH, centred on (X, Y) (default the "
                            "axis), at the panel size, graded to its rim by --edge-grading; repeat for more")
        s.add_argument("--plate", nargs="+", type=float, action="append", metavar="L W DEPTH [X Y]",
                       help="add a rectangular dipole sheet L x W at z = -DEPTH, centred on (X, Y); repeat for more")
        return s

    s = shape("cylinder", "vertical surface-piercing cylinder")
    s.add_argument("--radius", type=float, required=True)
    s.add_argument("--draft", type=float, required=True)
    s.add_argument("--ntheta", type=int, help="panels around (even)")
    s.add_argument("--nz", type=int, help="rows on the wall")
    s.add_argument("--nr", type=int, help="rows on the bottom")
    s = shape("box", "rectangular surface-piercing box, length along x")
    s.add_argument("--length", type=float, required=True)
    s.add_argument("--width", type=float, required=True)
    s.add_argument("--draft", type=float, required=True)
    for k in ("--nx", "--ny", "--nz"):
        s.add_argument(k, type=int)
    s = shape("hemisphere", "lower half of a sphere centred on the waterline")
    s.add_argument("--radius", type=float, required=True)
    s.add_argument("--ntheta", type=int)
    s.add_argument("--nphi", type=int, help="rows from the waterline to the bottom")
    s = shape("sphere", "sphere below the free surface", waterline=False)
    s.add_argument("--radius", type=float, required=True)
    s.add_argument("--centre-depth", type=float, required=True, help="depth of the centre (> radius)")
    s.add_argument("--ntheta", type=int)
    s.add_argument("--nphi", type=int, help="rows from pole to pole")
    s = shape("spar", "body of revolution of a polyline profile (truncated cone, stepped or tapered spar)")
    s.add_argument("--profile", type=float, nargs="+", required=True, metavar="R Z",
                   help="vertices r z from the waterline (z = 0) down; the bottom disc to the axis is added")
    s.add_argument("--ntheta", type=int)
    s.add_argument("--counts", type=int, nargs="+", help="rows per segment, the bottom disc included")
    s = shape("cone", "truncated cone from a top radius at the waterline to a bottom radius")
    s.add_argument("--top-radius", type=float, required=True)
    s.add_argument("--bottom-radius", type=float, required=True)
    s.add_argument("--draft", type=float, required=True)
    s.add_argument("--ntheta", type=int)
    s.add_argument("--nz", type=int, help="rows on the side")
    s.add_argument("--nr", type=int, help="rows on the bottom")
    s = shape("heave-plate", "cylinder column with a dipole annulus at its bottom")
    s.add_argument("--column-radius", type=float, required=True)
    s.add_argument("--draft", type=float, required=True)
    s.add_argument("--plate-radius", type=float, required=True)
    s.add_argument("--ntheta", type=int)
    s.add_argument("--nz", type=int)
    s.add_argument("--nr", type=int)
    s.add_argument("--nplate", type=int, help="rings on the annulus")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    body = build(args)
    positions = []
    for at in args.at or [[0.0, 0.0]]:
        if len(at) not in (2, 3):
            raise SystemExit("--at takes X Y or X Y ANGLE")
        positions.append((at[0], at[1], math.radians(at[2]) if len(at) == 3 else 0.0))
    if args.mass is not None and args.cog is None:
        raise SystemExit("--mass needs --cog")
    results, lines, warnings = write_input(body, args.out_dir, positions, args.lid, args.centre_z, args.mass, args.cog,
                                           args.gyration, args.lid_panel_size or args.panel_size)
    mesh = body.mesh
    print("%s: %d hull panels, %d dipole panels, %d nodes, written to %s" % (
        body.name, len(mesh.panels), len(mesh.dipoles), len(mesh.nodes), os.path.abspath(args.out_dir)))
    h = results[0]
    print("body 1: volume %.6g m^3, waterplane %.6g m^2, mass %.6g kg, CG (%.4g, %.4g, %.4g), C33 %.6g, C44 %.6g, C55 %.6g"
          % (h["volume"], h["waterplane"], h["mass"], h["cog"][0], h["cog"][1], h["cog"][2], h["restoring"][2, 2],
             h["restoring"][3, 3], h["restoring"][4, 4]))
    for w in warnings:
        print("Warning: " + w)
    if mesh.dipoles:
        print("Dipole panels need ISOLV = 1, no symmetry and no active SecondOrder.in (docs/input_files.md).")
    print("ControlFile.in lines:")
    print("\n".join(lines))
    if args.plot:
        plot([(body, lid(body, args.lid_panel_size or args.panel_size) if args.lid else None)], args.plot)


if __name__ == "__main__":
    main()
