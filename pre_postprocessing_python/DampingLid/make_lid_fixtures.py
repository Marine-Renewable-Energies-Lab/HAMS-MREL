"""Writes the damping lid fixtures under test/data/damping_lid.

Each case copies the hull, the hydrostatics and the control file of test/data/image_term/image_cos (cylinder r = 5, draft 5,
depth 50, 1.35 rad/s), sets Damping_lid, and writes DampingLid.pnl from annulus(). Run from anywhere: python -I make_lid_fixtures.py
"""
import math
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "test", "data", "image_term", "image_cos", "input")
SOURCE_TWO = os.path.join(ROOT, "test", "data", "image_term", "image_cos_two", "input")
TARGET = os.path.join(ROOT, "test", "data", "damping_lid")
EOL = "\r\n"


def annulus(r0, r1, nt, nr, tri=False, reverse=False, z=0.0, isx=0, isy=0, z_node=None, sectors=None, centre=(0.0, 0.0),
            dampers=None):
    """Panels of the annulus r0 <= r <= r1 about centre on z = z, nt around and nr radial.

    r0 = 0 gives a disc: the first row is nt triangles about a centre node. tri splits every quad in two. reverse flips the
    vertex order of every panel. z_node = (index, z) moves one node off the plane. sectors = [(a0, a1), ...] in degrees keeps
    only the panels whose middle angle lies in one of them; the other nodes stay in the file, unused. dampers as in lid_file.
    Returns the text of a .pnl file.
    """
    radii = [r0 + (r1 - r0) * i / nr for i in range(nr + 1)]
    nodes = []
    ring = []  # ring[i][j] is the node number of ring i at angle j, 1-based
    for i, r in enumerate(radii):
        if i == 0 and r0 == 0.0:
            nodes.append((centre[0], centre[1], z))
            ring.append([len(nodes)] * nt)
            continue
        row = []
        for j in range(nt):
            a = 2.0 * math.pi * j / nt
            nodes.append((centre[0] + r * math.cos(a), centre[1] + r * math.sin(a), z))
            row.append(len(nodes))
        ring.append(row)
    panels = []
    for i in range(nr):
        for j in range(nt):
            middle = 360.0 * (j + 0.5) / nt
            if sectors is not None and not any(a0 <= middle <= a1 for a0, a1 in sectors):
                continue
            k = (j + 1) % nt
            a, b, c, d = ring[i][j], ring[i][k], ring[i + 1][k], ring[i + 1][j]
            if a == b:
                quads = [[a, c, d]]
            elif tri:
                quads = [[a, b, c], [a, c, d]]
            else:
                quads = [[a, b, c, d]]
            panels.extend(quads)
    if reverse:
        panels = [p[::-1] for p in panels]
    if z_node is not None:
        index, zz = z_node
        x, y, _ = nodes[index - 1]
        nodes[index - 1] = (x, y, zz)
    return lid_file(nodes, panels, isx, isy, dampers)


def lid_file(nodes, panels, isx=0, isy=0, dampers=None):
    """The text of a DampingLid.pnl file with these nodes and panels (lists of 1-based node numbers).

    dampers is a list of (ip1, ip2, eps), or (ip1, ip2, "text") to write an unreadable value; the default is one damper over
    all panels with eps = 0.2 m/s. The string "none" omits the block. A damper count of 0 is written by an empty list."""
    lines = ["    --------------Damping Lid Mesh File---------------", "",
             "    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry",
             "%11d%10d%10d%10d" % (len(panels), len(nodes), isx, isy), "",
             "    # Start Definition of Node Coordinates     ! node_number   x   y   z"]
    for n, (x, y, zz) in enumerate(nodes, 1):
        lines.append("%10d%20.10f%20.10f%20.10f" % (n, x, y, zz))
    lines += ["    # End Definition of Node Coordinates", "",
              "    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ..."]
    for n, p in enumerate(panels, 1):
        lines.append("%10d%5d" % (n, len(p)) + "".join("%8d" % v for v in p))
    lines += ["    # End Definition of Node Relations", ""]
    if isinstance(dampers, str) and dampers == "none":
        pass
    else:
        if dampers is None:
            dampers = [(1, len(panels), 0.2)]
        lines += ["    # Start Definition of Damping   ! number of dampers, then: first panel  last panel  epsilon (m/s)",
                  "%10d" % len(dampers)]
        for ip1, ip2, eps in dampers:
            value = eps if isinstance(eps, str) else "%.10f" % eps
            lines.append("%10d%10d%16s" % (ip1, ip2, value if len(value) < 16 else " " + value))  # a blank before a long value
        lines += ["    # End Definition of Damping", ""]
    lines += ["    --------------End Damping Lid Mesh File---------------", ""]
    return EOL.join(lines)


def repeat_vertex(lid_text):
    """The lid text with the third vertex of panel 1 replaced by its second, a repeated node."""
    lines = lid_text.split(EOL)
    k = lines.index("    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ...") + 1
    f = [int(v) for v in lines[k].split()]
    f[4] = f[3]
    lines[k] = "%10d%5d" % (f[0], f[1]) + "".join("%8d" % v for v in f[2:])
    return EOL.join(lines)


def wall(a, h, nt, nz):
    """The vertical wall r = a from z = 0 to the sea bed z = -h, with no bottom: nt around, nz rows at
    z = -(h/2)(1 - cos(pi i/nz)), as the rows of image_cos. Node and vertex order as in image_cos, so the normals point the same
    way. Returns the text of HullMesh.pnl."""
    zs = [-(h / 2.0) * (1.0 - math.cos(math.pi * i / nz)) for i in range(nz + 1)]
    zs[0], zs[-1] = 0.0, -h
    nodes = []
    for j in range(nt):
        t = 2.0 * math.pi * j / nt
        for z in zs:
            nodes.append((a * math.cos(t), a * math.sin(t), z))

    def node(j, i):
        return (j % nt) * (nz + 1) + i + 1

    panels = [[node(j, i + 1), node(j + 1, i + 1), node(j + 1, i), node(j, i)] for j in range(nt) for i in range(nz)]
    lines = ["    --------------Hull Mesh File---------------", "",
             "    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry",
             "%12d%10d%10d%10d" % (len(panels), len(nodes), 0, 0), "",
             "    # Start Definition of Node Coordinates     ! node_number   x   y   z"]
    for n, (x, y, z) in enumerate(nodes, 1):
        lines.append("%10d%20.10f%20.10f%20.10f" % (n, x, y, z))
    lines += ["    # End Definition of Node Coordinates", "",
              "    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ..."]
    for n, p in enumerate(panels, 1):
        lines.append("%10d%5d" % (n, len(p)) + "".join("%8d" % v for v in p))
    lines += ["    # End Definition of Node Relations", "", "    --------------End Hull Mesh File---------------", ""]
    return EOL.join(lines)


def wall_hydrostatics(a, h, rho=1000.0, g=9.80665):
    """Hydrostatic.in of the wall: only the motions depend on it, not the excitation. Mass rho pi a**2 h (the water the column
    would displace), centre of gravity at the centre of buoyancy z = -h/2, inertia of a solid column about the origin, and the
    waterplane restoring of a free-floating column of that shape: C33 = rho g pi a**2, C44 = C55 = rho g pi a**4/4 (the
    rho g V z_B and -m g z_G terms cancel). rho and g of src/WavDynMods.f90."""
    m = rho * math.pi * a * a * h
    zg = -h / 2.0
    ixx = m * (3 * a * a + h * h) / 12.0 + m * zg * zg
    izz = m * a * a / 2.0
    mass = [[m, 0, 0, 0, m * zg, 0], [0, m, 0, -m * zg, 0, 0], [0, 0, m, 0, 0, 0], [0, -m * zg, 0, ixx, 0, 0],
            [m * zg, 0, 0, 0, ixx, 0], [0, 0, 0, 0, 0, izz]]
    c33 = rho * g * math.pi * a * a
    c44 = rho * g * math.pi * a ** 4 / 4.0
    restoring = [[0.0] * 6 for _ in range(6)]
    restoring[2][2], restoring[3][3], restoring[4][4] = c33, c44, c44
    zero = [[0.0] * 6 for _ in range(6)]

    def block(title, mat):
        return [" " + title] + ["".join("%14.5E" % v for v in row) for row in mat]

    lines = [" Center of Gravity:", "  %.12E  %.12E  %.12E" % (0.0, 0.0, zg)]
    lines += block("Body Mass Matrix:", mass) + block("External Linear Damping Matrix:", zero)
    lines += block("External Quadratic Damping Matrix:", zero) + block("Hydrostatic Restoring Matrix:", restoring)
    lines += block("External Restoring Matrix:", zero)
    return EOL.join(lines) + EOL


def rectangle(x0, x1, y0, y1, nx, ny, dampers=None):
    """The lid x0 <= x <= x1, y0 <= y <= y1 on z = 0, nx by ny quads. Returns the text of a .pnl file."""
    nodes = [(x0 + (x1 - x0) * i / nx, y0 + (y1 - y0) * j / ny, 0.0) for j in range(ny + 1) for i in range(nx + 1)]

    def node(i, j):
        return j * (nx + 1) + i + 1

    panels = [[node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)] for j in range(ny) for i in range(nx)]
    return lid_file(nodes, panels, dampers=dampers)


def box(length, width, draft, nx, ny, nz):
    """The wetted surface of the box |x| <= length/2, |y| <= width/2, -draft <= z <= 0, open at z = 0: nx, ny, nz panels along
    x, y, z. The vertex order of every panel makes (p3 - p1) x (p4 - p2) point out of the body, as in image_cos. Returns the
    text of HullMesh.pnl."""
    a, b = length / 2.0, width / 2.0
    xs = [-a + length * i / nx for i in range(nx + 1)]
    ys = [-b + width * j / ny for j in range(ny + 1)]
    zs = [-draft * k / nz for k in range(nz + 1)]
    nodes = []
    index = {}

    def node(p):
        key = tuple(round(c, 9) for c in p)
        if key not in index:
            nodes.append(p)
            index[key] = len(nodes)
        return index[key]

    quads = []  # (corners, outward normal)
    for i in range(nx):
        for j in range(ny):
            quads.append(([(xs[i], ys[j], -draft), (xs[i + 1], ys[j], -draft), (xs[i + 1], ys[j + 1], -draft),
                           (xs[i], ys[j + 1], -draft)], (0.0, 0.0, -1.0)))
    for k in range(nz):
        for i in range(nx):
            for y, s in ((-b, -1.0), (b, 1.0)):
                quads.append(([(xs[i], y, zs[k]), (xs[i + 1], y, zs[k]), (xs[i + 1], y, zs[k + 1]), (xs[i], y, zs[k + 1])],
                              (0.0, s, 0.0)))
        for j in range(ny):
            for x, s in ((-a, -1.0), (a, 1.0)):
                quads.append(([(x, ys[j], zs[k]), (x, ys[j + 1], zs[k]), (x, ys[j + 1], zs[k + 1]), (x, ys[j], zs[k + 1])],
                              (s, 0.0, 0.0)))
    panels = []
    for corners, outward in quads:
        d1 = [corners[2][c] - corners[0][c] for c in range(3)]
        d2 = [corners[3][c] - corners[1][c] for c in range(3)]
        n = (d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0])
        if sum(n[c] * outward[c] for c in range(3)) < 0.0:
            corners = corners[::-1]
        panels.append([node(p) for p in corners])
    lines = ["    --------------Hull Mesh File---------------", "",
             "    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry",
             "%12d%10d%10d%10d" % (len(panels), len(nodes), 0, 0), "",
             "    # Start Definition of Node Coordinates     ! node_number   x   y   z"]
    for n, (x, y, z) in enumerate(nodes, 1):
        lines.append("%10d%20.10f%20.10f%20.10f" % (n, x, y, z))
    lines += ["    # End Definition of Node Coordinates", "",
              "    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ..."]
    for n, p in enumerate(panels, 1):
        lines.append("%10d%5d" % (n, len(p)) + "".join("%8d" % v for v in p))
    lines += ["    # End Definition of Node Relations", "", "    --------------End Hull Mesh File---------------", ""]
    return EOL.join(lines)


def box_hydrostatics(length, width, draft, rho=1000.0, g=9.80665):
    """Hydrostatic.in of the floating box, as wall_hydrostatics: mass rho L B T, centre of gravity at the centre of buoyancy
    z = -T/2, inertia of a solid box about the origin, C33 = rho g L B, C44 = rho g L B**3/12, C55 = rho g B L**3/12 (the
    rho g V z_B and -m g z_G terms cancel)."""
    m = rho * length * width * draft
    zg = -draft / 2.0
    ixx = m * (width ** 2 + draft ** 2) / 12.0 + m * zg * zg
    iyy = m * (length ** 2 + draft ** 2) / 12.0 + m * zg * zg
    izz = m * (length ** 2 + width ** 2) / 12.0
    mass = [[m, 0, 0, 0, m * zg, 0], [0, m, 0, -m * zg, 0, 0], [0, 0, m, 0, 0, 0], [0, -m * zg, 0, ixx, 0, 0],
            [m * zg, 0, 0, 0, iyy, 0], [0, 0, 0, 0, 0, izz]]
    restoring = [[0.0] * 6 for _ in range(6)]
    restoring[2][2] = rho * g * length * width
    restoring[3][3] = rho * g * length * width ** 3 / 12.0
    restoring[4][4] = rho * g * width * length ** 3 / 12.0
    zero = [[0.0] * 6 for _ in range(6)]

    def block(title, mat):
        return [" " + title] + ["".join("%14.5E" % v for v in row) for row in mat]

    lines = [" Center of Gravity:", "  %.12E  %.12E  %.12E" % (0.0, 0.0, zg)]
    lines += block("Body Mass Matrix:", mass) + block("External Linear Damping Matrix:", zero)
    lines += block("External Quadratic Damping Matrix:", zero) + block("Hydrostatic Restoring Matrix:", restoring)
    lines += block("External Restoring Matrix:", zero)
    return EOL.join(lines) + EOL


def control_text_two(points, irsp=0, origins=None, heading=None, lid_flag=1, frequency=None):
    """The image_cos_two control file (two bodies, the second 5 km away along y) with Damping_lid lid_flag inserted after 'Number of
    threads' and the field points replaced. origins = [(x, y), (x, y)] moves the local frames and the reference centres of the
    two bodies, heading replaces the wave heading, frequency the one frequency (rad/s)."""
    with open(os.path.join(SOURCE_TWO, "ControlFile.in"), "rb") as f:
        text = f.read().decode("ascii")
    out = []
    skip = False
    heading_next = False
    for line in text.split(EOL):
        if skip:
            if "#End Definition of Pressure" in line:
                skip = False
            else:
                continue
        if frequency is not None and line.strip() == "1.3500":
            line = "  %.4f" % frequency
        if heading_next:
            line = "  %.4f" % heading
            heading_next = False
        if heading is not None and "Number_of_headings" in line:
            heading_next = True
        if origins is not None:
            for k, (x, y) in enumerate(origins, 1):
                if line.strip().startswith("LCS_%d " % k):
                    line = line[:26] + "%12.3f%12.3f%12.3f%12.3f" % (x, y, 0.0, 0.0)
                if line.strip().startswith("Reference_body_center_%d " % k):
                    line = line[:28] + "%12.3f%12.3f%12.3f" % (x, y, 0.0)
        if irsp and "If_remove_irr_freq" in line:
            line = "    If_remove_irr_freq       %d" % irsp
        if "Number_of_field_points" in line:
            out.append("    Number_of_field_points     %d" % len(points))
            out.extend("    %.12f  %.12f  %.12f" % p for p in points)
            skip = True
            continue
        out.append(line)
        if "Number of threads" in line:
            out.append("    Damping_lid              %d" % lid_flag)
    return EOL.join(out).replace("generated by pre_postprocessing_python/SecondOrder/make_drift_fixtures.py",
                                 "generated by pre_postprocessing_python/DampingLid/make_lid_fixtures.py")


def write_case_two(name, lid_text, points, hulls=None, water_plane=None, irsp=0, origins=None, heading=None, lid_flag=1,
                   frequency=None):
    """A two-body case with a lid (Damping_lid lid_flag). hulls = [(mesh text, hydrostatics text)] * 2 replaces the hulls of image_cos_two,
    water_plane is the path of a WaterPlaneMesh.pnl copied for both bodies, used with irsp = 1."""
    folder = os.path.join(TARGET, name, "input")
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    for k in (1, 2):
        for m, f in enumerate(("HullMesh_%d.pnl" % k, "Hydrostatic_%d.in" % k)):
            if hulls is not None:
                with open(os.path.join(folder, f), "wb") as out:
                    out.write(hulls[k - 1][m].encode("ascii"))
            else:
                shutil.copy(os.path.join(SOURCE_TWO, f), os.path.join(folder, f))
        if water_plane is not None:
            shutil.copy(water_plane, os.path.join(folder, "WaterPlaneMesh_%d.pnl" % k))
    with open(os.path.join(folder, "ControlFile.in"), "wb") as f:
        f.write(control_text_two(points, irsp, origins, heading, lid_flag, frequency).encode("ascii"))
    with open(os.path.join(folder, "DampingLid.pnl"), "wb") as f:
        f.write(lid_text.encode("ascii"))


def centroids(lid_text, count):
    """The vertex averages of the first count panels of a .pnl text, as HAMS computes the panel centroids."""
    lines = lid_text.split(EOL)
    nodes = {}
    start = lines.index("    # Start Definition of Node Coordinates     ! node_number   x   y   z") + 1
    k = start
    while not lines[k].startswith("    # End"):
        f = lines[k].split()
        nodes[int(f[0])] = tuple(float(v) for v in f[1:4])
        k += 1
    k += 3
    out = []
    while len(out) < count and not lines[k].startswith("    # End"):
        f = [int(v) for v in lines[k].split()]
        ids = f[2:2 + f[1]]
        out.append(tuple(sum(nodes[i][c] for i in ids) / len(ids) for c in range(3)))
        k += 1
    return out


def off_centroid_points(lid_text):
    """Field points on z = 0 of the annulus r 6 to 10, 40 x 4: a lid node, the middle of a radial edge, a point inside the first
    panel off its centroid, the 10 first centroids, then (0, 8, 0) and (20, 0, 0)."""
    a = math.radians(2.0)
    return ([(7.0, 0.0, 0.0), (7.5, 0.0, 0.0), (6.25 * math.cos(a), 6.25 * math.sin(a), 0.0)] + centroids(lid_text, 10)
            + [(0.0, 8.0, 0.0), (20.0, 0.0, 0.0)])


def control_text(lid, solver=None, depth=None, limits=None, points=None, irsp=0, wave_number=None, zero_frequency=False):
    """The image_cos control file with Damping_lid inserted after 'Number of threads', and Solver_type when given. depth
    replaces the water depth, limits the 0_inf_frequency_limits switch, points the field points, and wave_number the frequency
    by a wave number (Input_frequency_type 2). zero_frequency puts 0 before the frequency (two frequencies)."""
    with open(os.path.join(SOURCE, "ControlFile.in"), "rb") as f:
        text = f.read().decode("ascii")
    out = []
    skip = False
    replaced = 0
    for line in text.split(EOL):
        if skip:
            if "#End Definition of Pressure" in line:
                skip = False
            else:
                continue
        if depth is not None and line.strip().startswith("Waterdepth"):
            line = "   Waterdepth  %.1fD0" % depth
        if limits is not None and "0_inf_frequency_limits" in line:
            line = "    0_inf_frequency_limits      %d" % limits
        if wave_number is not None and "Input_frequency_type" in line:
            line = "    Input_frequency_type        2"
        if wave_number is not None and line.strip() == "1.3500":
            line = "  %.4f" % wave_number
            replaced += 1
        if zero_frequency and "Number_of_frequencies" in line:
            line = "    Number_of_frequencies      2"
        if zero_frequency and line.strip() == "1.3500":
            line = "  0.0000  1.3500"
        if irsp and "If_remove_irr_freq" in line:
            line = "    If_remove_irr_freq       %d" % irsp
        if points is not None and "Number_of_field_points" in line:
            out.append("    Number_of_field_points     %d" % len(points))
            out.extend("    %.12f  %.12f  %.12f" % p for p in points)
            skip = True
            continue
        out.append(line)
        if "Number of threads" in line:
            if solver is not None:
                out.append("    Solver_type             %d" % solver)
            out.append("    Damping_lid              %d" % lid)
    if wave_number is not None and replaced != 1:
        raise RuntimeError("control_text: the frequency line '1.3500' of image_cos was not found once")
    return EOL.join(out).replace("generated by pre_postprocessing_python/SecondOrder/make_drift_fixtures.py",
                                 "generated by pre_postprocessing_python/DampingLid/make_lid_fixtures.py")


def write_case(name, lid_text, lid_flag=1, solver=None, depth=None, limits=None, points=None, hull=None, water_plane=None,
               irsp=0, wave_number=None, hull_files=None, zero_frequency=False):
    """hull replaces the folder of the hull and hydrostatics (default image_cos), water_plane is the path of a WaterPlaneMesh.pnl
    to copy, used with irsp = 1. hull_files = (mesh text, hydrostatics text) writes them in place of the copies."""
    folder = os.path.join(TARGET, name, "input")
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    for k, f in enumerate(("HullMesh.pnl", "Hydrostatic.in")):
        if hull_files is not None:
            with open(os.path.join(folder, f), "wb") as out:
                out.write(hull_files[k].encode("ascii"))
        else:
            shutil.copy(os.path.join(hull or SOURCE, f), os.path.join(folder, f))
    if water_plane is not None:
        shutil.copy(water_plane, os.path.join(folder, "WaterPlaneMesh.pnl"))
    with open(os.path.join(folder, "ControlFile.in"), "wb") as f:
        f.write(control_text(lid_flag, solver, depth, limits, points, irsp, wave_number, zero_frequency).encode("ascii"))
    if lid_text is not None:
        with open(os.path.join(folder, "DampingLid.pnl"), "wb") as f:
            f.write(lid_text.encode("ascii"))


def write_solve_case(name, lid_text, far=(), **kwargs):
    """A lid case that runs to the end, with its field points at the centroids of the first 10 lid panels, then the points in
    far."""
    write_case(name, lid_text, points=centroids(lid_text, 10) + list(far), **kwargs)


def main():
    outside = dict(r0=6.0, r1=10.0, nt=20, nr=2)
    write_case("lid_stop_missing", None)
    write_case("lid_stop_offplane", annulus(**outside, z_node=(5, 1e-3)))
    write_case("lid_stop_symmetry", annulus(**outside, isy=1))
    write_case("lid_stop_gmres", annulus(**outside), solver=2)
    write_case("lid_stop_inside", annulus(0.0, 2.5, 12, 2))
    write_case("lid_stop_value", annulus(**outside), lid_flag=2)
    write_case("lid_stop_repeated", repeat_vertex(annulus(**outside)))
    write_case("lid_flag_off", annulus(**outside), lid_flag=0)

    # Stops of the damping block, on the 40-panel annulus; each message names the line or panel (README).
    write_case("lid_stop_damp_missing", annulus(**outside, dampers="none"))
    write_case("lid_stop_damp_count", annulus(**outside, dampers=[]))
    write_case("lid_stop_damp_range", annulus(**outside, dampers=[(1, 41, 0.2)]))
    write_case("lid_stop_damp_order", annulus(**outside, dampers=[(30, 10, 0.2)]))
    write_case("lid_stop_damp_overlap", annulus(**outside, dampers=[(1, 25, 0.2), (20, 40, 0.1)]))
    write_case("lid_stop_damp_gap", annulus(**outside, dampers=[(1, 20, 0.2), (25, 40, 0.1)]))
    write_case("lid_stop_damp_negative", annulus(**outside, dampers=[(1, 40, -0.1)]))
    write_case("lid_stop_damp_value", annulus(**outside, dampers=[(1, 40, "abc")]))
    write_case("lid_stop_zero_frequency", annulus(**outside), zero_frequency=True)

    # Single-body solves (test_integration_lid.f90). Annulus r = 6 to 10 unless said.
    write_solve_case("lid_ann_40", annulus(6.0, 10.0, 20, 2))
    # lid_ann_160 and its epsilon variants have field points at all 160 lid centroids, then two points on z = 0 off the lid and
    # one in the fluid, for the lid power and the comparisons with lid_nolid_160 (the same file with Damping_lid 0).
    ann160 = annulus(6.0, 10.0, 40, 4)
    points160 = centroids(ann160, 160) + [(12.0, 0.0, 0.0), (20.0, 0.0, 0.0), (12.0, 0.0, -3.0)]
    write_case("lid_ann_160", ann160, points=points160)
    write_case("lid_ann_160_eps0", annulus(6.0, 10.0, 40, 4, dampers=[(1, 160, 0.0)]), points=points160)
    write_case("lid_ann_160_rigid", annulus(6.0, 10.0, 40, 4, dampers=[(1, 160, 1.0e6)]), points=points160)
    write_case("lid_ann_160_two_dampers", annulus(6.0, 10.0, 40, 4, dampers=[(1, 80, 0.2), (81, 160, 0.2)]), points=points160)
    write_case("lid_ann_160_split", annulus(6.0, 10.0, 40, 4, dampers=[(1, 80, 0.2), (81, 160, 1.0)]), points=points160)
    write_case("lid_nolid_160", ann160, lid_flag=0, points=points160)
    write_solve_case("lid_ann_640", annulus(6.0, 10.0, 80, 8), far=[(12.0, 0.0, 0.0), (20.0, 0.0, 0.0)])
    write_solve_case("lid_ann_160_rev", annulus(6.0, 10.0, 40, 4, reverse=True))
    write_solve_case("lid_ann_tri", annulus(6.0, 10.0, 40, 4, tri=True))
    write_solve_case("lid_touch", annulus(5.0, 9.0, 40, 4))
    write_solve_case("lid_patches", annulus(6.0, 10.0, 40, 4, sectors=[(0.0, 60.0), (180.0, 240.0)]))
    write_solve_case("lid_shrink_1", annulus(6.0, 6.5, 40, 1))
    write_solve_case("lid_shrink_2", annulus(6.0, 6.25, 40, 1))
    write_solve_case("lid_deep", annulus(6.0, 10.0, 40, 4), depth=-1.0)
    write_solve_case("lid_limits", annulus(6.0, 10.0, 20, 2), limits=1)
    write_solve_case("lid_limits_off", annulus(6.0, 10.0, 20, 2), limits=1, lid_flag=0)

    # Field output (test_integration_lid.f90). lid_hull_160 keeps the hull centroids of image_cos as field points.
    ann = annulus(6.0, 10.0, 40, 4)
    write_case("lid_hull_160", ann)
    write_case("lid_field_off_centroid", ann, points=off_centroid_points(ann))

    # Symmetry and IRSP (test_integration_lid.f90). The half cases use the half hull of image_cos_half (y >= 0, ISY = 1) and the
    # half annulus of the same 40 x 4 mesh, whose panels at angles 0 to 180 degrees are the first 80 of the full one. The points
    # are the 10 first centroids, then points at y < 0 on the lid plane, outside it, and in the fluid.
    half_hull = os.path.join(ROOT, "test", "data", "image_term", "image_cos_half", "input")
    plane = os.path.join(ROOT, "test", "data", "image_term", "image_cos_lid", "input", "WaterPlaneMesh.pnl")
    full = annulus(6.0, 10.0, 40, 4)
    points = centroids(full, 10) + [(8.0, -3.0, 0.0), (-7.0, -4.0, 0.0), (12.0, -3.0, 0.0), (4.0, -2.0, -2.0), (-3.0, -4.0, -8.0),
                                    (7.0, -5.0, -3.0), (-9.0, 2.0, -1.0), (6.5, 1.0, -0.5)]
    write_case("lid_ann_160_pts", full, points=points)
    write_case("lid_ann_160_half", annulus(6.0, 10.0, 40, 4, isy=1, sectors=[(0.0, 180.0)]), points=points, hull=half_hull)
    write_solve_case("lid_ann_160_irsp", full, irsp=1, water_plane=plane)

    # Independent reference (test_integration_lid.f90, pre_postprocessing_python/DampingLid/annular_lid_reference.py): a
    # cylinder of radius 2 standing on the sea bed in 10 m of water, kh = 1.5, the lid r 2 to 4 at three meshes, epsilon
    # 0.2 m/s, and the finest mesh also at 1 and 1e6 m/s. The 48 x 8 cases have the 384 lid centroids as field points after
    # the two far points, for the lid power.
    cylinder = (wall(2.0, 10.0, 24, 10), wall_hydrostatics(2.0, 10.0))
    far = [(0.0, 8.0, 0.0), (8.0, 0.0, 0.0)]
    reference = dict(depth=10.0, wave_number=0.15, hull_files=cylinder)
    write_case("mm_cyl_nolid", None, lid_flag=0, points=far, **reference)
    write_case("mm_cyl_lid_a", annulus(2.0, 4.0, 24, 2), points=far, **reference)
    write_case("mm_cyl_lid_b", annulus(2.0, 4.0, 48, 4), points=far, **reference)
    lid_c = annulus(2.0, 4.0, 48, 8)
    points_c = far + centroids(lid_c, 384)
    write_case("mm_cyl_lid_c", lid_c, points=points_c, **reference)
    write_case("mm_cyl_lid_c_eps1", annulus(2.0, 4.0, 48, 8, dampers=[(1, 384, 1.0)]), points=points_c, **reference)
    write_case("mm_cyl_lid_c_rigid", annulus(2.0, 4.0, 48, 8, dampers=[(1, 384, 1.0e6)]), points=points_c, **reference)

    # Two bodies (test_integration_lid.f90). lid_two_far: the hulls of image_cos_two, the second 5 km away along y, with the lid
    # and field points of lid_ann_160 about the first. The stops: a lid disc inside the second body, and IRSP. lid_gap: two
    # boxes 20 x 8 x 4 m with 1 m panels, side by side 4 m apart (centres at y = -6 and 6), heading 90 degrees, with the lid on
    # the whole gap, 20 x 4 panels; field points at the 80 lid centroids, then three off the lid.
    write_case_two("lid_two_far", ann, centroids(ann, 10))
    write_case_two("lid_stop_inside_two", annulus(0.0, 2.5, 12, 2, centre=(0.0, 5000.0)), centroids(ann, 10))
    write_case_two("lid_stop_irsp_two", annulus(**outside), centroids(ann, 10), irsp=1, water_plane=plane)
    gap = rectangle(-10.0, 10.0, -2.0, 2.0, 20, 4)
    hull = (box(20.0, 8.0, 4.0, 20, 8, 4), box_hydrostatics(20.0, 8.0, 4.0))
    gap_points = centroids(gap, 80) + [(0.0, -15.0, 0.0), (0.0, 15.0, 0.0), (15.0, 0.0, 0.0)]
    gap_case = dict(hulls=[hull, hull], origins=[(0.0, -6.0), (0.0, 6.0)], heading=90.0)
    write_case_two("lid_gap", gap, gap_points, **gap_case)
    # The same with epsilon 0 and 1 m/s, and without the lid (Damping_lid 0, the lid file present and unused).
    write_case_two("lid_gap_eps0", rectangle(-10.0, 10.0, -2.0, 2.0, 20, 4, dampers=[(1, 80, 0.0)]), gap_points, **gap_case)
    write_case_two("lid_gap_eps1", rectangle(-10.0, 10.0, -2.0, 2.0, 20, 4, dampers=[(1, 80, 1.0)]), gap_points, **gap_case)
    write_case_two("lid_gap_nolid", gap, gap_points, lid_flag=0, **gap_case)
    # The same three for lid_gap_damps at 1.29 rad/s, the largest gap elevation of lid_gap_nolid in a sweep of 0.8 to 2.0 rad/s
    # (README, "Gap resonance"); the fixtures above stay at 1.35 rad/s.
    write_case_two("lid_gap_res", gap, gap_points, frequency=1.29, **gap_case)
    write_case_two("lid_gap_res_eps1", rectangle(-10.0, 10.0, -2.0, 2.0, 20, 4, dampers=[(1, 80, 1.0)]), gap_points, frequency=1.29,
                   **gap_case)
    write_case_two("lid_gap_res_nolid", gap, gap_points, lid_flag=0, frequency=1.29, **gap_case)


if __name__ == "__main__":
    main()
