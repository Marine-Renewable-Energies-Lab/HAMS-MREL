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
TARGET = os.path.join(ROOT, "test", "data", "damping_lid")
EOL = "\r\n"


def annulus(r0, r1, nt, nr, tri=False, reverse=False, z=0.0, isx=0, isy=0, z_node=None, sectors=None):
    """Panels of the annulus r0 <= r <= r1 on z = z, nt around and nr radial.

    r0 = 0 gives a disc: the first row is nt triangles about a centre node. tri splits every quad in two. reverse flips the
    vertex order of every panel. z_node = (index, z) moves one node off the plane. sectors = [(a0, a1), ...] in degrees keeps
    only the panels whose middle angle lies in one of them; the other nodes stay in the file, unused. Returns the text of a
    .pnl file.
    """
    radii = [r0 + (r1 - r0) * i / nr for i in range(nr + 1)]
    nodes = []
    ring = []  # ring[i][j] is the node number of ring i at angle j, 1-based
    for i, r in enumerate(radii):
        if i == 0 and r0 == 0.0:
            nodes.append((0.0, 0.0, z))
            ring.append([len(nodes)] * nt)
            continue
        row = []
        for j in range(nt):
            a = 2.0 * math.pi * j / nt
            nodes.append((r * math.cos(a), r * math.sin(a), z))
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
    lines += ["    # End Definition of Node Relations", "", "    --------------End Damping Lid Mesh File---------------", ""]
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


def control_text(lid, solver=None, depth=None, limits=None, points=None, irsp=0, wave_number=None):
    """The image_cos control file with Damping_lid inserted after 'Number of threads', and Solver_type when given. depth
    replaces the water depth, limits the 0_inf_frequency_limits switch, points the field points, and wave_number the frequency
    by a wave number (Input_frequency_type 2)."""
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
               irsp=0, wave_number=None, hull_files=None):
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
        f.write(control_text(lid_flag, solver, depth, limits, points, irsp, wave_number).encode("ascii"))
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
    write_case("lid_flag_off", annulus(**outside), lid_flag=0)

    # Single-body solves (test_integration_lid.f90). Annulus r = 6 to 10 unless said.
    write_solve_case("lid_ann_40", annulus(6.0, 10.0, 20, 2))
    write_solve_case("lid_ann_160", annulus(6.0, 10.0, 40, 4))
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
    # cylinder of radius 2 standing on the sea bed in 10 m of water, kh = 1.5, the lid r 2 to 4 at two meshes.
    cylinder = (wall(2.0, 10.0, 24, 10), wall_hydrostatics(2.0, 10.0))
    far = [(0.0, 8.0, 0.0), (8.0, 0.0, 0.0)]
    reference = dict(depth=10.0, wave_number=0.15, points=far, hull_files=cylinder)
    write_case("mm_cyl_nolid", None, lid_flag=0, **reference)
    write_case("mm_cyl_lid_a", annulus(2.0, 4.0, 24, 2), **reference)
    write_case("mm_cyl_lid_b", annulus(2.0, 4.0, 48, 4), **reference)
    write_case("mm_cyl_lid_c", annulus(2.0, 4.0, 48, 8), **reference)


if __name__ == "__main__":
    main()
