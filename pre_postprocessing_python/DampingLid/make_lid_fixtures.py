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


def control_text(lid, solver=None, depth=None, limits=None, points=None, irsp=0):
    """The image_cos control file with Damping_lid inserted after 'Number of threads', and Solver_type when given. depth
    replaces the water depth, limits the 0_inf_frequency_limits switch, and points the field points."""
    with open(os.path.join(SOURCE, "ControlFile.in"), "rb") as f:
        text = f.read().decode("ascii")
    out = []
    skip = False
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
    return EOL.join(out)


def write_case(name, lid_text, lid_flag=1, solver=None, depth=None, limits=None, points=None, hull=None, water_plane=None,
               irsp=0):
    """hull replaces the folder of the hull and hydrostatics (default image_cos), water_plane is the path of a WaterPlaneMesh.pnl
    to copy, used with irsp = 1."""
    folder = os.path.join(TARGET, name, "input")
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    for f in ("HullMesh.pnl", "Hydrostatic.in"):
        shutil.copy(os.path.join(hull or SOURCE, f), os.path.join(folder, f))
    if water_plane is not None:
        shutil.copy(water_plane, os.path.join(folder, "WaterPlaneMesh.pnl"))
    with open(os.path.join(folder, "ControlFile.in"), "wb") as f:
        f.write(control_text(lid_flag, solver, depth, limits, points, irsp).encode("ascii"))
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


if __name__ == "__main__":
    main()
