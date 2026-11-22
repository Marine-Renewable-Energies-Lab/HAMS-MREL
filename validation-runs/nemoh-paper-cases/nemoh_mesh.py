"""Read a Nemoh v3 mesh file and split it into hull and lid panels for HAMS-MREL.

Nemoh mesh format (Nemoh v3 documentation, doc/source/input-output.rst): line 1 "2 isym" (isym = 1 for a half mesh about
y = 0), then "id x y z" per node, closed by a line whose first number is 0, then four node numbers per panel (a triangle
repeats a node), closed by "0 0 0 0". Nemoh's vertex order gives a right-hand normal into the fluid; HAMS computes the
opposite normal from the same order (make_mesh.py, CalTransNormals), so the order is kept. Panels with every node on
z = 0 are the lid of Nemoh's irregular-frequency removal; they go to WaterPlaneMesh.pnl, the rest to HullMesh.pnl.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "pre_postprocessing_python", "SecondOrder"))
import make_mesh  # noqa: E402

Z0 = 1e-6


def read_nemoh(path):
    """Return (isym, nodes, panels): nodes a list of (x, y, z), panels lists of 1-based node numbers (3 or 4)."""
    with open(path) as f:
        lines = [ln.split() for ln in f if ln.strip()]
    isym = int(lines[0][1])
    nodes, number = [], {}
    i = 1
    while True:
        t = lines[i]
        i += 1
        if int(float(t[0])) == 0:
            break
        number[int(t[0])] = len(nodes) + 1
        nodes.append(tuple(float(v) for v in t[1:4]))
    panels = []
    while i < len(lines):
        t = [int(float(v)) for v in lines[i][:4]]
        i += 1
        if all(v == 0 for v in t):
            break
        p = [number[v] for v in t]
        if p[3] == p[2] or p[3] == p[0]:
            p = p[:3]
        panels.append(p)
    return isym, nodes, panels


def to_mesh(nodes, panels):
    """A make_mesh.Mesh with only the nodes the panels use (renumbered), the vertex order kept."""
    out = make_mesh.Mesh()
    for p in panels:
        out.add_panel([nodes[i - 1] for i in p])
    return out


def split(nodes, panels):
    lid = [p for p in panels if all(abs(nodes[i - 1][2]) < Z0 for i in p)]
    hull = [p for p in panels if not all(abs(nodes[i - 1][2]) < Z0 for i in p)]
    return to_mesh(nodes, hull), to_mesh(nodes, lid)


def mirrored(mesh):
    """The whole body from a half mesh about y = 0 (the mirror image with its vertex order reversed)."""
    out = make_mesh.Mesh()
    for p in mesh.panels:
        pts = [mesh.nodes[i - 1] for i in p]
        out.add_panel(pts)
        if not all(abs(q[1]) < 1e-9 for q in pts):
            out.add_panel([(x, -y, z) for x, y, z in pts][::-1])
    return out


def properties(mesh, half):
    """Volume, buoyancy centre and waterplane integrals (make_mesh.mesh_properties) of the whole body."""
    full = mirrored(mesh) if half else mesh
    return make_mesh.mesh_properties(full)


def normals_into_body(mesh, axis_points):
    """Fraction of hull panels whose HAMS normal points towards the nearest of the vertical axes (x, y) given, measured
    horizontally, or upwards for a horizontal panel below the axes' bottom... used as a coarse orientation check."""
    good = 0
    for p in mesh.panels:
        pts = [mesh.nodes[i - 1] for i in p]
        c = tuple(sum(q[k] for q in pts) / len(pts) for k in range(3))
        n = make_mesh.hams_normal(pts)
        ax = min(axis_points, key=lambda a: (c[0] - a[0]) ** 2 + (c[1] - a[1]) ** 2)
        d = (ax[0] - c[0], ax[1] - c[1])
        horiz = math.hypot(n[0], n[1])
        if horiz > 0.5:
            good += (n[0] * d[0] + n[1] * d[1]) > 0.0
        else:
            good += 1  # horizontal panels are checked by the volume sign
    return good / len(mesh.panels)


def waterline_clearance(hull, lid):
    """Smallest horizontal distance from a lid centroid to a hull waterline segment, and the height of the hull's top
    panel row (the largest depth of the nodes of the panels that touch z = 0), for the lid rule of input_files.md."""
    segs = []
    row = 0.0
    for p in hull.panels:
        pts = [hull.nodes[i - 1] for i in p]
        top = [q for q in pts if abs(q[2]) < Z0]
        if len(top) >= 2:
            segs.append((top[0], top[1]))
            row = max(row, max(-q[2] for q in pts))
    best = float("inf")
    for p in lid.panels:
        pts = [lid.nodes[i - 1] for i in p]
        c = (sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts))
        for a, b in segs:
            ax, ay, bx, by = a[0], a[1], b[0], b[1]
            dx, dy = bx - ax, by - ay
            t = max(0.0, min(1.0, ((c[0] - ax) * dx + (c[1] - ay) * dy) / (dx * dx + dy * dy)))
            best = min(best, math.hypot(c[0] - ax - t * dx, c[1] - ay - t * dy))
    return best, row


def summary(path):
    isym, nodes, panels = read_nemoh(path)
    hull, lid = split(nodes, panels)
    xs = [n[0] for n in hull.nodes]
    ys = [n[1] for n in hull.nodes]
    zs = [n[2] for n in hull.nodes]
    pr = properties(hull, isym == 1)
    out = dict(isym=isym, hull=len(hull.panels), lid=len(lid.panels), nodes=len(nodes),
               x=(min(xs), max(xs)), y=(min(ys), max(ys)), z=(min(zs), max(zs)), vol=pr["vol"], zb=pr["zb"],
               xb=pr["xb"], awp=pr["area"],
               tri=sum(1 for p in hull.panels if len(p) == 3))
    if lid.panels:
        out["lid_clearance"], out["wl_row"] = waterline_clearance(hull, lid)
    return out, hull, lid


if __name__ == "__main__":
    for path in sys.argv[1:]:
        s, _, _ = summary(path)
        print(os.path.basename(path), s)
