"""
Generate a HAMS hull mesh and a matching Hydrostatic.in for a surface-piercing vertical cylinder or box.

Usage:
    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --out HullMesh.pnl
    python make_mesh.py box --lx 20 --ly 10 --draft 5 --nx 10 --ny 5 --nz 4 --out HullMesh.pnl

The hull has no lid, so the top row of nodes lies on z = 0. --x0 and --y0 shift the body in the horizontal plane. Every
generated panel normal is checked to point out of the fluid (into the body), computed as HAMS does in CalTransNormals.

Hydrostatic.in is written next to the mesh (or to --hydro-out) for a freely floating, uniform solid prism, using the
discretized waterplane polygon: mass = rho*V, centre of gravity at half draft, all matrices about the point (x0, y0, 0),
which must be given as the rotation centre in ControlFile.in. External damping and stiffness are zero.
"""

import argparse
import math
import os

RHO = 1000.0     # as Const_mod in src/WavDynMods.f90
G = 9.80665


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


def cylinder(radius, draft, ntheta, nz, nr, graded=False):
    """Side wall and flat bottom, with triangles at the centre of the bottom. Nodes lie on y = 0, so the mesh is
    symmetric about y = 0. With graded, rows are cosine spaced, finer at the waterline and the bilge. Returns the mesh
    and its waterline polygon, counter-clockwise from above."""
    if ntheta % 2 != 0:
        raise ValueError("ntheta must be even, to keep the mesh symmetric about both axes")
    mesh = Mesh()
    angles = [2.0 * math.pi * j / ntheta for j in range(ntheta + 1)]
    angles[-1] = 0.0

    def ring(r, j, z):
        return (r * math.cos(angles[j]), r * math.sin(angles[j]), z)

    def depth(i):
        return -draft * (0.5 * (1.0 - math.cos(math.pi * i / nz)) if graded else i / nz)

    def bottom_radius(i):
        return radius * (math.cos(0.5 * math.pi * i / nr) if graded else (nr - i) / nr)

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


def hydrostatics(polygon, draft):
    """Mass and restoring matrices about the waterplane centre of a uniform solid prism floating at its draft."""
    area, ixx, iyy = polygon_moments(polygon)
    vol = area * draft
    mass = RHO * vol
    zg = -draft / 2.0
    zb = -draft / 2.0
    # Inertia of the prism about the origin, by parallel axes from its centre of gravity
    i44 = RHO * draft * ixx + mass * draft * draft / 12.0 + mass * zg * zg
    i55 = RHO * draft * iyy + mass * draft * draft / 12.0 + mass * zg * zg
    i66 = RHO * draft * (ixx + iyy)

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


def write_mesh(path, mesh):
    with open(path, "w") as f:
        f.write("    --------------Hull Mesh File---------------\n\n")
        f.write("    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n")
        f.write("    %8d  %8d  %8d  %8d\n\n" % (len(mesh.panels), len(mesh.nodes), 0, 0))
        f.write("    # Start Definition of Node Coordinates     ! node_number   x   y   z\n")
        for number, (x, y, z) in enumerate(mesh.nodes, start=1):
            f.write("    %6d  %18.10f  %18.10f  %18.10f\n" % (number, x, y, z))
        f.write("    # End Definition of Node Coordinates\n\n")
        f.write("    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID ...\n")
        for number, panel in enumerate(mesh.panels, start=1):
            f.write("    %6d  %3d  " % (number, len(panel)) + "  ".join("%6d" % i for i in panel) + "\n")
        f.write("    # End Definition of Node Relations\n\n")
        f.write("    --------------End Hull Mesh File---------------\n")


def write_matrix(f, matrix):
    for row in matrix:
        f.write("".join("  %12.5E" % value for value in row) + "\n")


def write_hydrostatic(path, hydro, x0, y0):
    zero = [[0.0] * 6 for _ in range(6)]
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
        write_matrix(f, zero)


def main():
    parser = argparse.ArgumentParser(description="Generate a HAMS hull mesh and Hydrostatic.in.")
    sub_parsers = parser.add_subparsers(dest="shape")
    sub_parsers.required = True

    cyl = sub_parsers.add_parser("cylinder")
    cyl.add_argument("--radius", type=float, required=True)
    cyl.add_argument("--ntheta", type=int, required=True)
    cyl.add_argument("--nr", type=int, required=True)
    cyl.add_argument("--graded", action="store_true", help="cosine spacing, finer at the waterline and the bilge")

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
        p.add_argument("--out", required=True, help="path of the HullMesh file")
        p.add_argument("--hydro-out", help="path of the hydrostatic file (default: Hydrostatic.in next to --out)")

    args = parser.parse_args()
    if args.shape == "cylinder":
        mesh, polygon = cylinder(args.radius, args.draft, args.ntheta, args.nz, args.nr, args.graded)
    else:
        mesh, polygon = box(args.lx, args.ly, args.draft, args.nx, args.ny, args.nz)
    mesh = mesh.shifted(args.x0, args.y0)
    check_normals(mesh, (args.x0, args.y0, -args.draft / 2.0))

    hydro = hydrostatics(polygon, args.draft)
    hydro_path = args.hydro_out or os.path.join(os.path.dirname(args.out), "Hydrostatic.in")
    write_mesh(args.out, mesh)
    write_hydrostatic(hydro_path, hydro, args.x0, args.y0)

    print("%d panels, %d nodes, normals checked" % (len(mesh.panels), len(mesh.nodes)))
    print("V = %.6f, Awp = %.6f, zg = zb = %.6f" % (hydro["vol"], hydro["area"], hydro["zg"]))
    print("M = %.6E, I44 = %.6E, I55 = %.6E, I66 = %.6E" % (hydro["mass"][0][0], hydro["mass"][3][3],
                                                             hydro["mass"][4][4], hydro["mass"][5][5]))
    print("C33 = %.6E, C44 = %.6E, C55 = %.6E" % (hydro["restoring"][2][2], hydro["restoring"][3][3],
                                                   hydro["restoring"][4][4]))


if __name__ == "__main__":
    main()
