"""HAMS case of the VolturnUS-S semisubmersible (IEA Wind 15 MW reference turbine), matched to the second-order WAMIT run
of the IEA Wind Task 37 repository, with the free-surface integral on.

Reference data (kept outside the HAMS-MREL repository; the script reads them from a folder given on the command line):
IEA Wind Systems, IEA-15-240-RWT, release v1.1.18 (Barter et al. 2026, doi via GitHub release), Apache License 2.0,
OpenFAST/IEA-15-240-RWT-UMaineSemi/HydroData/:
  - wamit_inputs_1stOrder/IEA-15-240-RWT.gdf: the low-order panel mesh (2766 panels, half in y, "VolturnUS-S mesh from
    09-Mar-2020 rotated 180 deg 10-Sep-2022"). It is the mesh Roland et al. (2026) gave to every solver but Nemoh.
  - IEA-15-240-RWT-UMaineSemi_2ndOrder.out.wamit: the log of the second-order WAMIT run (v6.107S, 29 Feb 2024), whose
    .12d is the reference. The run used a higher-order (NURBS, 64 patches) geometry of the same hull and orientation, with a
    free-surface mesh (FDF: RINNER = 200 m), mass, external damping and mooring stiffness as below.

The model, as the WAMIT run (values from the .out.wamit and the 2nd-order .frc):
  - depth 200 m, rho = 1025 kg/m^3, g = 9.80665 m/s^2, heading 0;
  - mass 2.01e7 kg, centre of gravity (0, 0, -2.32) m, body mass matrix about the origin (IMASS, .frc);
  - external linear damping: diagonal 1.6299e5, 1.6299e5, 1.4632e6 N s/m, 1.2627e9, 1.2627e9, 0 N m s (IDAMP);
  - external stiffness (mooring): the full matrix of ISTIF;
  - hydrostatic and gravitational restoring of the WAMIT run (its .out: C33 = 446.700, C35 = 5.57620e-3,
    C44 = 267224, C46 = -0.180757, C55 = 267226, all over rho g), so that the restoring is that of the WAMIT hull, not of the
    low-order mesh (whose .hst gives C33 = 443.05, 0.8 % lower).
HAMS has rho = 1000 kg/m^3 built in: the mass, the damping and the mooring are scaled by 1000/1025 and the restoring is
taken as the WAMIT value over rho g times 1000 g. That keeps the motions and the loads over rho g of the 1025 kg/m^3 system
[inferred: every term of the equation of motion scales with rho], as for the OC4 WAMIT cases (make_mesh.oc4w_matrices).

Meshes written:
  - HullMesh.pnl: the IEA low-order mesh without its 228 interior waterplane (lid) panels, mirrored in y = 0 (full mesh,
    5076 panels), normals into the body. It is non-conforming where the columns meet the pontoons, which HAMS allows.
  - CSMesh_1.pnl: the control surface C, a vertical cylinder about the origin of radius C_RADIUS with a bottom at
    -(20 + C_BOTTOM_MARGIN) m, panels of about C_PANEL m (the waterlines are four circles, not star-shaped, so C and F come
    from files, as for OC4).
  - FSMesh_1.pnl: F, the free surface between the four column waterlines and the top ring of C (make_mesh.oc4_free_surface).
  - FSInnerMesh.pnl: the inner free surface of qtf_terms = full, between the waterlines and a circle of radius FS_RADIUS
    about the centre of the waterline bounding box (as ControlSurfaceMesh.f90 needs), triangles of about FS_PANEL m.

Frequencies: 0.25 to 1.90 rad/s in steps of 0.05 (34), on the 0.05 grid of the WAMIT run (0.25 to 2.50). The grid stops
below the first irregular frequency of the outer columns, about 1.94 rad/s [derived: k = j_01 / a with a = 6.25 m and
omega^2 = g k coth(k T), T = 20 m], since velocity_method = source rejects irregular-frequency removal.

Usage: python -I make_volturnus_case.py IEA_HYDRODATA_DIR OUT_CASE_DIR [THREADS]
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "pre_postprocessing_python", "SecondOrder"))
import make_mesh as mm  # noqa: E402

WAMIT_RHO = 1025.0
G = mm.G
DEPTH = 200.0
DRAFT = 20.0
FREQUENCIES = tuple(round(0.25 + 0.05 * i, 4) for i in range(34))
MASS = [[2.01e7, 0.0, 0.0, 0.0, -4.66e7, 0.0],
        [0.0, 2.01e7, 0.0, 4.66e7, 0.0, 0.0],
        [0.0, 0.0, 2.01e7, 0.0, 0.0, 0.0],
        [0.0, 4.66e7, 0.0, 4.26e10, 0.0, 0.0],
        [-4.66e7, 0.0, 0.0, 0.0, 4.26e10, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 2.37e10]]
ZG = -2.32
DAMPING = [1.6299e5, 1.6299e5, 1.4632e6, 1.2627e9, 1.2627e9, 0.0]
MOORING = [[7.192e4, 0.0, -3.587e-1, 0.0, 1.145e6, 0.0],
           [0.0, 7.192e4, 0.0, -1.145e6, 0.0, 2.165e1],
           [-3.587e-1, 0.0, 6.076e4, 0.0, -3.198e0, 0.0],
           [0.0, -1.145e6, 0.0, 2.587e8, 0.0, 3.995e2],
           [1.145e6, 0.0, -3.198e0, 0.0, 2.587e8, 0.0],
           [0.0, 2.165e1, 0.0, -1.908e1, 0.0, 2.524e8]]
# Hydrostatic and gravitational restoring of the second-order WAMIT run, over rho g (its .out.wamit); symmetric
RESTORING_OVER_RHO_G = {(2, 2): 446.700, (2, 4): 5.57620e-3, (3, 3): 267224.0, (3, 5): -0.180757, (4, 4): 267226.0}
C_RADIUS = 70.0          # m, C about the origin; the hull reaches 58.0 m from it
C_BOTTOM_MARGIN = 8.0    # m below the keel (z = -20 m)
C_PANEL = 4.0            # m, as OC4_CS
FS_RADIUS = 72.0         # m, inner free surface about the waterline box centre; the hull reaches about 65.5 m from it
FS_PANEL = 2.0           # m, as OC4_FS_PANEL (lambda / 8.5 at 1.90 rad/s)


def read_gdf(path):
    """Panels of a low-order WAMIT .gdf (ULEN = 1, ISX = 0) and its ISY flag; repeated vertices collapse to triangles."""
    lines = open(path).read().splitlines()
    ulen = float(lines[1].split()[0])
    isx, isy = (int(v) for v in lines[2].split()[:2])
    n = int(lines[3].split()[0])
    if ulen != 1.0 or isx != 0:
        raise ValueError("expected ULEN = 1 and ISX = 0")
    values = [float(v) for line in lines[4:4 + 4 * n] for v in line.split()[:3]]
    if len(values) < 12 * n:
        raise ValueError("the .gdf has fewer vertices than 4 x %d" % n)
    panels = []
    for k in range(n):
        pts = [tuple(values[12 * k + 3 * j:12 * k + 3 * j + 3]) for j in range(4)]
        unique = []
        for p in pts:
            if not unique or max(abs(a - b) for a, b in zip(p, unique[-1])) > 1e-9:
                unique.append(p)
        if len(unique) > 1 and max(abs(a - b) for a, b in zip(unique[0], unique[-1])) < 1e-9:
            unique.pop()
        if len(unique) < 3:
            raise ValueError("panel %d of the .gdf is degenerate" % (k + 1))
        panels.append(unique)
    return panels, isy


def tolerant_open_edges(mesh, tol=1e-9):
    """As make_mesh.open_edges, but an edge away from z = 0 used by one panel only is allowed: the IEA mesh is
    non-conforming where the columns meet the pontoons (z = -13 and -20 m), so panel edges there end on the side of a
    larger panel. HAMS treats such edges as non-conforming (SurfaceGeometry.f90, build_neighbours). Closure is checked
    separately by the projected areas (closure_check). An edge used twice in the same direction still stops."""
    use = {}
    for panel in mesh.panels:
        for k in range(len(panel)):
            a, b = panel[k], panel[(k + 1) % len(panel)]
            use.setdefault((min(a, b), max(a, b)), []).append((a, b))
    edges = []
    for (a, b), directed in use.items():
        if len(directed) > 2 or (len(directed) == 2 and directed[0] == directed[1]):
            raise ValueError("the mesh is not consistently oriented at the edge of nodes %d and %d" % (a, b))
        if len(directed) == 1 and abs(mesh.nodes[a - 1][2]) < tol and abs(mesh.nodes[b - 1][2]) < tol:
            edges.append(directed[0])
    return edges


def closure_check(mesh, props):
    """A closed hull (open at z = 0 only) has sum n A = (0, 0, waterplane area) for HAMS normals [divergence theorem]. Returns
    the largest |sum n_x A|, |sum n_y A| and |sum n_z A - area| over the waterplane area."""
    s = [0.0, 0.0, 0.0]
    for panel in mesh.panels:
        p = [mesh.nodes[i - 1] for i in panel]
        for k in range(1, len(p) - 1):
            c = mm.cross(mm.sub(p[k], p[0]), mm.sub(p[k + 1], p[0]))
            area = 0.5 * math.sqrt(sum(v * v for v in c))
            n = mm.hams_normal(p)
            for j in range(3):
                s[j] += n[j] * area
    return max(abs(s[0]), abs(s[1]), abs(abs(s[2]) - props["area"])) / props["area"]


def hull_mesh(gdf):
    """Full hull with HAMS normals (into the body): the .gdf panels without the interior waterplane panels (the lid of
    WAMIT's irregular-frequency removal, all vertices on z = 0), and their mirror in y = 0. The vertex order is reversed
    if needed so that the volume comes out positive (WAMIT normals point into the fluid). Returns the mesh, its properties
    and the number of lid panels removed from the half mesh."""
    panels, isy = read_gdf(gdf)
    if isy != 1:
        raise ValueError("expected a half mesh in y (ISY = 1)")
    lid = [p for p in panels if all(abs(v[2]) < 1e-9 for v in p)]
    panels = [p for p in panels if not all(abs(v[2]) < 1e-9 for v in p)]
    full = panels + [[(x, -y, z) for x, y, z in reversed(p)] for p in panels]
    for flip in (False, True):
        mesh = mm.Mesh()
        for p in full:
            mesh.add_panel(list(reversed(p)) if flip else p)
        props = mm.mesh_properties(mesh)
        if props["vol"] > 0.0:
            return mesh, props, len(lid)
    raise ValueError("no orientation gives a positive volume")


def waterline_loops(mesh):
    """The waterline edges (z = 0) chained into closed loops, each as counter-clockwise (x, y) points, with the centre and
    the mean radius of its nodes."""
    edges = tolerant_open_edges(mesh)
    nxt = {}
    for a, b in edges:
        nxt.setdefault(a, []).append(b)
    loops, seen = [], set()
    for a, _ in edges:
        if a in seen:
            continue
        loop, k = [], a
        while k not in seen:
            seen.add(k)
            loop.append(k)
            k = nxt[k][0]
        pts = [mesh.nodes[i - 1][:2] for i in loop]
        cx, cy = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
        pts.sort(key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
        r = sum(math.hypot(p[0] - cx, p[1] - cy) for p in pts) / len(pts)
        loops.append((pts, (cx, cy, r)))
    loops.sort(key=lambda item: (round(math.hypot(item[1][0], item[1][1]), 3), item[1][1]))
    return loops


def control_surface():
    depth = DRAFT + C_BOTTOM_MARGIN
    n = 16 * int(math.ceil(2.0 * math.pi * C_RADIUS / C_PANEL / 16.0))
    out = mm.Mesh()
    mm._wall(out, 0.0, 0.0, C_RADIUS, n, mm._wall_levels(depth, C_PANEL))
    rows = int(math.ceil(C_RADIUS / C_PANEL))
    radii = [C_RADIUS * (1.0 - i / rows) for i in range(rows + 1)]
    mm._disc(out, 0.0, 0.0, radii, mm.ring_counts(radii, n), -depth, up=True)
    return out, n


def write_control_file(path, threads):
    freqs = " ".join("%.4f" % w for w in FREQUENCIES)
    text = """   --------------HAMS Control file---------------

   Waterdepth  %.1fD0

   #Start Definition of Wave Frequencies
    0_inf_frequency_limits      0
    Input_frequency_type        3
    Output_frequency_type       3
    Number_of_frequencies      %d
  %s
   #End Definition of Wave Frequencies

   #Start Definition of Wave Headings
    Number_of_headings          1
  0.0000
   #End Definition of Wave Headings

    Reference_body_center          0.000       0.000       0.000
    Reference_body_length   1.0D0
    Wave_diffrac_solution    1
    If_remove_irr_freq       0
    Number of threads       %d

   #Start Definition of Pressure and/or Elevation (PE)
    Number_of_field_points     0
   #End Definition of Pressure and/or Elevation

   ----------End HAMS Control file---------------

   Note: generated by validation-runs/roland-benchmark/make_volturnus_case.py.
""" % (DEPTH, len(FREQUENCIES), freqs, threads)
    open(path, "w").write(text)


SECOND_ORDER = """! VolturnUS-S (IEA 15 MW), matched to the second-order WAMIT run of IEA-15-240-RWT v1.1.18; difference QTF with the
! free-surface integral, indirect and direct methods
mean_drift = 1
drift_methods = far near
motion_state = free
write_terms = 1
velocity_method = source
waterline_method = extrapolate
qtf = 1
qtf_terms = full
qtf_direct = 1
cs_shape = file
"""


def main(argv):
    hydro_dir, out = argv[1], argv[2]
    threads = int(argv[3]) if len(argv) > 3 else 24
    hull, props, nlid = hull_mesh(os.path.join(hydro_dir, "wamit_inputs_1stOrder", "IEA-15-240-RWT.gdf"))
    closure = closure_check(hull, props)
    if closure > 1e-6:
        raise ValueError("the hull is not closed: projected areas off by %.2e of the waterplane area" % closure)
    loops = waterline_loops(hull)
    if len(loops) != 4:
        raise ValueError("expected four column waterlines, found %d" % len(loops))
    polygons = [pts for pts, _ in loops]
    axes = [ax for _, ax in loops]
    reach_origin = max(math.hypot(x, y) for x, y, _ in hull.nodes)
    xs = [p[0] for poly in polygons for p in poly]
    ys = [p[1] for poly in polygons for p in poly]
    centre = (0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)))
    reach_centre = max(math.hypot(x - centre[0], y - centre[1]) for x, y, _ in hull.nodes)
    if not (C_RADIUS > reach_origin and FS_RADIUS > reach_centre and abs(centre[1]) < 1e-9):
        raise ValueError("C or the inner free surface cuts the hull: %.3f / %.3f m" % (reach_origin, reach_centre))
    control, c_count = control_surface()
    if len(mm.open_edges(control)) != c_count:
        raise ValueError("C is not closed below z = 0")
    free = mm.oc4_free_surface(polygons, C_RADIUS, c_count, C_PANEL, axes=axes)
    inner = mm.oc4_free_surface(polygons, FS_RADIUS, 16 * int(math.ceil(2 * math.pi * FS_RADIUS / FS_PANEL / 16.0)),
                                FS_PANEL, centre=centre, axes=axes)
    clearance = mm.check_control_clearance(hull, control)

    scale = mm.RHO / WAMIT_RHO
    mass = [[scale * v for v in row] for row in MASS]
    restoring = [[0.0] * 6 for _ in range(6)]
    for (i, j), v in RESTORING_OVER_RHO_G.items():
        restoring[i][j] = restoring[j][i] = v * mm.RHO * G
    mooring = [[scale * v for v in row] for row in MOORING]
    damping = {k: scale * v for k, v in enumerate(DAMPING) if v}

    inp = os.path.join(out, "input")
    os.makedirs(inp, exist_ok=True)
    mm.write_mesh(os.path.join(inp, "HullMesh.pnl"), hull)
    mm.write_hydrostatic(os.path.join(inp, "Hydrostatic.in"), dict(zg=ZG, mass=mass, restoring=restoring), 0.0, 0.0,
                         mooring, damping)
    mm.write_mesh(os.path.join(inp, "CSMesh_1.pnl"), control)
    mm.write_mesh(os.path.join(inp, "FSMesh_1.pnl"), free, title="Free Surface Mesh File")
    mm.write_mesh(os.path.join(inp, "FSInnerMesh.pnl"), inner, title="Free Surface Mesh File")
    # check_written_surfaces finds the waterline with open_edges, which needs a conforming hull: use the tolerant one
    strict = mm.open_edges
    mm.open_edges = tolerant_open_edges
    try:
        mm.check_written_surfaces(inp)
    finally:
        mm.open_edges = strict
    write_control_file(os.path.join(inp, "ControlFile.in"), threads)
    open(os.path.join(inp, "SecondOrder.in"), "w").write(SECOND_ORDER)
    # Unix line ends, whatever the platform (the inputs go to Snellius as they are)
    for name in os.listdir(inp):
        path = os.path.join(inp, name)
        data = open(path, "rb").read().replace(b"\r\n", b"\n")
        open(path, "wb").write(data)

    wamit = dict(vol=20206.3, area=446.700, zb=-13.625722)   # second-order WAMIT run, .out.wamit
    print("hull: %d lid panels of the half mesh removed; projected-area closure %.1e" % (nlid, closure))
    print("hull: %d panels, %d nodes; volume %.1f m^3 (WAMIT %.1f, %+.2f %%), waterplane %.2f m^2 (WAMIT %.2f, %+.2f %%), "
          "zb %.3f m (WAMIT %.3f)" % (len(hull.panels), len(hull.nodes), props["vol"], wamit["vol"],
                                      100 * (props["vol"] / wamit["vol"] - 1), props["area"], wamit["area"],
                                      100 * (props["area"] / wamit["area"] - 1), props["zb"], wamit["zb"]))
    for (pts, (x, y, r)) in loops:
        print("waterline: centre (%.3f, %.3f), radius %.3f m, %d nodes" % (x, y, r, len(pts)))
    print("C: %d panels, radius %.1f m (hull reach %.2f m), clearance %.2f; F: %d triangles; inner FS: %d triangles about "
          "(%.3f, 0), R %.1f m (hull reach %.2f m), smallest angle %.1f deg"
          % (len(control.panels), C_RADIUS, reach_origin, clearance, len(free.panels), len(inner.panels), centre[0],
             FS_RADIUS, reach_centre, mm.triangle_quality(inner)))


if __name__ == "__main__":
    main(sys.argv)
