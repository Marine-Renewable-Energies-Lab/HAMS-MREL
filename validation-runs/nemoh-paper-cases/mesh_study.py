"""Inputs of the heave-plate mesh study (Kurnia and Ducrozet 2023, Sec. 4.2.3, Fig. 8): the lid test D3 on the mesh m1,
a mesh of the paper's coarseness (p0) and its one-part-at-a-time refinements, for HAMS-MREL and Nemoh v3.0.3.

Usage: python -I mesh_study.py NEMOH_TESTCASES OUT_DIR {d3|paper}

Nothing here runs HAMS or Nemoh. OUT_DIR gets one folder per case (<case>/input for HAMS, nemoh_<case>/ for Nemoh) and
mesh_study.json with the panel counts and sizes.

The paper's mesh (its Fig. 3; GMSH, unstructured quadrilaterals): 3464 panels, 262 of them on the lid, so 3202 on the hull.
The text gives no split over the surfaces. The figure shows near-square panels of one size on the wall and the plate: the
hull's wetted area, 1731 m^2, over 3202 panels is 0.735 m per panel, about 14 rows on the 9.9 m wall; the lid, 314 m^2
over 262 panels, is 1.09 m per panel. p0 copies that with rings of near-square panels: every ring takes its own count
round(2 pi r / h), so neighbouring rings do not share nodes (a constant-panel method needs no conforming mesh). The plate
rim (0.1 m high) has one row. The top annulus and the bottom of the plate share their ring radii and counts on 10 m <= r
<= 15 m ("aligned": each top panel right above a bottom panel), unless offset is set (the bottom turned by half a panel
and its rings moved by half a ring).
"""
import json
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_cases as C  # noqa: E402
from nemoh_mesh import make_mesh as M  # noqa: E402

R_COL, R_PLATE, DRAFT, THICK, ZG = 10.0, 15.0, 10.0, 0.1, -7.5
TOP = DRAFT - THICK


def uniform(a, b, h):
    n = max(1, int(round(abs(b - a) / h)))
    return [a + (b - a) * i / n for i in range(n + 1)]


def ring_panels(mesh, ra, rb, n, z, up, offset=0.0):
    """One ring of n panels between radii ra and rb at z (quadrilaterals, or a fan if rb = 0); HAMS normal +z if up."""
    want = (0.0, 0.0, 1.0 if up else -1.0)

    def p(r, k):
        a = offset + 2.0 * math.pi * k / n
        return (r * math.cos(a), r * math.sin(a), z)
    for k in range(n):
        if rb <= 0.0:
            M._add_oriented(mesh, [p(ra, k), p(ra, k + 1), (0.0, 0.0, z)], want)
        else:
            M._add_oriented(mesh, [p(ra, k), p(ra, k + 1), p(rb, k + 1), p(rb, k)], want)


def rings(mesh, radii, h, z, up, offset=0.0, counts=None, nmin=6):
    """Rings between successive radii, each with round(2 pi r_mid / h) panels (or the given counts)."""
    used = []
    for i in range(len(radii) - 1):
        ra, rb = radii[i], radii[i + 1]
        n = counts[i] if counts else max(nmin, int(round(2.0 * math.pi * 0.5 * (ra + rb) / h)))
        dth = math.pi / n if offset else 0.0
        ring_panels(mesh, ra, rb, n, z, up, dth)
        used.append(n)
    return used


def wall(mesh, r, n, levels):
    M._wall(mesh, 0.0, 0.0, r, n, levels)


def plate_body(h=0.735, h_wall=None, h_wl=None, h_top=None, h_bot=None, h_edge=None, offset=False):
    """The heave-plate cylinder from rings of near-square panels.
    h: base panel size. h_wall: wall panel size (default h). h_wl: first wall row below the waterline, the rows growing by
    1.25 to h_wall (default: uniform rows). h_top, h_bot: plate top annulus and bottom disc panel size (default h).
    h_edge: first radial row at the plate edge, top and bottom, growing by 1.25 to the plate size (default: uniform)."""
    h_wall = h_wall or h
    h_top = h_top or h
    h_bot = h_bot or h
    mesh = M.Mesh()
    n_wall = int(round(2.0 * math.pi * R_COL / h_wall))
    levels = C.graded(0.0, -TOP, h_wl, hmax=h_wall) if h_wl else uniform(0.0, -TOP, h_wall)
    wall(mesh, R_COL, n_wall, levels)
    # top annulus (15 -> 10)
    top_r = C.graded(R_PLATE, R_COL, h_edge, hmax=h_top) if h_edge else uniform(R_PLATE, R_COL, h_top)
    top_counts = rings(mesh, top_r, h_top, -TOP, False)
    # bottom under the annulus: the same rings (aligned) if the sizes agree and offset is not set, else its own rings;
    # offset: ring edges half-way between the top's (the plate edge and r = 10 kept) and turned by half a panel
    if offset:
        bot_r = [top_r[0]] + [0.5 * (top_r[i] + top_r[i + 1]) for i in range(1, len(top_r) - 2)] + [top_r[-1]]
        rings(mesh, bot_r, h_bot, -DRAFT, True, offset=True)
    elif h_bot == h_top:
        rings(mesh, top_r, h_bot, -DRAFT, True, counts=top_counts)
    else:
        bot_r = C.graded(R_PLATE, R_COL, h_edge, hmax=h_bot) if h_edge else uniform(R_PLATE, R_COL, h_bot)
        rings(mesh, bot_r, h_bot, -DRAFT, True)
    rings(mesh, uniform(R_COL, 0.0, h_bot), h_bot, -DRAFT, True)
    # rim: one row, the count of the outermost top ring
    wall(mesh, R_PLATE, top_counts[0], [-TOP, -DRAFT])
    return mesh


def lid(h_lid=1.09, r=R_COL):
    mesh = M.Mesh()
    rings(mesh, uniform(r, 0.0, h_lid), h_lid, 0.0, False)
    return mesh


def sizes(mesh, select):
    out = []
    for p in mesh.panels:
        pts = [mesh.nodes[i - 1] for i in p]
        if not select(pts):
            continue
        edges = [math.sqrt(sum((a - b) ** 2 for a, b in zip(pts[k], pts[(k + 1) % len(pts)]))) for k in range(len(pts))]
        edges = [e for e in edges if e > 1e-9]   # a triangle repeats a node
        out.append((max(edges), min(edges)))
    if not out:
        return None
    return dict(n=len(out), longest=round(max(e[0] for e in out), 3), median_long=round(sorted(e[0] for e in out)[len(out) // 2], 3),
                max_aspect=round(max(e[0] / e[1] for e in out), 1))


def parts(mesh):
    def r(pts):
        return max(math.hypot(q[0], q[1]) for q in pts)
    z = lambda pts: [q[2] for q in pts]  # noqa: E731
    return dict(
        wall=sizes(mesh, lambda p: abs(r(p) - R_COL) < 1e-6 and max(z(p)) > -TOP + 1e-6 and min(z(p)) < -1e-9 and max(abs(q[2] - p[0][2]) for q in p) > 1e-9),
        top=sizes(mesh, lambda p: all(abs(q[2] + TOP) < 1e-9 for q in p)),
        bottom=sizes(mesh, lambda p: all(abs(q[2] + DRAFT) < 1e-9 for q in p)),
        rim=sizes(mesh, lambda p: abs(r(p) - R_PLATE) < 1e-6 and min(z(p)) < -TOP - 1e-9 and max(z(p)) > -DRAFT + 1e-9
                  and max(abs(q[2] - p[0][2]) for q in p) > 1e-9))


def hydro_m1(root):
    return os.path.join(root, "hp_m1", "input", "Hydrostatic.in")


def write_case(root, name, hull, lidm, freqs, note, irsp=1):
    folder = C.write_hams_case(root, name, hull, lidm, 0,
                               dict(depth=-1.0, freqs=freqs, headings=[0.0], centre=(0.0, 0.0, ZG), irsp=irsp, note=note), None)
    shutil.copy(hydro_m1(root), os.path.join(folder, "Hydrostatic.in"))
    props = M.mesh_properties(hull)
    info = dict(hull=len(hull.panels), lid=len(lidm.panels) if lidm else 0, vol=round(props["vol"], 2), nfreq=len(freqs),
                wmin=freqs[0], wmax=freqs[-1], parts=parts(hull))
    if lidm:
        info["lid_sizes"] = sizes(lidm, lambda p: True)
    return info


def write_nemoh(root, tc, name, hull, lidm, nfreq, wmin, wmax):
    nd = os.path.join(root, "nemoh_" + name)
    os.makedirs(nd, exist_ok=True)
    both = M.Mesh()
    for msh in (hull, lidm):
        for p in msh.panels:
            both.add_panel([msh.nodes[i - 1] for i in p])   # same orders as nemoh_hp_m1 (checked against test cases 8a, 8b, 9)
    M.write_nemoh_mesh(os.path.join(nd, name + ".dat"), both, 0)
    C.write_nemoh_cal(os.path.join(nd, "Nemoh.cal"), name + ".dat", len(both.nodes), len(both.panels), 0.0,
                      (0.0, 0.0, ZG), nfreq, wmin, wmax, (1, 0.0, 0.0))
    C.write_mesh_cal(os.path.join(nd, "Mesh.cal"), name, (0.0, 0.0, ZG), len(both.panels))
    shutil.copy(os.path.join(tc, "10_QTF_SOFTWIND", "mesh_1872half_QTF_Floating", "input_solver.txt"), nd)
    return dict(panels=len(both.panels), nfreq=nfreq)


def d3(tc, root):
    """hp_m1 at 2.0-4.0 rad/s by 0.05: (a) its lid (1 m rings, 744 panels), (b) no lid (IRSP = 0), (c) a lid of 0.35 m."""
    kw = dict(ntheta=96, nz=14, nr_plate=7, nr_bottom=12)
    hull = C.heave_plate(**kw)
    freqs = [round(2.0 + 0.05 * k, 4) for k in range(41)]
    lid_a, _ = C.lid_for(hull, 1.0)
    facts = {}
    facts["d3_lid_m1"] = write_case(root, "d3_lid_m1", hull, lid_a, freqs, "D3: hp_m1, its lid (1 m rings).")
    facts["d3_nolid"] = write_case(root, "d3_nolid", hull, None, freqs, "D3: hp_m1 without a lid.", irsp=0)
    facts["d3_lid035"] = write_case(root, "d3_lid035", hull, lid(0.35), freqs, "D3: hp_m1, lid of 0.35 m panels.")
    return facts


H0, HLID0 = 0.75, 1.09   # p0: hull 3230 panels (the paper's 3202), lid 259 (262)
FREQS = [round(0.2 + 0.05 * k, 4) for k in range(77)]   # 0.2-4.0 rad/s by 0.05
NEMOH_CHUNKS = 7   # Nemoh runs on one thread: 7 runs of 11 frequencies side by side


def with_nemoh(root, tc, name, hull, lidm, facts):
    """Nemoh inputs of the case in NEMOH_CHUNKS folders nemoh_<name>_c<k>, 11 frequencies each."""
    per = len(FREQS) // NEMOH_CHUNKS
    for k in range(NEMOH_CHUNKS):
        f = FREQS[k * per:(k + 1) * per]
        facts["nemoh_%s_c%d" % (name, k + 1)] = write_nemoh(root, tc, "%s_c%d" % (name, k + 1), hull, lidm, len(f), f[0], f[-1])


VARIANTS = {
    # name: (plate_body keywords, lid size, what it tests)
    "p0": (dict(), HLID0, "the paper's coarseness: near-square panels of 0.75 m, lid 1.09 m"),
    "p_lid": (dict(), 0.375, "lid halved to 0.375 m"),
    "p_wl": (dict(h_wl=0.15), HLID0, "wall rows graded from 0.15 m at the waterline to 0.75 m"),
    "p_wall": (dict(h_wall=0.375), HLID0, "wall halved to 0.375 m"),
    "p_plate": (dict(h_top=0.375, h_bot=0.375), HLID0, "plate top and bottom halved to 0.375 m (3.75 thicknesses)"),
    "p_edge": (dict(h_edge=0.1), HLID0, "plate rows graded from 0.1 m (one thickness) at the edge to 0.75 m, top and bottom"),
    "p_offset": (dict(offset=True), HLID0, "top and bottom panels of the plate offset by half a panel, not aligned"),
}


def paper(tc, root):
    """p0, the paper's coarseness, on 0.2-4.0 rad/s by 0.05 for HAMS and Nemoh."""
    facts = {}
    kw, hl, what = VARIANTS["p0"]
    hull, lidm = plate_body(H0, **kw), lid(hl)
    facts["p0"] = write_case(root, "p0", hull, lidm, FREQS, "Heave plate, %s." % what)
    with_nemoh(root, tc, "p0", hull, lidm, facts)
    return facts


def refine(tc, root):
    """The one-part-at-a-time refinements of p0, for HAMS and Nemoh (written for the cost table; staged only after the
    approval)."""
    facts = {}
    for name, (kw, hl, what) in VARIANTS.items():
        if name == "p0":
            continue
        hull, lidm = plate_body(H0, **kw), lid(hl)
        facts[name] = write_case(root, name, hull, lidm, FREQS, "Heave plate, %s." % what)
        facts[name]["what"] = what
        with_nemoh(root, tc, name, hull, lidm, facts)
    return facts


def irr8a(tc, root):
    """D4: the irregular-frequency removal of HAMS-MREL on Nemoh's test case 8a (cylinder R 6 m, T 14 m, deep water, 502
    hull panels; HydroStar with IRR to 6.28 rad/s), 100 frequencies 0.0628-6.28 rad/s: with the test case's own lid (132
    panels), with our ring lid, and without a lid. A, B and X only; the mass and restoring are placeholders (rho V, the
    hull's hydrostatics) and enter only the motions."""
    import nemoh_mesh as N
    path = os.path.join(tc, "8a_Cylinder_irregfreq", "mesh634", "OC5_cylinder634_Nem.dat")
    _, nodes, panels = N.read_nemoh(path)
    hull, lid_n = N.split(nodes, panels)
    s, _, _ = N.summary(path)
    lid_o, _ = C.lid_for(hull, s["wl_row"])
    freqs = C.grid(100)
    vol, zb, awp = s["vol"], s["zb"], s["awp"]
    m = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        m[i][i] = C.RHO_H * vol
    m[3][3] = m[4][4] = m[5][5] = C.RHO_H * vol * 25.0
    k = [[0.0] * 6 for _ in range(6)]
    k[2][2] = C.RHO_H * C.G_H * awp
    k[3][3] = k[4][4] = C.RHO_H * C.G_H * (math.pi * 6.0 ** 4 / 4.0 + vol * zb)
    hydro = dict(cg=(0.0, 0.0, 0.0), mass=m, restoring=k)
    facts = {}
    for name, lidm, irsp, note in (("d4_lid_nemoh", lid_n, 1, "test case 8a with its own lid"),
                                   ("d4_lid_ours", lid_o, 1, "test case 8a with our ring lid"),
                                   ("d4_nolid", None, 0, "test case 8a without a lid")):
        C.write_hams_case(root, name, hull, lidm, 0, dict(depth=-1.0, freqs=freqs, headings=[0.0], centre=(0.0, 0.0, 0.0),
                                                       irsp=irsp, note="D4: " + note + "."), hydro)
        facts[name] = dict(hull=len(hull.panels), lid=len(lidm.panels) if lidm else 0, nfreq=len(freqs))
        if lidm:
            facts[name]["lid_clearance_rows"] = round(N.waterline_clearance(hull, lidm)[0] / s["wl_row"], 3)
    return facts


def main():
    tc, root, what = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(root, exist_ok=True)
    facts = {"d3": d3, "paper": paper, "refine": refine, "irr8a": irr8a}[what](tc, root)
    path = os.path.join(root, "mesh_study_%s.json" % what)
    with open(path, "w") as f:
        json.dump(facts, f, indent=1)
    for k, v in facts.items():
        print(k, json.dumps(v))


if __name__ == "__main__":
    main()
