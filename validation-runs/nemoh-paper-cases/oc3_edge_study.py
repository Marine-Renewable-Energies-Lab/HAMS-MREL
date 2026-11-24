"""Inputs of two OC3 Hywind checks of the difference-QTF rows omega_1 <= 0.25 rad/s (Nemoh test case 11, h = 320 m),
where HAMS-MREL differs from HydroStar by up to 65 % of max|QTF| (J5B and J5C).

(a) nemoh_oc3: Nemoh v3.0.3 on its own test case (mesh 4842, with its lid), as shipped except the frequency list: 20
    frequencies 0.0628-1.256 rad/s (step 0.0628, the grid of HydroStar and J5C) for the first order and the QTF. Every
    difference of the list is on the list, as the potential part (HASBO) needs.
(b) oc3_e0, oc3_e1, oc3_e2: HAMS-MREL on structured meshes of the same hull, refined only at the bottom edge (z = -120 m,
    r = 4.7 m), where the body velocity that the near-field terms II and III use is singular. e0 is uniform (rows of
    1.5 m on the wall, rings of about 0.53 m on the bottom, 56 panels round); e1 and e2 grade the wall rows and the
    bottom rings from an edge size of 0.375 m and 0.094 m (1/4 and 1/16 of the wall row) by 1.25 per row up to the
    uniform size. The rest of the hull is the same in the three. Same frequency list as (a). Settings of J5C
    (SecondOrder.in of oc3_diff_C), lid from the waterline as there, Hydrostatic.in of J5C (hydrosCal of the Nemoh mesh).
    The control surface C is held fixed at J5C's: cs_margin, cs_bottom_margin and cs_panel_size are set to the values
    HAMS-MREL derives from the median hull panel size of the Nemoh mesh, so that only the hull changes between e0-e2.

Usage: python -I oc3_edge_study.py NEMOH_TESTCASES J5C_INPUT_DIR OUT_DIR
Writes OUT_DIR/<case>/input (HAMS), OUT_DIR/nemoh_oc3 (Nemoh), OUT_DIR/oc3_edge_study.json and oc3_edge_meshes.png.
Nothing here runs HAMS or Nemoh.
"""
import json
import math
import os
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_cases as C  # noqa: E402
import nemoh_mesh  # noqa: E402
from nemoh_mesh import make_mesh as M  # noqa: E402

R_TOP, R_BOT, Z_T1, Z_T2, DRAFT = 3.25, 4.7, -4.0, -12.0, 120.0   # OC3 Hywind: the test case's mesh (rings, levels)
NTH = 56
H_TOP, H_TAPER, H_WALL, H_BOT = 0.8, 1.0, 1.5, 0.53
EDGE = (("oc3_e0", None), ("oc3_e1", H_WALL / 4), ("oc3_e2", H_WALL / 16))
FREQS = C.grid(20)                                                  # 0.0628 ... 1.256
GROWTH = 1.25


def cone(mesh, r1, z1, r2, z2, n, rows):
    """Frustum between (r1, z1) and (r2, z2), z1 > z2, in rows of quadrilaterals; HAMS normals into the body."""
    for i in range(rows):
        ra, za = r1 + (r2 - r1) * i / rows, z1 + (z2 - z1) * i / rows
        rb, zb = r1 + (r2 - r1) * (i + 1) / rows, z1 + (z2 - z1) * (i + 1) / rows
        for k in range(n):
            a0, a1 = 2 * math.pi * k / n, 2 * math.pi * (k + 1) / n
            pts = [(rb * math.cos(a0), rb * math.sin(a0), zb), (rb * math.cos(a1), rb * math.sin(a1), zb),
                   (ra * math.cos(a1), ra * math.sin(a1), za), (ra * math.cos(a0), ra * math.sin(a0), za)]
            cx, cy = sum(p[0] for p in pts) / 4, sum(p[1] for p in pts) / 4
            M._add_oriented(mesh, pts, (-cx, -cy, 0.0))


def hull(h_edge):
    mesh = M.Mesh()
    M._wall(mesh, 0.0, 0.0, R_TOP, NTH, uniform(0.0, Z_T1, H_TOP))
    cone(mesh, R_TOP, Z_T1, R_BOT, Z_T2, NTH, int(round(math.hypot(R_BOT - R_TOP, Z_T2 - Z_T1) / H_TAPER)))
    if h_edge:
        wall = C.graded(-DRAFT, Z_T2, h_edge, GROWTH, H_WALL)[::-1]          # top first
        bot = C.graded(R_BOT, 0.0, h_edge, GROWTH, H_BOT)
    else:
        wall = uniform(Z_T2, -DRAFT, H_WALL)
        bot = uniform(R_BOT, 0.0, H_BOT)
    M._wall(mesh, 0.0, 0.0, R_BOT, NTH, wall)
    M._disc(mesh, 0.0, 0.0, bot, M.ring_counts(bot, NTH), -DRAFT, up=True)
    return mesh, wall, bot


def uniform(a, b, h):
    n = max(1, int(round(abs(b - a) / h)))
    return [a + (b - a) * i / n for i in range(n + 1)]


def longest_edges(mesh):
    out = []
    for p in mesh.panels:
        pts = [mesh.nodes[i - 1] for i in p]
        out.append(max(float(np.linalg.norm(np.subtract(pts[k], pts[(k + 1) % len(pts)]))) for k in range(len(pts))))
    return np.array(out)


def median_size(mesh):
    """The median panel size of HAMS-MREL (the longest side, CalPanelChartLength in NormalProcess.f90), as the control
    surface takes it (ControlSurfaceMesh.f90: median of mesh%psize)."""
    return float(np.median(longest_edges(mesh)))


def second_order(j5c_input, med):
    keep = [ln.rstrip("\r\n") for ln in open(os.path.join(j5c_input, "SecondOrder.in"))
            if ln.strip() and not ln.startswith("!")]
    keep = [ln for ln in keep if not ln.startswith("qtf_quad_terms")]
    return keep + ["qtf_quad_terms = 1", "cs_margin = %.10f" % (4 * med), "cs_bottom_margin = %.10f" % (4 * med),
                   "cs_panel_size = %.10f" % med]


def nemoh_case(tc, out):
    src = os.path.join(tc, "11_QTF_OC3_Hywind", "mesh4842_Floating")
    dst = os.path.join(out, "nemoh_oc3")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    # every file with LF line ends (the test case ships CRLF; Nemoh runs on Linux), the Makefile left out
    os.makedirs(dst)
    for fn in os.listdir(src):
        if fn == "Makefile":
            continue
        text = open(os.path.join(src, fn)).read()       # universal newlines: CRLF read as \n
        with open(os.path.join(dst, fn), "w", newline="\n") as f:
            f.write(text)
    cal = open(os.path.join(src, "Nemoh.cal")).read().split("\n")
    w0, w1 = FREQS[0], FREQS[-1]
    n = 0
    for i, ln in enumerate(cal):
        t = ln.split()
        if ln.lstrip().startswith("1 65") and "Freq type" in ln:          # load cases: rad/s, 65 frequencies
            cal[i] = "1 %d\t%.4f\t\t%.4f\t\t\t! Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies, Min, and Max values" % (len(FREQS), w0, w1)
            n += 1
        elif len(t) > 2 and t[0] == "65" and "QTF computation" in ln:     # QTF frequencies
            cal[i] = "%d\t%.4f %.4f \t! Number of radial frequencies, Min, and Max values for the QTF computation" % (len(FREQS), w0, w1)
            n += 1
    if n != 2:
        raise RuntimeError("Nemoh.cal of test case 11: expected 2 frequency lines, changed %d" % n)
    with open(os.path.join(dst, "Nemoh.cal"), "w", newline="\n") as f:
        f.write("\n".join(cal))
    return dst


def main(tc, j5c_input, out):
    os.makedirs(out, exist_ok=True)
    facts = dict(freqs=FREQS)
    # J5C's hull (the Nemoh mesh without its lid) and the median panel size HAMS takes from it
    isym, nodes, panels = nemoh_mesh.read_nemoh(os.path.join(tc, "11_QTF_OC3_Hywind", "mesh4842_Floating", "OC3_hywind-4842_Nem.dat"))
    nh, _ = nemoh_mesh.split(nodes, panels)
    med = median_size(nh)
    facts["j5c"] = dict(hull=len(nh.panels), median_longest_side=round(med, 6))
    keys = second_order(j5c_input, med)
    meshes = {}
    for name, he in EDGE:
        h, wall, bot = hull(he)
        lid, lidinfo = C.lid_for(h, H_TOP)
        folder = C.write_hams_case(out, name, h, lid, 0,
                                   dict(depth=320.0, freqs=FREQS, headings=[0.0], centre=(0.0, 0.0, 0.0), irsp=1,
                                        note="OC3 Hywind, bottom-edge study (oc3_edge_study.py), edge size %s m."
                                        % ("uniform" if he is None else "%.4f" % he)),
                                   None, keys)
        shutil.copy(os.path.join(j5c_input, "Hydrostatic.in"), os.path.join(folder, "Hydrostatic.in"))
        props = M.mesh_properties(h)
        e = longest_edges(h)
        meshes[name] = (h, wall, bot)
        facts[name] = dict(hull=len(h.panels), lid=lidinfo["panels"], edge_size=he, wall_rows=len(wall) - 1,
                           bottom_rings=len(bot) - 1, first_wall_row=round(abs(wall[-1] - wall[-2]), 4),
                           first_bottom_ring=round(abs(bot[0] - bot[1]), 4), median_longest_side=round(float(np.median(e)), 4),
                           vol=round(props["vol"], 1), vol_nemoh=7982.4)
    facts["nemoh_oc3"] = nemoh_case(tc, out)
    with open(os.path.join(out, "oc3_edge_study.json"), "w") as f:
        json.dump(facts, f, indent=1)
    plot(meshes, nh, facts, out)
    print(json.dumps(facts, indent=1))


def plot(meshes, nh, facts, out):
    fig, ax = plt.subplots(1, 4, figsize=(17, 6.4), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.4]))
    cols = dict(oc3_e0="tab:blue", oc3_e1="tab:orange", oc3_e2="tab:red")
    # section through the bottom edge: the hull (black), wall row boundaries (ticks out from the wall) and bottom ring
    # boundaries (ticks down from the bottom)
    for k, (name, (h, wall, bot)) in enumerate(meshes.items()):
        a = ax[k]
        a.plot([0, R_BOT, R_BOT], [-DRAFT, -DRAFT, -114], color="k", lw=1.6,
               label="hull section through the bottom edge (r = 4.7 m, z = -120 m)")
        lab = "%s: panel boundaries (wall rows, bottom rings)" % name
        for z in wall:
            if z > -114:
                continue
            a.plot([R_BOT, R_BOT + 0.35], [z, z], color=cols[name], lw=0.9, label=lab)
            lab = None
        for r in bot:
            a.plot([r, r], [-DRAFT, -DRAFT - 0.35], color=cols[name], lw=0.9)
        a.set_xlim(0, 5.4)
        a.set_ylim(-120.8, -114)
        a.set_aspect("equal")
        f = facts[name]
        a.set_title("%s, HAMS-MREL: %d hull panels\nrow at the edge %.3f m (wall), %.3f m (bottom)"
                    % (name, f["hull"], f["first_wall_row"], f["first_bottom_ring"]), fontsize=9)
        a.set_xlabel("r [m]")
        a.set_ylabel("z [m]")
        a.legend(fontsize=6.3, loc="upper left")
    # panel extent normal to the edge (height of a wall panel, radial width of a bottom panel) against distance
    a = ax[3]
    for name, (h, wall, bot) in meshes.items():
        d, s = dist_size(h)
        a.loglog(d, s, ".", ms=3, color=cols[name], alpha=0.7, label="HAMS-MREL %s (%d panels)" % (name, facts[name]["hull"]))
    d, s = dist_size(nh)
    a.loglog(d, s, "x", ms=2.5, color="0.45", alpha=0.5, label="Nemoh test-case mesh, used in J5C (%d panels)" % len(nh.panels))
    a.set_xlabel("distance of the panel centroid from the bottom edge [m]")
    a.set_ylabel("panel extent normal to the edge [m] (wall: height, bottom: radial width)")
    a.set_xlim(0.02, 20)
    a.set_ylim(0.03, 3)
    a.grid(alpha=0.3, which="both")
    a.legend(fontsize=7, markerscale=3)
    a.set_title("Panel size against distance from the bottom edge (within 20 m)", fontsize=9)
    fig.suptitle("OC3 Hywind bottom-edge mesh study (HAMS-MREL inputs, not run): structured hull, rows graded by 1.25 "
                 "only next to the edge; 56 panels round in all three", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "oc3_edge_meshes.png"), dpi=100)


def dist_size(mesh):
    """Distance of each panel centroid from the bottom edge, and the panel extent normal to the edge: z extent on the
    wall and cone, radial extent on the bottom."""
    d, s = [], []
    for p in mesh.panels:
        pts = np.array([mesh.nodes[i - 1] for i in p])
        c = pts.mean(0)
        if c[2] < -DRAFT + 1e-6:
            r = np.hypot(pts[:, 0], pts[:, 1])
            s.append(np.ptp(r))
        else:
            s.append(np.ptp(pts[:, 2]))
        d.append(math.hypot(min(math.hypot(c[0], c[1]) - R_BOT, 0.0), c[2] + DRAFT))
    return np.array(d), np.array(s)


if __name__ == "__main__":
    main(*sys.argv[1:4])
