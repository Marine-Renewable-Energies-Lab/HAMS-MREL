"""Figure of the VolturnUS-S case of make_volturnus_case.py: the hull (3D view), and in plan the inner free-surface mesh,
F and the top ring of C. Usage: python -I plot_volturnus_case.py CASE_INPUT_DIR OUT_PNG"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "pre_postprocessing_python", "SecondOrder"))
import make_mesh as mm  # noqa: E402

INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
folder, out = sys.argv[1], sys.argv[2]
hull, control, free, inner = (mm.read_mesh(os.path.join(folder, n))
                              for n in ("HullMesh.pnl", "CSMesh_1.pnl", "FSMesh_1.pnl", "FSInnerMesh.pnl"))
polys = lambda m: [[m.nodes[i - 1] for i in p] for p in m.panels]
fig = plt.figure(figsize=(13, 5.6), facecolor=SURFACE)
ax = fig.add_subplot(1, 3, 1, projection="3d")
ax.add_collection3d(Poly3DCollection(polys(hull), facecolor="#86b6ef", edgecolor="#1c5cab", linewidths=0.15))
ax.set_xlim(-60, 35); ax.set_ylim(-52, 52); ax.set_zlim(-25, 5)
ax.set_box_aspect((95, 104, 30))
ax.view_init(25, -60)
ax.set_title("Hull: %d panels (IEA mesh, lid removed)" % len(hull.panels), loc="left", fontsize=9)
for k, (m, title) in enumerate(((inner, "Inner free surface: %d triangles, R = 72 m" % len(inner.panels)),
                                (free, "F: %d triangles; C: %d panels, R = 70 m" % (len(free.panels), len(control.panels))))):
    a = fig.add_subplot(1, 3, k + 2)
    a.add_collection(PolyCollection([[p[:2] for p in q] for q in polys(m)], facecolor="none", edgecolor=INK2, linewidths=0.2))
    wl = [[p[:2] for p in q] for q in polys(hull) if max(abs(v[2]) for v in q) < 1e-6 or min(v[2] for v in q) > -1.5]
    a.add_collection(PolyCollection(wl, facecolor="#2a78d6", edgecolor="#2a78d6", linewidths=0.4))
    a.set_aspect("equal"); a.autoscale()
    a.set_title(title, loc="left", fontsize=9)
    a.set_xlabel("x (m)"); a.set_ylabel("y (m)")
fig.suptitle("VolturnUS-S case of HAMS, matched to the second-order WAMIT run of IEA-15-240-RWT v1.1.18 (Apache-2.0)",
             x=0.01, ha="left", color=INK)
fig.tight_layout()
fig.savefig(out, dpi=150)
