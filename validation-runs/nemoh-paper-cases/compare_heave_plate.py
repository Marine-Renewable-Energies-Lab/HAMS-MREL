"""Heave-plate cylinder of Kurnia and Ducrozet (2023), Fig. 8: HAMS-MREL (own meshes hp_m1, hp_m2) and Nemoh v3.0.3 on hp_m1,
against the five curves of the paper (WAMIT, HydroStar, DIODORE, Nemoh v2, Nemoh v3 on the paper's 3464-panel mesh), read
from the PDF by extract_kurnia_fig8.py.

Usage: python -I compare_heave_plate.py FIG8_CSV RESULTS_DIR OUT_DIR

HAMS values are made dimensional with rho = 1025 kg/m^3 and g = 9.81 m/s^2 (the paper's): A = rho a, B = rho omega b,
|X| = rho g |x|; the rotation centre is the CG (0, 0, -7.5). Nemoh's RadiationCoefficients.tec and ExcitationForce.tec
are dimensional with rho = 1025 already. Errors: |HAMS - ref| / max |ref| over 0.2-4.0 rad/s, the reference interpolated
linearly to the HAMS frequencies (the paper's curves are polylines; WAMIT is points); max and median, in percent. The
paper's meshes are not ours, so these differences mix mesh and solver; Nemoh on hp_m1 is the same-mesh comparison.
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import readers as R  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RHO, G = 1025.0, 9.81
QUANT = ["A11", "A33", "A55", "B11", "B33", "B55", "X1", "X3", "X5"]
SYM = {"A11": "A_11 [kg]", "A33": "A_33 [kg]", "A55": "A_55 [kg m²]", "B11": "B_11 [N s/m]", "B33": "B_33 [N s/m]",
       "B55": "B_55 [N m s]", "X1": "|X_1| [N/m]", "X3": "|X_3| [N/m]", "X5": "|X_5| [N m/m]"}
PAPER = "Kurnia and Ducrozet (2023), Fig. 8"
# References in black and grey (markers and dashed lines); our runs as solid coloured lines
STYLE = {"WAMIT": dict(ls="", marker="+", color="k", ms=5), "HydroStar": dict(color="k", ls="--", lw=1.0),
         "DIODORE": dict(color="0.55", ls=":", lw=1.2), "Nemoh v2": dict(color="0.65", ls=(0, (6, 2)), lw=1.0),
         "Nemoh v3": dict(color="0.35", ls="-.", lw=1.0)}


def hams(folder):
    f = R.hams_first(os.path.join(folder, "output", "Wamit_format"))
    out = {}
    for i in (1, 3, 5):
        h = f["AB"][(i, i)]
        out["A%d%d" % (i, i)] = (h[:, 0], h[:, 1] * RHO)
        out["B%d%d" % (i, i)] = (h[:, 0], h[:, 2] * RHO * h[:, 0])
        w, c = f["X"][(0.0, i)]
        out["X%d" % i] = (w, abs(c) * RHO * G)
    return out


def tec_zones(path):
    zones, cur = [], None
    for ln in open(path):
        if ln.lstrip().startswith("Zone"):
            cur = []
            zones.append(cur)
        elif cur is not None and ln.strip():
            cur.append([float(v) for v in ln.split()])
    return [np.array(z) for z in zones]


def nemoh(folder):
    rad = tec_zones(os.path.join(folder, "results", "RadiationCoefficients.tec"))
    exc = tec_zones(os.path.join(folder, "results", "ExcitationForce.tec"))[0]
    out = {}
    for i in (1, 3, 5):
        z = rad[i - 1]
        out["A%d%d" % (i, i)] = (z[:, 0], z[:, 1 + 2 * (i - 1)])
        out["B%d%d" % (i, i)] = (z[:, 0], z[:, 2 + 2 * (i - 1)])
        out["X%d" % i] = (exc[:, 0], exc[:, 1 + 2 * (i - 1)])
    return out


def err(w, v, wr, vr):
    o = np.argsort(wr)
    wr, vr = np.asarray(wr)[o], np.asarray(vr)[o]
    m = (w >= wr.min() - 1e-9) & (w <= wr.max() + 1e-9)
    ri = np.interp(w[m], wr, vr)
    d = abs(v[m] - ri) / np.max(abs(vr)) * 100
    return d.max(), np.median(d)


def main(fig8, res, out):
    ref = defaultdict(list)
    for r in csv.DictReader(open(fig8)):
        ref[(r["quantity"], r["solver"])].append((float(r["omega_rad_s"]), float(r["value"])))
    h1, h2, nm = hams(os.path.join(res, "hp_m1")), hams(os.path.join(res, "hp_m2")), nemoh(os.path.join(res, "nemoh_hp_m1"))
    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    lines = ["Error |ours - ref| / max|ref| over 0.2-4.0 rad/s, max / median in %%; ref: the curves of %s" % PAPER]
    lines.append("%-5s %-22s %-16s %-16s %-16s" % ("", "reference", "HAMS m1", "HAMS m2", "Nemoh v3.0.3 m1"))
    for k, q in enumerate(QUANT):
        a = ax[k % 3][k // 3]
        for s, st in STYLE.items():
            p = sorted(ref[(q, s)])
            a.plot([x for x, _ in p], [y for _, y in p], label="%s (Kurnia & Ducrozet 2023, Fig. 8)" % s, **st)
        for lab, d, st in (("Nemoh v3.0.3 on our mesh m1 (this work)", nm, dict(color="tab:orange", lw=1.6)),
                           ("HAMS-MREL m1, 3096 panels (this work)", h1, dict(color="tab:red", lw=1.8)),
                           ("HAMS-MREL m2, 5632 panels (this work)", h2, dict(color="tab:blue", lw=1.3))):
            a.plot(d[q][0], d[q][1], label=lab, **st)
        a.set_title(SYM[q], fontsize=10)
        a.set_xlim(0, 4)
        a.set_xlabel("ω [rad/s]", fontsize=8)
        a.grid(alpha=0.3)
        a.legend(fontsize=5.5, loc="best")
        for s in ("WAMIT", "HydroStar", "Nemoh v3"):
            p = sorted(ref[(q, s)])
            wr, vr = [x for x, _ in p], [y for _, y in p]
            row = "%-5s %-22s" % (q, s)
            for d in (h1, h2, nm):
                e = err(np.asarray(d[q][0]), np.asarray(d[q][1]), wr, vr)
                row += " %6.2f / %5.2f   " % e
            lines.append(row)
    fig.suptitle("Truncated cylinder (R 10 m, T 10 m) with a thin heave plate (R 15 m, 0.1 m), deep water, CG (0, 0, -7.5) m:"
                 " HAMS-MREL and Nemoh v3.0.3 on our meshes against the curves of %s (the paper's 3464-panel mesh)" % PAPER,
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "heave_plate_fig8.png"), dpi=95)
    with open(os.path.join(out, "heave_plate_summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
