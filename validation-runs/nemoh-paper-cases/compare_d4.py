"""D4: HAMS-MREL on Nemoh's test case 8a (cylinder R 6 m, T 14 m, deep water) with the test case's lid, our ring lid and
no lid, against HydroStar with irregular-frequency removal (the test case's reference, 0.0628-6.28 rad/s).

Usage: python -I compare_d4.py RESULTS_DIR HYDROSTAR_DIR OUT_DIR
Writes d4_irr8a.png and d4_irr8a.txt: |HAMS - HydroStar| / max|HydroStar| (max and median, %) per quantity and case, over
0.06-3 and 3-6.28 rad/s. HAMS is made dimensional with rho = 1025, g = 9.81 (A = rho a, B = rho omega b, X = rho g x).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import readers as R  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RHO, G = 1025.0, 9.81
CASES = (("d4_lid_nemoh", "HAMS-MREL, the test case's lid (132 panels)", dict(color="tab:red", lw=1.6)),
         ("d4_lid_ours", "HAMS-MREL, our ring lid (96 panels)", dict(color="tab:blue", lw=1.6)),
         ("d4_nolid", "HAMS-MREL, no lid (no irregular-frequency removal)", dict(color="tab:orange", lw=1.2, ls="--")))
SRC = "HydroStar with irregular-frequency removal (Nemoh 3.0.3 test case 8a)"
QS = (("A11", "AddedMass_11.dat", "A_11 [kg]"), ("A33", "AddedMass_33.dat", "A_33 [kg]"), ("A55", "AddedMass_55.dat", "A_55 [kg m²]"),
      ("B11", "Damping_11.dat", "B_11 [N s/m]"), ("B33", "Damping_33.dat", "B_33 [N s/m]"), ("B55", "Damping_55.dat", "B_55 [N m s]"),
      ("X1", "fxf1st.rao", "|X_1| [N/m]"), ("X3", "fzf1st.rao", "|X_3| [N/m]"), ("X5", "myf1st.rao", "|X_5| [N m/m]"))


def ref(rdir, fn):
    if fn.endswith(".rao"):
        w, a, _, _, _ = R.hs_rao(os.path.join(rdir, fn))
        return w, a[:, 0]
    rows = [ln.split() for ln in open(os.path.join(rdir, fn)) if ln.strip() and not ln.startswith("#")]
    d = np.array([[float(v) for v in r[:2]] for r in rows])
    return d[:, 0], d[:, 1]


def hams(folder, q):
    f = R.hams_first(os.path.join(folder, "output", "Wamit_format"))
    i = int(q[1])
    if q[0] == "A":
        h = f["AB"][(i, i)]
        return h[:, 0], h[:, 1] * RHO
    if q[0] == "B":
        h = f["AB"][(i, i)]
        return h[:, 0], h[:, 2] * RHO * h[:, 0]
    w, c = f["X"][(0.0, i)]
    return w, abs(c) * RHO * G


def main(res, rdir, out):
    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    lines = ["|HAMS - HydroStar| / max|HydroStar|, max / median in %, over 0.06-3.0 rad/s | 3.0-6.28 rad/s"]
    for k, (q, fn, sym) in enumerate(QS):
        a = ax[k // 3][k % 3]
        wr, vr = ref(rdir, fn)
        a.plot(wr, vr, color="k", ls="--", marker="o", ms=2.5, lw=1.0, label=SRC)
        row = "%-4s" % q
        for c, lab, st in CASES:
            p = os.path.join(res, c)
            if not os.path.isdir(p):
                continue
            w, v = hams(p, q)
            a.plot(w, v, label=lab, **st)
            ri = np.interp(w, wr, vr)
            d = abs(v - ri) / np.max(abs(vr)) * 100
            lo, hi = d[w <= 3.0], d[w > 3.0]
            row += "   %s %6.2f / %5.2f | %6.2f / %5.2f" % (c, lo.max(), np.median(lo), hi.max(), np.median(hi))
        lines.append(row)
        a.set_title(sym, fontsize=9)
        a.set_xlabel("ω [rad/s]", fontsize=8)
        a.grid(alpha=0.3)
        a.legend(fontsize=6)
    fig.suptitle("D4: irregular-frequency removal of HAMS-MREL on Nemoh's test case 8a (cylinder R 6 m, T 14 m, 502 hull "
                 "panels, deep water) against %s" % SRC, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "d4_irr8a.png"), dpi=95)
    with open(os.path.join(out, "d4_irr8a.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
