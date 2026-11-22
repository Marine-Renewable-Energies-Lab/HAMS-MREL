"""OC4 (Nemoh test case 9): where are the largest QTF errors inside dw <= 0.6 rad/s? For each heading pair and mode
(QTF_1, QTF_3, QTF_5), the worst pairs against HydroStar, with the first-order motions at both frequencies, and a map of the
error over (omega_1, dw) with the motion resonances and the frequency above which the mesh has fewer than 6 waterline
panels per wavelength.

Usage: python -I oc4_worst_pairs.py HAMS_WAMIT_FORMAT_DIR HYDROSTAR_DIR OUT_DIR
Writes oc4_worst_pairs.txt and oc4_worst_pairs.png. Error as in compare_case.py: |HAMS - HydroStar| / max|HydroStar|
over the pairs with dw <= 0.6 rad/s, in percent.
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

RG = R.RHO_HS * R.G_HS
GRID = [round(0.0628 * (k + 1), 4) for k in range(50)]
TAG = {1: "Fx", 3: "Fz", 5: "My"}
WL_PANEL = 3.1   # m, the waterline panel of the test-case mesh (2200 panels), from nemoh_mesh.summary
G = 9.81


def main(hdir, rdir, out):
    f = R.hams_first(hdir)
    lines = []
    # resonances: peaks of the HAMS motions at heading 0
    res = {}
    for i, nm in ((1, "surge"), (3, "heave"), (5, "pitch")):
        w, c = f["RAO"][(0.0, i)]
        k = np.argmax(abs(c))
        res[nm] = (w[k], abs(c[k]))
    lines.append("HAMS motion peaks (heading 0): " + ", ".join("%s %.3f rad/s (|RAO| %.3g)" % (k, v[0], v[1]) for k, v in res.items()))
    w6 = np.sqrt(2 * np.pi * G / (6 * WL_PANEL))
    lines.append("6 waterline panels (3.1 m) per wavelength at omega = %.2f rad/s" % w6)
    W1, W2 = np.meshgrid(GRID, GRID, indexing="ij")
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for r, heads in enumerate(((0.0, 0.0), (30.0, 0.0))):
        for c, i in enumerate((1, 3, 5)):
            q = R.hs_qtf(os.path.join(rdir, "MQTF%s.rao" % TAG[i]))
            ref = R.hs_pair_matrix(q, heads, GRID) / RG
            ham = R.hams_qtf(os.path.join(hdir, "QTF_Indirect.11d"), i, heads[0], heads[1], GRID)
            sel = np.isfinite(ref) & np.isfinite(ham) & (W2 >= W1) & (W2 - W1 <= 0.6 + 1e-6)
            sel &= ~((ref == 0) & (W2 > W1))
            sc = np.nanmax(abs(ref[sel]))
            e = np.where(sel, abs(ham - ref) / sc * 100, np.nan)
            order = np.argsort(np.where(sel, -e, np.inf), axis=None)[:6]
            lines.append("QTF-_%d (%g, %g): median %.2f %%, scale %.4g; worst pairs:" % (i, heads[0], heads[1], np.nanmedian(e), sc))
            for o in order:
                a, b = np.unravel_index(o, e.shape)
                xs = []
                for j in (a, b):
                    vals = []
                    for m in (1, 3, 5):
                        w, cc = f["RAO"][(heads[0] if j == a else heads[1], m)]
                        vals.append(abs(np.interp(GRID[j], w, cc.real) + 1j * np.interp(GRID[j], w, cc.imag)))
                    xs.append("|xi1,3,5|(%.3f) = %.3g, %.3g, %.3g" % (GRID[j], *vals))
                lines.append("   (%.3f, %.3f) dw %.3f: err %6.2f %%, |HydroStar| %.4g, |HAMS| %.4g; %s; %s"
                             % (GRID[a], GRID[b], GRID[b] - GRID[a], e[a, b], abs(ref[a, b]), abs(ham[a, b]), xs[0], xs[1]))
            a_ = ax[r][c]
            dw = W2 - W1
            m = np.isfinite(e)
            sc_ = a_.scatter(W1[m], dw[m], c=e[m], s=16, cmap="magma_r", vmin=0, vmax=max(5, np.nanpercentile(e, 99)))
            plt.colorbar(sc_, ax=a_, label="|HAMS-MREL − HydroStar| / max|HydroStar| [%]")
            for n_, (nm, (wr, _)) in enumerate(res.items()):
                a_.axvline(wr, color="tab:blue", lw=0.8, ls="--", label="ω_1 or ω_2 at a HAMS-MREL motion peak" if n_ == 0 else None)
                a_.plot(GRID, [wr - x for x in GRID], color="tab:blue", lw=0.8, ls="--")
                a_.text(wr, 0.61, nm, color="tab:blue", fontsize=7, ha="center")
            a_.axvline(w6, color="tab:green", lw=1.0, label="ω_1 or ω_2 = %.2f rad/s: 6 waterline panels per wavelength" % w6)
            a_.plot(GRID, [w6 - x for x in GRID], color="tab:green", lw=1.0)
            for n_, o in enumerate(order[:3]):
                aa, bb = np.unravel_index(o, e.shape)
                a_.plot(GRID[aa], GRID[bb] - GRID[aa], "o", mfc="none", mec="tab:cyan", ms=10, label="worst 3 pairs" if n_ == 0 else None)
            a_.set_ylim(-0.02, 0.72)
            a_.set_xlim(0, 3.2)
            a_.set_xlabel("ω_1 [rad/s]", fontsize=8)
            a_.set_ylabel("Δω = ω_2 − ω_1 [rad/s]", fontsize=8)
            a_.legend(fontsize=6.5, loc="upper right")
            a_.set_title("QTF⁻_%d, (β_1, β_2) = (%g°, %g°)" % (i, heads[0], heads[1]), fontsize=9)
    fig.suptitle("OC4 (Nemoh test case 9): error of the HAMS-MREL difference QTF against HydroStar inside Δω ≤ 0.6 rad/s; "
                 "colour = |HAMS-MREL − HydroStar| / max|HydroStar| [%]", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "oc4_worst_pairs.png"), dpi=95)
    with open(os.path.join(out, "oc4_worst_pairs.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
