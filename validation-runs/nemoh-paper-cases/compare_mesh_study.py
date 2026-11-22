"""Heave-plate mesh study: HAMS-MREL and Nemoh v3.0.3 on the mesh p0 of the paper's coarseness (and its refinements),
against the curves of Kurnia and Ducrozet (2023), Fig. 8 (their 3464-panel mesh) and HAMS-MREL on m3.

Usage: python -I compare_mesh_study.py FIG8_CSV RESULTS_DIR OUT_DIR TAG CASE [CASE ...]
Each CASE is a HAMS-MREL case folder in RESULTS_DIR; Nemoh's results of the same mesh are read from the folders
nemoh_<case>_c1 .. c7 when present (the 7 frequency chunks). Writes <TAG>.png and <TAG>.txt: per quantity, the error
|ours - WAMIT| / max|WAMIT| over 0.2-3.2 rad/s (max and median, %), and the mean over 3.3-4.0 rad/s divided by the mean of
Nemoh v3 of the paper there (the fall-off).
"""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_heave_plate as H  # noqa: E402
from same_mesh_m1 import read_fig8  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLOURS = ["tab:red", "tab:blue", "tab:orange", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:cyan", "tab:gray",
           "gold", "navy", "lime", "maroon"]


def nemoh_chunks(res, case):
    parts = sorted(glob.glob(os.path.join(res, "nemoh_%s_c*" % case)), key=lambda p: int(p.rsplit("_c", 1)[1]))
    if not parts:
        return None
    out = {}
    for p in parts:
        d = H.nemoh(p)
        for q, (w, v) in d.items():
            w0, v0 = out.get(q, (np.array([]), np.array([])))
            out[q] = (np.concatenate([w0, w]), np.concatenate([v0, v]))
    return out


def label(res, case):
    p = os.path.join(res, case, "label.txt")
    return open(p).read().strip() if os.path.exists(p) else case


def main(csv_path, res, out, tag, cases):
    ref = read_fig8(csv_path)
    hm3 = H.hams(os.path.join(res, "hp_m3"))
    runs = [(c, H.hams(os.path.join(res, c)), nemoh_chunks(res, c)) for c in cases]
    fig, ax = plt.subplots(3, 3, figsize=(16, 12))
    lines = ["error against WAMIT over 0.2-3.2 rad/s (max / median, %% of max|WAMIT|), and mean over 3.3-4.0 rad/s / that of "
             "Nemoh v3 (%s)" % H.PAPER]
    for k, q in enumerate(H.QUANT):
        a = ax[k // 3][k % 3]
        for src in ("WAMIT", "HydroStar", "Nemoh v3"):
            if (q, src) in ref:
                wr, vr = ref[(q, src)]
                a.plot(wr, vr, label="%s (%s)" % (src, H.PAPER), **H.STYLE[src])
        a.plot(hm3[q][0], hm3[q][1], color="tab:green", lw=1.0, marker="s", ms=2.5, label="HAMS-MREL m3 (11400 + 1240 lid)")
        ww, vw = ref[(q, "WAMIT")]
        wn, vn = ref[(q, "Nemoh v3")]
        nhi = vn[wn >= 3.3].mean()
        row = "%-4s" % q
        for n, (c, h, nm) in enumerate(runs):
            col = COLOURS[n % len(COLOURS)]
            for lab, d, st in (("HAMS-MREL", h, dict(color=col, lw=1.7)),
                               ("Nemoh v3.0.3 (our run)", nm, dict(color=col, lw=1.1, ls=":", marker=".", ms=3))):
                if d is None:
                    continue
                w, v = d[q]
                o = np.argsort(w)
                w, v = w[o], v[o]
                a.plot(w, v, label="%s, %s" % (lab, label(res, c)), **st)
                m = w <= 3.2
                e = abs(v[m] - np.interp(w[m], ww, vw)) / np.max(abs(vw)) * 100
                row += "  %s %s %5.1f/%4.1f fall %.2f" % (c, "H" if lab.startswith("HAMS") else "N", e.max(), np.median(e),
                                                       v[w >= 3.3].mean() / nhi)
        lines.append(row)
        a.set_title(H.SYM[q], fontsize=9)
        a.set_xlim(0, 4.05)
        a.set_xlabel("ω [rad/s]", fontsize=8)
        a.grid(alpha=0.3)
        a.legend(fontsize=5.2)
    fig.suptitle("Heave-plate cylinder, mesh study (%s): HAMS-MREL (solid) and Nemoh v3.0.3 (dotted, same colour = same mesh) "
                 "against %s" % (tag, H.PAPER), fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "%s.png" % tag), dpi=95)
    with open(os.path.join(out, "%s.txt" % tag), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
