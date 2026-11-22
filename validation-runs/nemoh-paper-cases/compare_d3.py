"""D3, the lid test of the heave-plate cylinder above 2 rad/s: HAMS-MREL on mesh m1 with its lid (744 panels, 1 m rings),
without a lid, and with a 0.35 m lid (2602 panels), 2.0-4.0 rad/s by 0.05, against the curves of Kurnia and Ducrozet
(2023), Fig. 8 (their 3464-panel mesh), Nemoh v3.0.3 on m1 (our run) and HAMS-MREL on m3.

Usage: python -I compare_d3.py FIG8_CSV RESULTS_DIR OUT_DIR
Writes d3_lid.png and d3_lid.txt (per quantity and case: mean over 2.0-3.2 and 3.3-4.0 rad/s, relative to WAMIT's mean
over 2.0-3.2 rad/s, the range WAMIT covers).
"""
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

CASES = (("d3_lid_m1", "HAMS-MREL m1, its lid (744 panels, 1 m)", dict(color="tab:red", lw=1.8)),
         ("d3_nolid", "HAMS-MREL m1, no lid (no irregular-frequency removal)", dict(color="tab:orange", lw=1.3, ls="--")),
         ("d3_lid035", "HAMS-MREL m1, lid of 0.35 m (2602 panels)", dict(color="tab:blue", lw=1.8)))


def main(csv_path, res, out):
    ref = read_fig8(csv_path)
    runs = {c: H.hams(os.path.join(res, c)) for c, _, _ in CASES if os.path.isdir(os.path.join(res, c))}
    nm1 = H.nemoh(os.path.join(res, "nemoh_hp_m1"))
    hm3 = H.hams(os.path.join(res, "hp_m3"))
    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    lines = ["mean value over 2.0-3.2 and 3.3-4.0 rad/s, divided by WAMIT's mean over 2.0-3.2 rad/s (%s)" % H.PAPER]
    for k, q in enumerate(H.QUANT):
        a = ax[k // 3][k % 3]
        for src in ("WAMIT", "HydroStar", "DIODORE", "Nemoh v3"):
            if (q, src) in ref:
                wr, vr = ref[(q, src)]
                a.plot(wr, vr, label="%s (%s)" % (src, H.PAPER), **H.STYLE[src])
        a.plot(nm1[q][0], nm1[q][1], color="tab:purple", lw=1.2, marker="o", ms=3, label="Nemoh v3.0.3 (our run), m1 and its lid")
        a.plot(hm3[q][0], hm3[q][1], color="tab:green", lw=1.2, marker="s", ms=3, label="HAMS-MREL m3 (11400 + 1240 lid), 0.2 rad/s grid")
        wr, vr = ref[(q, "WAMIT")]
        wam = np.mean(vr[(wr >= 2.0) & (wr <= 3.2)]) if np.any((wr >= 2.0) & (wr <= 3.2)) else np.nan
        row = "%-4s" % q
        for c, lab, st in CASES:
            if c not in runs:
                continue
            w, v = runs[c][q]
            a.plot(w, v, label=lab, **st)
            lo, hi = v[(w >= 2.0) & (w <= 3.2)].mean(), v[(w >= 3.3)].mean()
            row += "  %s %.3f / %.3f" % (c, lo / wam, hi / wam)
        lines.append(row)
        a.set_xlim(1.9, 4.05)
        a.set_title(H.SYM[q], fontsize=9)
        a.set_xlabel("ω [rad/s]", fontsize=8)
        a.grid(alpha=0.3)
        a.legend(fontsize=5.5)
    fig.suptitle("Heave-plate cylinder, D3: HAMS-MREL on mesh m1 with three lids, 2.0-4.0 rad/s by 0.05, against %s (their "
                 "3464-panel mesh), Nemoh v3.0.3 on m1 and HAMS-MREL on m3" % H.PAPER, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "d3_lid.png"), dpi=95)
    with open(os.path.join(out, "d3_lid.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
