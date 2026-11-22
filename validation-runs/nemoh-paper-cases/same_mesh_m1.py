"""Same-mesh check of the heave-plate mesh m1: Nemoh v3.0.3 against HAMS-MREL on m1 (hull 3096 + lid 744), with the
paper's WAMIT and Nemoh v3 (Kurnia and Ducrozet 2023, Fig. 8, their 3464-panel mesh) and HAMS-MREL on m3 for scale.

Usage: python -I same_mesh_m1.py FIG8_CSV RESULTS_DIR OUT_DIR
Writes same_mesh_m1.png and same_mesh_m1.txt: per quantity, |Nemoh m1 - HAMS m1| / max|HAMS m1| (max, median, %) and
the same for each code on m1 against HAMS m3, and the frequencies where a damping is negative.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_heave_plate as H  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_fig8(path):
    import csv
    ref = {}
    for r in csv.DictReader(open(path)):
        ref.setdefault((r["quantity"], r["solver"]), []).append((float(r["omega_rad_s"]), float(r["value"])))
    return {k: (np.array([x for x, _ in sorted(v)]), np.array([y for _, y in sorted(v)])) for k, v in ref.items()}


def main(csv_path, res, out):
    ref = read_fig8(csv_path)
    hm1, hm3 = H.hams(os.path.join(res, "hp_m1")), H.hams(os.path.join(res, "hp_m3"))
    nm1 = H.nemoh(os.path.join(res, "nemoh_hp_m1"))
    lines = ["quantity  |Nemoh m1 - HAMS m1|  |HAMS m1 - HAMS m3|  |Nemoh m1 - HAMS m3|   (max / median, % of max|.|)"]
    fig, ax = plt.subplots(3, 3, figsize=(15, 11))
    for k, q in enumerate(H.QUANT):
        a = ax[k // 3][k % 3]
        w1, v1 = hm1[q]
        wn, vn = nm1[q]
        w3, v3 = hm3[q]
        vn_on1 = np.interp(w1, wn, vn)
        sc = np.max(abs(v1))
        d1 = abs(vn_on1 - v1) / sc * 100
        d2 = abs(v1 - np.interp(w1, w3, v3)) / np.max(abs(v3)) * 100
        d3 = abs(vn_on1 - np.interp(w1, w3, v3)) / np.max(abs(v3)) * 100
        neg = ""
        if q.startswith("B"):
            neg = " negative Nemoh B at %s" % (", ".join("%.1f" % x for x in wn[vn < 0]) or "none")
            neg += "; HAMS m1 at %s" % (", ".join("%.1f" % x for x in w1[v1 < 0]) or "none")
        lines.append("%-6s %8.2f / %5.2f       %8.2f / %5.2f       %8.2f / %5.2f   %s"
                     % (q, d1.max(), np.median(d1), d2.max(), np.median(d2), d3.max(), np.median(d3), neg))
        if ref is not None:
            for src in ("WAMIT", "Nemoh v3"):
                if (q, src) in ref:
                    wr, vr = ref[(q, src)]
                    a.plot(wr, vr, label="%s (%s)" % (src, H.PAPER), **H.STYLE[src])
        a.plot(w3, v3, color="tab:green", lw=1.6, label="HAMS-MREL, mesh m3 (11400 + 1240 lid)")
        a.plot(w1, v1, color="tab:red", lw=1.6, label="HAMS-MREL, mesh m1 (3096 + 744 lid)")
        a.plot(wn, vn, color="tab:blue", lw=1.6, ls="-", marker="o", ms=3, label="Nemoh v3.0.3 (our run), mesh m1")
        a.set_title(H.SYM[q], fontsize=9)
        a.set_xlabel("ω [rad/s]", fontsize=8)
        a.grid(alpha=0.3)
        a.legend(fontsize=6)
    fig.suptitle("Heave-plate cylinder: Nemoh v3.0.3 and HAMS-MREL on the same mesh m1, against the paper's WAMIT and "
                 "Nemoh v3 (their 3464-panel mesh) and HAMS-MREL on m3", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "same_mesh_m1.png"), dpi=95)
    with open(os.path.join(out, "same_mesh_m1.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
