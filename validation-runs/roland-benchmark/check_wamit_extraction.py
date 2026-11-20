"""Check of the extraction: the WAMIT curves of Figs. 5 and 7 against the WAMIT .12d files they were plotted from.
Usage: python -I check_wamit_extraction.py CURVES_CSV FIGURE WAMIT_12D [RHO]
The .12d values are F / (rho g A^2) with ULEN = 1 (WAMIT manual); rho = 1025 and g = 9.80665 in both WAMIT runs."""
import csv
import math
import sys

def read12(path):
    out = {}
    for line in open(path).read().splitlines()[1:]:
        p = line.split()
        if len(p) < 9 or int(p[4]) != 1:
            continue
        wi, wj = round(2 * math.pi / float(p[0]), 2), round(2 * math.pi / float(p[1]), 2)
        out[(wi, wj)] = complex(float(p[7]), float(p[8]))
    return out

curves, fig, path = sys.argv[1], sys.argv[2], sys.argv[3]
rho = float(sys.argv[4]) if len(sys.argv) > 4 else 1025.0
q = read12(path)
scale = rho * 9.80665
worst = 0.0
for dw in ("0.00", "0.05"):
    rows = [r for r in csv.DictReader(open(curves)) if r["figure"] == fig and r["solver"] == "WAMIT" and r["diagonal_dw"] == dw]
    errs = []
    for r in rows:
        w = float(r["omega1"]); w2 = round(w - float(dw), 2)
        v = q.get((w, w2), q.get((w2, w)))
        if v is None:
            continue
        ref = abs(v) * scale
        errs.append((abs(float(r["value_N_per_m2"]) - ref), ref, w))
    peak = max(e[1] for e in errs)
    e = max(errs)
    worst = max(worst, e[0] / peak)
    print("fig %s dw %s: %d points, largest |figure - file| = %.1f N/m^2 at %.2f rad/s, %.3f %% of the largest value %.4g"
          % (fig, dw, len(errs), e[0], e[2], 100 * e[0] / peak, peak))
print("worst %.3f %%" % (100 * worst))
