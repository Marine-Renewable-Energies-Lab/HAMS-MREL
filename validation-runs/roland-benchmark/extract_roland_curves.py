"""Read the curves of Figs. 5 and 7 of Roland et al. (2026) from the vector paths of the PDF.

Source: S. Roland, A. Pribadi, A. Al-Ghuwaidi, E. Lataire (2026). A benchmark of boundary element method solvers for
predicting second-order hydrodynamic loads on floating offshore wind turbines. J. Phys.: Conf. Ser. 3346, 012006,
doi:10.1088/1742-6596/3346/1/012006. Open access, CC BY 4.0: the data written here are derived from that work.

Fig. 5 (DeepCwind, page 8) and Fig. 7 (VolturnUS-S, page 10) plot |QTF| of the surge difference-frequency force in N/m^2
against omega_1, on the mean-drift diagonal (left) and on the diagonal omega_2 - omega_1 = 0.05 rad/s (right), for six
solver set-ups. The figures are matplotlib vector graphics, so each curve is a polyline whose vertices are the plotted
points: no tracing by eye. The axes are calibrated from the tick marks (x: 0.50 to 2.50 in steps of 0.25; y: 0 upwards in
the step of the figure). Each vertex x is snapped to the 0.05 rad/s grid of the paper (section 3.3); a snap error above
0.005 rad/s stops the script, except for Nemoh of Fig. 5, which is on its own grid (Al-Ghuwaidi and Pribadi 2025) and
is kept as read. A point on the 0.05 rad/s diagonal belongs to the pair (omega_1 - 0.05, omega_1) (check_wamit_extraction.py:
the WAMIT curves then match the WAMIT files to 0.04 %).

Usage: python -I extract_roland_curves.py ROLAND_PDF OUT_CSV
The PDF is the user's copy (kept out of the repository). OUT_CSV columns: figure, platform, diagonal_dw, solver, omega1,
value_N_per_m2.
"""
import csv
import sys

import fitz  # PyMuPDF

# Legend order of both figures (matplotlib tab10 colours, in the order of the legend entries)
SOLVERS = {(0.122, 0.467, 0.706): "WAMIT", (1.0, 0.498, 0.055): "Nemoh", (0.173, 0.627, 0.173): "OrcaWave PI",
           (0.839, 0.153, 0.157): "OrcaWave CS", (0.58, 0.404, 0.741): "HydroStar PI", (0.549, 0.337, 0.294): "HydroStar CS"}
# page (1-based), y band of the plots on the page (pt), y tick step (N/m^2), platform
FIGURES = {"5": (8, 590.0, 720.0, 0.5e5, "DeepCwind"), "7": (10, 395.0, 530.0, 0.25e5, "VolturnUS-S")}


def ticks(page, y0, y1):
    """Short tick marks: (x, y) of the x-axis ticks (vertical) and of the y-axis ticks (horizontal) in the band."""
    xt, yt = [], []
    for d in page.get_drawings():
        r = d["rect"]
        if r.y0 < y0 or r.y1 > y1 or len(d["items"]) != 1 or d["items"][0][0] != "l":
            continue
        a, b = d["items"][0][1], d["items"][0][2]
        if abs(a - b) > 3.0:
            continue
        if abs(a.x - b.x) < 1e-3:
            xt.append((round(a.x, 3), round(min(a.y, b.y), 3)))
        else:
            yt.append((round(max(a.x, b.x), 3), round(a.y, 3)))
    return sorted(set(xt)), sorted(set(yt))


def panels(xt, yt):
    """Split the ticks into the left and the right plot by x; return per plot the linear maps x -> omega, y -> value."""
    xs = sorted(x for x, _ in xt)
    split = max(range(1, len(xs)), key=lambda i: xs[i] - xs[i - 1])
    out = []
    for group in (xs[:split], xs[split:]):
        left, right = group[0], group[-1]
        ys = sorted(y for x, y in yt if left - 30.0 < x < left + 1.0)
        out.append(dict(x0=left, x1=right, w0=0.50, w1=0.50 + 0.25 * (len(group) - 1), ytop=ys[0], ybottom=ys[-1],
                        nyt=len(ys)))
    return out


def extract(pdf, out_csv):
    doc = fitz.open(pdf)
    rows = []
    for fig, (pno, y0, y1, ystep, platform) in FIGURES.items():
        page = doc[pno - 1]
        xt, yt = ticks(page, y0, y1)
        maps = panels(xt, yt)
        curves = []
        for d in page.get_drawings():
            r = d["rect"]
            colour = d.get("color")
            if colour is None or r.y0 < y0 or r.y1 > y1 or len(d["items"]) < 8:
                continue
            solver = SOLVERS.get(tuple(round(c, 3) for c in colour))
            if solver is None:
                continue
            side = 0 if r.x1 < 0.5 * (maps[0]["x1"] + maps[1]["x0"]) else 1
            m = maps[side]
            points = [d["items"][0][1]] + [it[2] for it in d["items"] if it[0] == "l"]
            omega = [m["w0"] + (p.x - m["x0"]) / (m["x1"] - m["x0"]) * (m["w1"] - m["w0"]) for p in points]
            curves.append((side, solver, omega, [p.y for p in points]))
        for side in (0, 1):
            m = maps[side]
            vtop = ystep * (m["nyt"] - 1)
            # The visible x range: from the left spine (where the y ticks are) to the right spine, which is the largest
            # vertex x of the plot (matplotlib clips each line there). Vertices outside are not in the figure; a vertex
            # on the right spine that is off the grid is the clip point, not a data point.
            spine = max(x for x, y in yt if m["x0"] - 30.0 < x < m["x0"] + 1.0)
            left = m["w0"] + (spine - m["x0"]) / (m["x1"] - m["x0"]) * (m["w1"] - m["w0"])
            right = max(w for s, _, om, _ in curves if s == side for w in om)
            for s, solver, omega, ys in curves:
                if s != side:
                    continue
                for w, y in zip(omega, ys):
                    ws = round(w / 0.05) * 0.05
                    if w < left - 1e-3 or (w > right - 1e-3 and abs(w - ws) > 1e-3):
                        continue
                    if abs(w - ws) <= 0.005:
                        text = "%.2f" % ws
                    elif solver == "Nemoh" and fig == "5":
                        text = "%.4f" % w   # the DeepCwind Nemoh run of Al-Ghuwaidi and Pribadi has its own grid
                    else:
                        raise ValueError("fig %s %s: omega %.4f is off the 0.05 grid" % (fig, solver, w))
                    v = (m["ybottom"] - y) / (m["ybottom"] - m["ytop"]) * vtop
                    rows.append((fig, platform, "0.00" if side == 0 else "0.05", solver, text, "%.6g" % v))
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("figure", "platform", "diagonal_dw", "solver", "omega1", "value_N_per_m2"))
        w.writerows(rows)
    print("wrote %d points to %s" % (len(rows), out_csv))


if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2])
