"""Read the curves of Fig. 8 of Kurnia and Ducrozet (2023) from the vector paths of the PDF.

Source: R. Kurnia, G. Ducrozet (2023). NEMOH: Open-source boundary element solver for computation of first- and
second-order hydrodynamic loads in the frequency domain. Comput. Phys. Commun. 292, 108885,
doi:10.1016/j.cpc.2023.108885. The PDF is the user's copy and stays out of the repository; so does the CSV written here.

Fig. 8 (page 10) plots, for the truncated cylinder (radius 10 m, draft 10 m) with a thin heave plate (radius 15 m,
0.1 m thick), deep water, CG (0, 0, -7.5) m: the added mass A_11, A_33, A_55 (left column), the damping B_11, B_33,
B_55 (right column) and the excitation |X_1|, |X_3|, |X_5| per unit amplitude (bottom), against omega in rad/s, for
WAMIT (black '+', no irregular-frequency removal), HydroStar (blue solid, with removal), DIODORE (cyan dotted, with
removal), Nemoh v2 (red dashed, no removal) and Nemoh v3 (green dash-dotted, with removal). The figure is MATLAB vector
graphics: each solver curve is a polyline whose vertices are the plotted points, and each WAMIT point is a '+' of two
strokes, whose crossing is the point. No tracing by eye.

Calibration: the x axis of every plot runs from 0 to 4 rad/s between the left and right spines (the tick marks at 0 and 4
sit on the spines; checked below). The y axis is calibrated from the tick marks on the left spine. The tick labels are
glyph outlines, not text, so their values are listed in Y_TICKS, read from the page rendered at 200 dpi; the script
stops if the number of tick marks on a spine differs from the number of listed values.

Usage: python -I extract_kurnia_fig8.py KURNIA_PDF OUT_CSV
OUT_CSV columns: quantity, solver, omega_rad_s, value (SI units: kg, kg m^2, N s/m, N m s, N/m, N).
"""
import csv
import sys

import fitz  # PyMuPDF

PAGE = 10
# Axes boxes of Fig. 8 (white-filled rectangles), identified by their top-left corner in points.
# quantity, (x0, y0) of the box, tick values on the left spine from bottom to top, scale of the axis
AXES = [
    ("A11", (134.76, 327.66), (0, 2, 4), 1e6),
    ("A33", (134.76, 392.02), (6.0, 6.5), 1e6),
    ("A55", (134.76, 456.23), (2.2, 2.3, 2.4), 1e8),
    ("B11", (333.96, 327.66), (0, 1, 2), 1e6),
    ("B33", (333.96, 392.02), (0, 10, 20), 1e4),
    ("B55", (333.96, 456.23), (0, 10, 20), 1e6),
    ("X1", (232.36, 527.40), (0, 2, 4), 1e6),
    ("X3", (232.36, 591.76), (0, 1, 2), 1e6),
    ("X5", (232.36, 655.97), (0, 2), 1e6),
]
SOLVERS = {(0.0, 0.0, 1.0): "HydroStar", (0.0, 1.0, 1.0): "DIODORE", (1.0, 0.0, 0.0): "Nemoh v2",
           (0.0, 1.0, 0.0): "Nemoh v3"}
TOL = 0.6  # pt


def find_boxes(drawings):
    boxes = {}
    for d in drawings:
        if d.get("fill") != (1.0, 1.0, 1.0) or len(d["items"]) != 1 or d["items"][0][0] != "re":
            continue
        r = d["rect"]
        for name, (x0, y0), _, _ in AXES:
            if abs(r.x0 - x0) < 0.05 and abs(r.y0 - y0) < 0.05 and r.width > 100:
                boxes[name] = r
    missing = [a[0] for a in AXES if a[0] not in boxes]
    if missing:
        raise ValueError("axes boxes not found: %s" % missing)
    return boxes


def black_strokes(drawings):
    out = []
    for d in drawings:
        if d.get("color") != (0.0, 0.0, 0.0) or len(d["items"]) != 1 or d["items"][0][0] != "l":
            continue
        a, b = d["items"][0][1], d["items"][0][2]
        out.append((a, b))
    return out


def calibrate(box, strokes, ticks):
    """Linear map y(pt) -> value from the short horizontal tick strokes starting on the left spine."""
    ys = sorted({round(a.y, 2) for a, b in strokes
                 if abs(a.y - b.y) < 1e-3 and min(abs(a.x - box.x0), abs(b.x - box.x0)) < 0.05
                 and abs(a.x - b.x) < 4.0 and box.y0 - TOL <= a.y <= box.y1 + TOL})
    if len(ys) != len(ticks):
        raise ValueError("box at %s: %d tick marks on the left spine, %d values listed" % (box, len(ys), len(ticks)))
    ys = ys[::-1]  # bottom first
    slope = (ticks[-1] - ticks[0]) / (ys[-1] - ys[0])
    return lambda y: ticks[0] + (y - ys[0]) * slope, ys


def xmap(box):
    return lambda x: 4.0 * (x - box.x0) / (box.x1 - box.x0)


def plus_centres(box, strokes):
    """Centres of the '+' markers inside the box: a horizontal and a vertical stroke of equal length crossing at
    their midpoints."""
    hor, ver = {}, []
    for a, b in strokes:
        mx, my = 0.5 * (a.x + b.x), 0.5 * (a.y + b.y)
        if not (box.x0 + TOL < mx < box.x1 - TOL and box.y0 - 2 < my < box.y1 + 2):
            continue
        length = abs(a.x - b.x) + abs(a.y - b.y)
        if length > 4.0:
            continue
        if abs(a.y - b.y) < 1e-3:
            hor[(round(mx, 2), round(my, 2))] = length
        elif abs(a.x - b.x) < 1e-3:
            ver.append((round(mx, 2), round(my, 2), length))
    out = []
    for mx, my, length in ver:
        if (mx, my) in hor and abs(hor[(mx, my)] - length) < 0.05:
            out.append((mx, my))
    return out


def extract(pdf, out_csv):
    doc = fitz.open(pdf)
    page = doc[PAGE - 1]
    drawings = page.get_drawings()
    boxes = find_boxes(drawings)
    strokes = black_strokes(drawings)
    # legend boxes: white-filled rectangles inside an axes box, smaller than it
    legends = [d["rect"] for d in drawings
               if d.get("fill") == (1.0, 1.0, 1.0) and len(d["items"]) == 1 and d["items"][0][0] == "re"
               and d["rect"].width < 100 and d["rect"].y0 > 320]
    rows = []
    for name, _, ticks, scale in AXES:
        box = boxes[name]
        fy, _ = calibrate(box, strokes, ticks)
        fx = xmap(box)
        n_curves = {}
        for d in drawings:
            col = d.get("color")
            if col is None or len(d["items"]) < 8:
                continue
            solver = SOLVERS.get(tuple(round(c, 3) for c in col))
            r = d["rect"]
            if solver is None or not (box.x0 - TOL <= r.x0 and r.x1 <= box.x1 + TOL and box.y0 - 30 <= r.y0
                                      and r.y1 <= box.y1 + 30):
                continue
            if any(lg.contains(r) for lg in legends):
                continue
            n_curves[solver] = n_curves.get(solver, 0) + 1
            points = [d["items"][0][1]] + [it[2] for it in d["items"] if it[0] == "l"]
            for p in points:
                rows.append((name, solver, "%.5f" % fx(p.x), "%.6g" % (fy(p.y) * scale)))
        if sorted(n_curves) != sorted(SOLVERS.values()) or any(v != 1 for v in n_curves.values()):
            raise ValueError("%s: curves found %s, expected one per solver" % (name, n_curves))
        pts = [c for c in plus_centres(box, strokes) if not any(lg.contains(fitz.Point(*c)) for lg in legends)]
        if len(pts) < 10:
            raise ValueError("%s: only %d WAMIT points" % (name, len(pts)))
        for x, y in sorted(pts):
            rows.append((name, "WAMIT", "%.5f" % fx(x), "%.6g" % (fy(y) * scale)))
        print("%s: %d WAMIT points, curves %s" % (name, len(pts), n_curves))
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("quantity", "solver", "omega_rad_s", "value"))
        w.writerows(rows)
    print("wrote %d points to %s" % (len(rows), out_csv))


if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2])
