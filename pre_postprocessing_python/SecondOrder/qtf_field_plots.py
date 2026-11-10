"""
Maps and field-point grids for QTF_FieldPartial.txt of HAMS-MREL (the partial second-order field).

Usage:
    python -I qtf_field_plots.py plot <QTF_FieldPartial.txt> --out <folder> [--kernel diff|sum|both]
                                      [--type pressure|elevation|both]
    python -I qtf_field_plots.py grid --x <x0> <x1> <nx> --y <y0> <y1> <ny> [--z 0] [--exclude <xc> <yc> <r>]...
    python -I qtf_field_plots.py check

plot: for every kernel, frequency pair, heading pair, type and depth level (pressure points are grouped by z), one PNG
with two panels, the modulus and the phase (degrees, cyclic colour map, at the points) over x, y. The values are those of the file:
pressure / (rho g A^2 / L) and elevation / (A^2 / L), e^{+i omega t} as the other WAMIT-format files. The field is
PARTIAL: the second-order diffracted and radiated waves are omitted, so the maps are for the far field and for checks.
grid: prints the field-point block of ControlFile.in (Number_of_field_points and one "x y z" line per point) for a
rectangular grid, leaving out the points inside discs around the bodies.
check: a synthetic file is written, read back, plotted and gridded; exits 1 on any failure.
"""
import argparse
import os
import sys
import tempfile

# python -I does not put the script's folder on sys.path; do it here so that sibling modules import.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLUMNS = ["per1", "per2", "beta1", "beta2", "k", "ipt", "typ", "mod", "pha", "re", "im",
           "re_quad", "im_quad", "re_inc", "im_inc"]
INTEGER = ("k", "ipt", "typ")
DTYPE = np.dtype([(c, "i8" if c in INTEGER else "f8") for c in COLUMNS])
TITLE = "PARTIAL second-order field (no second-order diffracted or radiated waves)"
KERNELS = {1: "diff", 2: "sum"}
TYPES = {1: "pressure", 2: "elevation"}
KERNEL_LABEL = {"diff": "difference frequency T-", "sum": "sum frequency T+"}
QUANTITY = {"pressure": "pressure / (rho g A^2 / L)", "elevation": "elevation / (A^2 / L)"}


def read_field_partial(path):
    """(points, rows): points {ipt: (x, y, z, type)}; rows a structured array of the 15 columns of the file."""
    points, rows = {}, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                w = s[1:].split()
                if len(w) == 6 and w[0] == "POINT" and w[1] != "IPT":
                    points[int(w[1])] = (float(w[2]), float(w[3]), float(w[4]), int(w[5]))
                continue
            w = s.split()
            if len(w) != 15:
                raise ValueError("%s: expected 15 columns, got %d: %s" % (path, len(w), s))
            rows.append(tuple(int(v) if COLUMNS[i] in INTEGER else float(v) for i, v in enumerate(w)))
    return points, np.array(rows, dtype=DTYPE)


def _num(v):
    return ("%g" % v).replace("-", "m")


def png_name(kernel, per1, per2, beta1, beta2, kind, z=None):
    name = "field_%s_%.4f_%.4f_%s_%s_%s" % (kernel, per1, per2, _num(beta1), _num(beta2), kind)
    return name + ("_z%s" % _num(z) if kind == "pressure" else "") + ".png"


def _collinear(x, y):
    if len(x) < 3:
        return True
    s = np.linalg.svd(np.column_stack([x - x.mean(), y - y.mean()]), compute_uv=False)
    return s[0] == 0.0 or s[1] <= 1e-9 * s[0]


def _markers_only(x, y, v):
    """True when tricontourf cannot be used: fewer than 3 non-collinear points, duplicate (x, y), or constant values."""
    if _collinear(x, y) or np.ptp(v) == 0.0:
        return True
    return len({(round(a, 9), round(b, 9)) for a, b in zip(x, y)}) < len(x)


def _panel(fig, ax, x, y, v, label, cmap, vmin=None, vmax=None, points_only=False):
    if points_only or _markers_only(x, y, v):
        sc = ax.scatter(x, y, c=v, cmap=cmap, vmin=vmin, vmax=vmax, s=90, edgecolors="k")
    else:
        sc = ax.tricontourf(x, y, v, levels=20, cmap=cmap)
        ax.plot(x, y, "k.", ms=3)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.colorbar(sc, ax=ax, label=label)


def plot(path, out, kernel="both", kind="both"):
    """Write the PNGs; returns the list of file names."""
    points, rows = read_field_partial(path)
    os.makedirs(out, exist_ok=True)
    written = []
    keys = sorted({(r["k"], r["per1"], r["per2"], r["beta1"], r["beta2"], r["typ"]) for r in rows})
    for k, per1, per2, beta1, beta2, typ in keys:
        kern, kd = KERNELS[k], TYPES[typ]
        if kernel not in ("both", kern) or kind not in ("both", kd):
            continue
        sel = rows[(rows["k"] == k) & (rows["per1"] == per1) & (rows["per2"] == per2) & (rows["beta1"] == beta1)
                   & (rows["beta2"] == beta2) & (rows["typ"] == typ)]
        levels = {}
        for r in sel:
            x, y, z, _ = points[int(r["ipt"])]
            levels.setdefault(round(z, 9) if typ == 1 else 0.0, []).append((x, y, r["mod"], r["pha"]))
        for z, pts in sorted(levels.items()):
            x, y, mod, pha = (np.array(c) for c in zip(*pts))
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            _panel(fig, axes[0], x, y, mod, "modulus, " + QUANTITY[kd], "viridis")
            _panel(fig, axes[1], x, y, pha, "phase (deg)", "twilight", -180.0, 180.0, points_only=True)
            axes[0].set_title("modulus")
            axes[1].set_title("phase at the points (not interpolated)")
            fig.suptitle("%s\n%s, %s%s, T1 = %.4g s, T2 = %.4g s, beta1 = %g, beta2 = %g deg"
                         % (TITLE, KERNEL_LABEL[kern], kd, " at z = %g m" % z if typ == 1 else "",
                            per1, per2, beta1, beta2), fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 0.9))
            name = png_name(kern, per1, per2, beta1, beta2, kd, z)
            fig.savefig(os.path.join(out, name), dpi=110)
            plt.close(fig)
            written.append(name)
    return written


def grid(x, y, z=0.0, exclude=()):
    """Lines of the ControlFile.in field-point block for the nx * ny grid, minus the discs (xc, yc, r)."""
    xs = np.linspace(x[0], x[1], int(x[2]))
    ys = np.linspace(y[0], y[1], int(y[2]))
    pts = [(a, b) for a in xs for b in ys
           if not any((a - xc) ** 2 + (b - yc) ** 2 <= r * r for xc, yc, r in exclude)]
    return ["Number_of_field_points %d" % len(pts)] + ["%.6g %.6g %.6g" % (a, b, z) for a, b in pts]


# ------------------------------------------------------------------------------------------------------------- check

def _write_synthetic(path):
    """Two pressure levels (z = -2, -6) on a 4 x 4 grid and an elevation line; one pair per kernel."""
    pts = {}
    ipt = 0
    for z in (-2.0, -6.0):
        for a in range(4):
            for b in range(4):
                ipt += 1
                pts[ipt] = (10.0 * a, 10.0 * b, z, 1)
    for a in range(4):
        ipt += 1
        pts[ipt] = (10.0 * a, 0.0, 0.0, 2)
    with open(path, "w") as f:
        f.write("# QTF_FieldPartial.txt: synthetic\n# POINT IPT X Y Z TYPE\n")
        for i, (x, y, z, t) in pts.items():
            f.write("# POINT %5d %13.6E %13.6E %13.6E %5d\n" % (i, x, y, z, t))
        for k in (1, 2):
            for i, (x, y, z, t) in pts.items():
                v = complex(0.01 * i, -0.002 * i * k)
                q, c = 0.75 * v, 0.25 * v
                vals = (abs(v), float(np.degrees(np.angle(v))), v.real, v.imag, q.real, q.imag, c.real, c.imag)
                f.write("%14.6E%14.6E%14.6E%14.6E%6d%6d%6d%s\n" % (
                    8.0, 8.0 + k, 0.0, 30.0, k, i, t, "".join("%14.6E" % u for u in vals)))
    return pts


def check():
    failed = 0

    def test(name, ok, info=""):
        nonlocal failed
        failed += not ok
        print("%-28s %s %s" % (name, "ok" if ok else "FAILED", info))

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "QTF_FieldPartial.txt")
        pts = _write_synthetic(path)
        points, rows = read_field_partial(path)
        test("read points", points == pts and len(points) == 36)
        test("read rows", len(rows) == 72 and rows.dtype.names == tuple(COLUMNS))
        r = rows[(rows["k"] == 2) & (rows["ipt"] == 5)][0]
        v = complex(0.05, -0.002 * 5 * 2)
        test("read values exact", abs(r["re"] - v.real) < 1e-12 and abs(r["im"] - v.imag) < 1e-12
             and abs(r["mod"] - abs(v)) < 1e-6 * abs(v) and r["typ"] == 1 and r["per2"] == 10.0)
        out = os.path.join(tmp, "png")
        names = sorted(plot(path, out))
        want = sorted(["field_%s_8.0000_%.4f_0_30_%s.png" % (k, 8.0 + i, s)
                       for i, k in ((1, "diff"), (2, "sum")) for s in ("pressure_zm2", "pressure_zm6", "elevation")])
        test("png names", names == want, "" if names == want else "%s != %s" % (names, want))
        test("png files exist", all(os.path.getsize(os.path.join(out, n)) > 0 for n in names))
        test("kernel/type filter", plot(path, os.path.join(tmp, "f"), "diff", "elevation")
             == ["field_diff_8.0000_9.0000_0_30_elevation.png"])
        g = grid((0, 30, 4), (0, 30, 4))
        test("duplicates -> markers", _markers_only(np.array([0., 1, 0, 1, 0]), np.array([0., 0, 1, 1, 0]),
                                                    np.arange(5.)) and not _markers_only(
            np.array([0., 1, 0, 1]), np.array([0., 0, 1, 1]), np.arange(4.)))
        test("grid count", g[0] == "Number_of_field_points 16" and len(g) == 17)
        g = grid((0, 30, 4), (0, 30, 4), exclude=[(0.0, 0.0, 12.0)])
        test("grid with exclusion", g[0] == "Number_of_field_points 13" and len(g) == 14, g[0])
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check")
    p = sub.add_parser("plot")
    p.add_argument("path")
    p.add_argument("--out", required=True)
    p.add_argument("--kernel", choices=("diff", "sum", "both"), default="both")
    p.add_argument("--type", choices=("pressure", "elevation", "both"), default="both", dest="kind")
    g = sub.add_parser("grid")
    g.add_argument("--x", nargs=3, type=float, required=True, metavar=("X0", "X1", "NX"))
    g.add_argument("--y", nargs=3, type=float, required=True, metavar=("Y0", "Y1", "NY"))
    g.add_argument("--z", type=float, default=0.0)
    g.add_argument("--exclude", nargs=3, type=float, action="append", default=[], metavar=("XC", "YC", "R"))
    args = parser.parse_args(argv)
    if args.command == "check":
        return check()
    if args.command == "plot":
        print("\n".join(plot(args.path, args.out, args.kernel, args.kind)))
        return 0
    print("\n".join(grid(args.x, args.y, args.z, args.exclude)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
