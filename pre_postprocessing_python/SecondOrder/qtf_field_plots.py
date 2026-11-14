"""
Maps and field-point grids for QTF_FieldPartial.txt and QTF_Field.txt of HAMS-MREL (the second-order field).

Usage:
    python -I qtf_field_plots.py plot <QTF_FieldPartial.txt or QTF_Field.txt> --out <folder> [--kernel diff|sum|both]
                                      [--type pressure|elevation|both]
    python -I qtf_field_plots.py grid --x <x0> <x1> <nx> --y <y0> <y1> <ny> [--z 0] [--exclude <xc> <yc> <r>]...
    python -I qtf_field_plots.py check

plot: for every kernel, frequency pair, heading pair, type and depth level (pressure points are grouped by z), one PNG
with two panels, the modulus and the phase (degrees, cyclic colour map, at the points) over x, y. The values are those of the file:
pressure / (rho g A^2 / L) and elevation / (A^2 / L), e^{+i omega t} as the other WAMIT-format files. The field of
QTF_FieldPartial.txt is PARTIAL: the second-order diffracted and radiated waves are omitted, so the maps are for the far
field and for checks. QTF_Field.txt (qtf_field = 1 with qtf_direct = 1) has the columns Re_dir Im_dir has_direct after those
of the partial file, read by the names of its "# Columns:" line; its PNGs have the modulus of the partial value beside that
of the complete value (partial + direct, where has_direct = 1) and the phase of the complete value.
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
DIRECT_COLUMNS = ["re_dir", "im_dir", "has_direct"]
INTEGER_ALL = INTEGER + ("has_direct",)
TITLE = "PARTIAL second-order field (no second-order diffracted or radiated waves)"
TITLE_COMPLETE = ("second-order field with the direct part (body forcing; no free-surface forcing, no radiated wave of "
                  "the second-order motion)")
KERNELS = {1: "diff", 2: "sum"}
TYPES = {1: "pressure", 2: "elevation"}
KERNEL_LABEL = {"diff": "difference frequency T-", "sum": "sum frequency T+"}
QUANTITY = {"pressure": "pressure / (rho g A^2 / L)", "elevation": "elevation / (A^2 / L)"}


def read_field_partial(path):
    """(points, rows): points {ipt: (x, y, z, type)}; rows a structured array of the columns of the file, named by its
    "# Columns:" line (lower case): the 15 columns of QTF_FieldPartial.txt, and with the direct part re_dir, im_dir and
    has_direct after them. Without such a line, the 15 columns."""
    points, rows = {}, []
    names = list(COLUMNS)
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                w = s[1:].split()
                if len(w) > 2 and w[0] == "Columns:":
                    names = [{"type": "typ"}.get(c.lower(), c.lower()) for c in w[1:]]
                    unknown = [c for c in names if c not in COLUMNS + DIRECT_COLUMNS]
                    if unknown:
                        raise ValueError("%s: unknown column(s) %s" % (path, unknown))
                if len(w) == 6 and w[0] == "POINT" and w[1] != "IPT":
                    points[int(w[1])] = (float(w[2]), float(w[3]), float(w[4]), int(w[5]))
                continue
            w = s.split()
            if len(w) != len(names):
                raise ValueError("%s: expected %d columns, got %d: %s" % (path, len(names), len(w), s))
            rows.append(tuple(int(v) if names[i] in INTEGER_ALL else float(v) for i, v in enumerate(w)))
    dtype = np.dtype([(c, "i8" if c in INTEGER_ALL else "f8") for c in names])
    return points, np.array(rows, dtype=dtype)


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
    complete = "re_dir" in rows.dtype.names
    title = TITLE_COMPLETE if complete else TITLE
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
            # The complete value: the partial one plus the direct part, which a row without it (has_direct = 0) has as 0
            re, im = r["re"], r["im"]
            if complete:
                re, im = re + r["re_dir"], im + r["im_dir"]
            levels.setdefault(round(z, 9) if typ == 1 else 0.0, []).append(
                (x, y, r["mod"], r["pha"], np.hypot(re, im), np.degrees(np.arctan2(im, re))))
        for z, pts in sorted(levels.items()):
            x, y, mod, pha, cmod, cpha = (np.array(c) for c in zip(*pts))
            ncol = 3 if complete else 2
            fig, axes = plt.subplots(1, ncol, figsize=(6 * ncol, 5))
            _panel(fig, axes[0], x, y, mod, "modulus, " + QUANTITY[kd], "viridis")
            axes[0].set_title("modulus, partial" if complete else "modulus")
            if complete:
                _panel(fig, axes[1], x, y, cmod, "modulus, " + QUANTITY[kd], "viridis")
                axes[1].set_title("modulus, complete (with the direct part)")
            _panel(fig, axes[-1], x, y, cpha if complete else pha, "phase (deg)", "twilight", -180.0, 180.0,
                   points_only=True)
            axes[-1].set_title("phase at the points (not interpolated)" + (", complete" if complete else ""))
            fig.suptitle("%s\n%s, %s%s, PER1 = %.4g, PER2 = %.4g, beta1 = %g, beta2 = %g deg"
                         % (title, KERNEL_LABEL[kern], kd, " at z = %g m" % z if typ == 1 else "",
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
    return ["    Number_of_field_points     %d" % len(pts)] + ["%.6g %.6g %.6g" % (a, b, z) for a, b in pts]


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


def _write_synthetic_complete(path):
    """The file of _write_synthetic as QTF_Field.txt: Re_dir Im_dir has_direct after each row (half the partial value for
    the difference rows, k = 1, and 0 for the sum rows), and the "# Columns:" line that names them."""
    _write_synthetic(path)
    with open(path) as f:
        lines = f.readlines()
    header = ("# Columns: PER1 PER2 BETA1 BETA2 K IPT TYPE Mod Pha Re Im Re_quad Im_quad Re_inc Im_inc Re_dir Im_dir "
              "has_direct\n")
    out = [lines[0], header]
    for line in lines[1:]:
        if not line.startswith("#") and line.strip():
            w = line.split()
            has = int(w[4] == "1")
            line = line.rstrip("\n") + "%14.6E%14.6E%6d\n" % (0.5 * float(w[9]) * has, 0.5 * float(w[10]) * has, has)
        out.append(line)
    with open(path, "w") as f:
        f.writelines(out)


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
        # The file with the direct part, read by column name
        cpath = os.path.join(tmp, "QTF_Field.txt")
        _write_synthetic_complete(cpath)
        cpoints, crows = read_field_partial(cpath)
        test("read direct columns", cpoints == pts and crows.dtype.names == tuple(COLUMNS + DIRECT_COLUMNS)
             and len(crows) == 72)
        r = crows[(crows["k"] == 1) & (crows["ipt"] == 5)][0]
        test("direct values", abs(r["re_dir"] - 0.5 * r["re"]) < 1e-12 and r["has_direct"] == 1
             and crows[(crows["k"] == 2)]["has_direct"].max() == 0)
        cnames = sorted(plot(cpath, os.path.join(tmp, "c")))
        test("complete png names", cnames == want, "" if cnames == want else "%s != %s" % (cnames, want))
        test("complete png files exist", all(os.path.getsize(os.path.join(tmp, "c", n)) > 0 for n in cnames))
        g = grid((0, 30, 4), (0, 30, 4))
        test("duplicates -> markers", _markers_only(np.array([0., 1, 0, 1, 0]), np.array([0., 0, 1, 1, 0]),
                                                    np.arange(5.)) and not _markers_only(
            np.array([0., 1, 0, 1]), np.array([0., 0, 1, 1]), np.arange(4.)))
        test("grid count", int(g[0][27:43]) == 16 and len(g) == 17)
        g = grid((0, 30, 4), (0, 30, 4), exclude=[(0.0, 0.0, 12.0)])
        test("grid with exclusion", int(g[0][27:43]) == 13 and len(g) == 14, g[0])
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
