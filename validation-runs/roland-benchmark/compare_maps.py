"""Maps of the surge difference-frequency QTF in the style of Roland et al. (2026), Figs. 4 and 6: magnitude (jet, 0 to
6e5 N/m^2) and phase, for HAMS and for the WAMIT file the paper plots, on the HAMS grid, and the modulus of their complex
difference.

Source of the style and the colour range: Roland, S., Pribadi, A., Al-Ghuwaidi, A., Lataire, E. (2026). J. Phys.: Conf.
Ser. 3346, 012006, doi:10.1088/1742-6596/3346/1/012006, CC BY 4.0. The WAMIT files are those of the paper (OpenFAST r-test
marin_semi.12d for DeepCwind; IEA-15-240-RWT v1.1.18 IEA-15-240-RWT-UMaineSemi.12d for VolturnUS-S); they stay outside the
repository.

Both files are in the WAMIT format and are read as written (one convention for both). Values are |T| rho g with rho = 1025
kg/m^3, g = 9.80665 m/s^2 (see compare_diagonals.py). Phase in radians, of the (omega_1, omega_2) entry; the lower triangle
is the conjugate of the upper one.

Measure: |HAMS - WAMIT| / max |WAMIT| over the pairs of the HAMS grid with both frequencies in [LO, HI], for all pairs, for
0 < omega_2 - omega_1 <= 0.30 rad/s, and without the pairs that have a frequency in an optional band.

Usage: python -I compare_maps.py HAMS_FILE WAMIT_12D OUT_PNG OUT_TXT TITLE [--range LO HI] [--band A B] [--mode I]
"""
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RHO_G = 1025.0 * 9.80665
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def read(path, mode, periods):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split()
        if len(p) != 9:
            continue
        try:
            vals = [float(v) for v in p]
        except ValueError:
            continue
        if vals[2] != 0.0 or vals[3] != 0.0 or int(vals[4]) != mode:
            continue
        a, b = (2 * math.pi / vals[0], 2 * math.pi / vals[1]) if periods else (vals[0], vals[1])
        a, b = round(a, 2), round(b, 2)
        v = complex(vals[7], vals[8]) * RHO_G
        out[(a, b)] = v
        out.setdefault((b, a), v.conjugate())
    return out


def main(argv):
    opts = dict(range=(0.30, 1.90), band=None, mode=1)
    for key, n in (("--range", 2), ("--band", 2), ("--mode", 1)):
        if key in argv:
            k = argv.index(key)
            vals = argv[k + 1:k + 1 + n]
            opts[key[2:]] = tuple(float(v) for v in vals) if n == 2 else int(vals[0])
            del argv[k:k + 1 + n]
    hams_path, wamit_path, out_png, out_txt, title = argv[1:6]
    lo, hi = opts["range"]
    h = read(hams_path, opts["mode"], periods=False)
    w = read(wamit_path, opts["mode"], periods=True)
    grid = sorted({a for a, _ in h if lo - 1e-9 <= a <= hi + 1e-9})
    n = len(grid)
    H = np.full((n, n), np.nan, complex)
    W = np.full((n, n), np.nan, complex)
    for i, a in enumerate(grid):
        for j, b in enumerate(grid):
            if (a, b) in h:
                H[j, i] = h[(a, b)]
            if (a, b) in w:
                W[j, i] = w[(a, b)]
    D = np.abs(H - W)
    scale = np.nanmax(np.abs(W))
    lines = [title, "|HAMS - WAMIT| / max |WAMIT| (max |WAMIT| = %.4g N/m^2), omega in [%.2f, %.2f], largest / median, %%"
             % (scale, lo, hi)]
    g = np.array(grid)
    dw = g[None, :] - g[:, None]   # [j, i] = omega_i - omega_j
    upper = np.abs(dw) > -1   # all
    sets = [("all pairs", np.ones_like(D, bool)), ("0 < dw <= 0.30", (np.abs(dw) > 1e-9) & (np.abs(dw) <= 0.30 + 1e-9)),
            ("diagonal", np.abs(dw) < 1e-9)]
    for name, m in sets:
        for band in ([None] if opts["band"] is None else [None, opts["band"]]):
            mm = m & np.isfinite(D)
            if band is not None:
                inb = (g >= band[0] - 1e-9) & (g <= band[1] + 1e-9)
                mm &= ~(inb[None, :] | inb[:, None])
            e = 100 * D[mm] / scale
            lines.append("  %-15s %-20s %6.2f / %5.2f  (%d pairs)" % (
                name, "" if band is None else "w/o %.2f-%.2f rad/s" % band, e.max(), np.median(e), mm.sum()))

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "figure.facecolor": SURFACE})
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.2), layout="constrained")
    ext = (lo - 0.025, hi + 0.025, lo - 0.025, hi + 0.025)
    for col, (name, M) in enumerate((("HAMS", H), ("WAMIT file", W))):
        im = axes[0, col].imshow(np.abs(M), origin="lower", extent=ext, cmap="jet", vmin=0, vmax=6e5)
        axes[0, col].set_title(name + ", magnitude", loc="left")
        fig.colorbar(im, ax=axes[0, col], label="|QTF| (N/m$^2$)", extend="max", shrink=0.85)
        ip = axes[1, col].imshow(np.angle(M), origin="lower", extent=ext, cmap="jet", vmin=-math.pi, vmax=math.pi)
        axes[1, col].set_title(name + ", phase", loc="left")
        fig.colorbar(ip, ax=axes[1, col], label="Phase (rad)", shrink=0.85)
    idiff = axes[0, 2].imshow(100 * D / scale, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=10)
    axes[0, 2].set_title("Difference", loc="left")
    fig.colorbar(idiff, ax=axes[0, 2], label="|HAMS $-$ WAMIT| / max|WAMIT| (%)", extend="max", shrink=0.85)
    axes[1, 2].axis("off")
    table = ["|HAMS - WAMIT| / max|WAMIT|", "largest / median, %", "max|WAMIT| = %.3g N/m^2" % scale,
             "omega in [%.2f, %.2f] rad/s" % (lo, hi), ""] + [l.strip() for l in lines[2:]]
    axes[1, 2].text(0.0, 1.0, "\n".join(table), va="top", ha="left", fontsize=7.5,
                    family="monospace", color=INK, transform=axes[1, 2].transAxes)
    for ax in axes.flat[:5]:
        ax.set_xlabel(r"$\omega_1$ (rad/s)")
        ax.set_ylabel(r"$\omega_2$ (rad/s)")
    fig.suptitle(title + "\nStyle and colour range of Roland et al. (2026), J. Phys.: Conf. Ser. 3346, 012006, Figs. 4 "
                 "and 6 (CC BY 4.0)", x=0.01, ha="left", color=INK, fontsize=10)
    fig.savefig(out_png, dpi=150)
    open(out_txt, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(list(sys.argv))
