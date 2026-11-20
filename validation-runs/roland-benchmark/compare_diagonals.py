"""HAMS against the six solver set-ups of Roland et al. (2026), Fig. 5 (DeepCwind) or Fig. 7 (VolturnUS-S): |QTF| of the
surge difference-frequency force on the mean-drift diagonal and on the diagonal of 0.05 rad/s.

The paper's curves come from extract_roland_curves.py (vector paths of the PDF; checked against the WAMIT files to 0.04 % by
check_wamit_extraction.py). Source: Roland, S., Pribadi, A., Al-Ghuwaidi, A., Lataire, E. (2026), J. Phys.: Conf. Ser.
3346, 012006, doi:10.1088/1742-6596/3346/1/012006, CC BY 4.0.

The 0.05 diagonal is plotted at the higher frequency of each pair, (omega_1 - 0.05, omega_1), as in the paper (that
reading matches the WAMIT file; the other does not). HAMS values are |T| rho g with rho = 1025 kg/m^3 and g = 9.80665 m/s^2,
the values of both WAMIT runs: HAMS writes T / (rho g), and its runs scale the mass and the restoring by 1000/1025, which
keeps T / (rho g) of the 1025 kg/m^3 system.

Measure, per diagonal, over omega_1 in [LO, HI] rad/s: |X - WAMIT| / max |WAMIT|, largest and median, in percent. Nemoh
of Fig. 5 has its own frequency grid; WAMIT is interpolated linearly onto it (marked in the table).

Usage: python -I compare_diagonals.py CURVES_CSV FIGURE OUT_PNG OUT_TXT LABEL=HAMS_11D [LABEL=HAMS_11D ...]
       [--range LO HI] [--title TEXT]
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RHO_G = 1025.0 * 9.80665
WAMIT_FILE = {"5": "OpenFAST r-test marin_semi.12d", "7": "IEA-15-240-RWT .12d, IEA Wind Task 37"}
PAPER = ("WAMIT", "OrcaWave CS", "HydroStar CS", "OrcaWave PI", "HydroStar PI", "Nemoh")
STYLE = {"WAMIT": dict(color="#0b0b0b", lw=1.6, ls="-"),
         "OrcaWave CS": dict(color="#9a9893", lw=1.0, ls="-"),
         "HydroStar CS": dict(color="#9a9893", lw=1.0, ls=(0, (6, 2))),
         "OrcaWave PI": dict(color="#c3c2b7", lw=1.0, ls="-"),
         "HydroStar PI": dict(color="#c3c2b7", lw=1.0, ls=(0, (6, 2))),
         "Nemoh": dict(color="#c3c2b7", lw=1.0, ls=(0, (1, 1.5)))}
HAMS_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#e87ba4")
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"


def read_hams(path, mode=1):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) != 9 or float(p[2]) != 0.0 or float(p[3]) != 0.0 or int(p[4]) != mode:
            continue
        a, b = round(float(p[0]), 2), round(float(p[1]), 2)
        out[(a, b)] = out[(b, a)] = abs(complex(float(p[7]), float(p[8]))) * RHO_G
    return out


def main(argv):
    lo, hi, title = 0.30, 1.90, None
    if "--range" in argv:
        k = argv.index("--range")
        lo, hi = float(argv[k + 1]), float(argv[k + 2])
        del argv[k:k + 3]
    if "--title" in argv:
        k = argv.index("--title")
        title = argv[k + 1]
        del argv[k:k + 2]
    curves, fig, out_png, out_txt = argv[1:5]
    hams = [(spec.split("=", 1)[0], read_hams(spec.split("=", 1)[1])) for spec in argv[5:]]
    rows = [r for r in csv.DictReader(open(curves)) if r["figure"] == fig]
    paper = {}
    for r in rows:
        paper.setdefault((r["diagonal_dw"], r["solver"]), []).append((float(r["omega1"]), float(r["value_N_per_m2"])))
    platform = rows[0]["platform"]

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.facecolor": SURFACE, "figure.facecolor": SURFACE, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.0), gridspec_kw=dict(height_ratios=(3, 2)), sharex=True)
    lines = ["%s, |QTF| of the surge difference-frequency force against WAMIT (Roland et al. 2026, Fig. %s)" % (platform, fig),
             "|X - WAMIT| / max |WAMIT| over omega_1 in [%.2f, %.2f] rad/s, largest (at omega_1) / median, percent" % (lo, hi), ""]
    for col, dw in enumerate(("0.00", "0.05")):
        top, bottom = axes[0, col], axes[1, col]
        wam = sorted(paper[(dw, "WAMIT")])
        ww, wv = np.array([w for w, _ in wam]), np.array([v for _, v in wam])
        keep = (ww >= lo - 1e-9) & (ww <= hi + 1e-9)
        scale = wv[keep].max()
        lines.append("diagonal %s rad/s (max |WAMIT| %.4g N/m^2)" % (dw, scale))
        for solver in PAPER:
            pts = sorted(paper[(dw, solver)])
            x, y = np.array([w for w, _ in pts]), np.array([v for _, v in pts])
            inside = (x >= lo - 1e-9) & (x <= hi + 1e-9)   # the plots and the measure keep the HAMS range only
            x, y = x[inside], y[inside]
            source = "Roland et al. (2026) Fig. %s, %s" % (fig, solver)
            if solver == "WAMIT":
                source += " (%s)" % WAMIT_FILE[fig]
            top.plot(x, y, label=source, **STYLE[solver])
            if solver == "WAMIT":
                continue
            ref = np.interp(x, ww, wv)
            m = np.ones_like(x, bool)
            e = 100 * np.abs(y - ref) / scale
            bottom.plot(x, 100 * (y - ref) / scale, **STYLE[solver])
            lines.append("  %-32s %7.2f at %.2f / %5.2f%s" % (solver + " (paper)", e[m].max(), x[m][np.argmax(e[m])],
                                                              np.median(e[m]),
                                                      "  (WAMIT interpolated)" if solver == "Nemoh" and fig == "5" else ""))
        for k, (label, data) in enumerate(hams):
            d = float(dw)
            x = np.array(sorted({a for a, _ in data if lo - 1e-9 <= a <= hi + 1e-9 and (round(a - d, 2), a) in data}))
            y = np.array([data[(round(a - d, 2), a)] for a in x])
            c = HAMS_COLOURS[k]
            top.plot(x, y, color=c, lw=2.0, label=label, zorder=5)
            ref = np.interp(x, ww, wv)
            bottom.plot(x, 100 * (y - ref) / scale, color=c, lw=2.0, zorder=5)
            e = 100 * np.abs(y - ref) / scale
            lines.append("  %-32s %7.2f at %.2f / %5.2f" % (label, e.max(), x[np.argmax(e)], np.median(e)))
        lines.append("")
        name = r"$\omega_2 = \omega_1$ (mean drift $\bar F_1$)" if dw == "0.00" else r"$\omega_1 - \omega_2 = 0.05$ rad/s"
        top.set_title(r"$QTF_{11}$, %s" % name, loc="left", color=INK)
        top.set_ylabel(r"$|QTF_{11}|$ (N/m$^2$)")
        bottom.set_ylabel(r"$(|QTF_{11}| - |QTF_{11}^{WAMIT}|)$ / max$|QTF_{11}^{WAMIT}|$ (%)", fontsize=8)
        bottom.set_xlabel(r"$\omega_1$ (rad/s)")
        bottom.axhline(0.0, color=INK2, lw=0.8)
        bottom.set_ylim(-60, 60)
        top.set_xlim(lo - 0.03, hi + 0.03)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.02), ncol=2, fontsize=7.5, frameon=False,
                  title="Paper curves read from the vector paths of the PDF; HAMS-MREL from its output files",
                  title_fontsize=7.5)
    figure.suptitle(title or "%s: HAMS against the solvers of Roland et al. (2026), Fig. %s" % (platform, fig),
                    x=0.01, ha="left", color=INK)
    figure.text(0.01, 0.002, "Paper curves read from the vector paths of Roland et al. (2026), J. Phys.: Conf. Ser. 3346, "
                "012006, Fig. %s (CC BY 4.0), shown over the HAMS range %.2f to %.2f rad/s only. Lower plots clipped at "
                "+/-60 %%." % (fig, lo, hi),
                fontsize=7.5, color=INK2)
    figure.tight_layout(rect=(0, 0.14, 1, 0.97))
    figure.savefig(out_png, dpi=160)
    open(out_txt, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(list(sys.argv))
