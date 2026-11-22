"""Compare a HAMS-MREL case of make_cases.py with the HydroStar reference files of the Nemoh v3.0.3 test case, and plot.

Usage: python -I compare_case.py CASE HAMS_WAMIT_FORMAT_DIR HYDROSTAR_DIR OUT_DIR [--qtf FILE_PREFIX] [--kernel diff|sum]
    [--hams-file QTF_Indirect.11d] [--heads B1,B2] [--wmax 3.14]

Writes to OUT_DIR: <case>_first_order.png (A_ii, B_ii, |X_i|, |xi_i| for i = 1, 3, 5), <case>_qtf_<kernel>_<b1>_<b2>.png
(maps of |QTF_i|/(rho g) of HydroStar and HAMS and their difference, i = 1, 3, 5, as Figs. 14, 16 and 18 of Kurnia and
Ducrozet 2023), <case>_qtf_<kernel>_<b1>_<b2>_offdiag.png (Re and Im on the first off-diagonal for QTF_1 and QTF_3 and
the third for QTF_5, as Figs. 15, 17 and 19, plus the diagonal), and <case>_summary.txt with the errors.

Conventions (checked on OC4, all four heading pairs, readers.py): the HAMS row (omega_i, omega_j, beta1, beta2) equals the
HydroStar block "Heading1 = beta1, Heading2 = beta2", row omega_i, column dw = omega_j - omega_i, with no conjugation;
HydroStar is divided by its rho g = 1025 x 9.81, HAMS (non-dimensional, rho g A^2) is taken as it is. First-order HAMS
values are made dimensional with rho = 1025 and g = 9.81 to compare with HydroStar's: A = rho a, B = rho omega b,
X = rho g x; rotations in degrees. Errors: |HAMS - HydroStar| / max |HydroStar| over the pairs compared (both values
present, omega_i <= omega_j <= wmax), max and median, in percent; also for dw <= 0.6 rad/s.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import readers as R  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RG = R.RHO_HS * R.G_HS
MODES = ((1, "Fx", "QTF{s}_1", "surge"), (3, "Fz", "QTF{s}_3", "heave"), (5, "My", "QTF{s}_5", "pitch"))
SRC = "HydroStar (reference files of Nemoh v3.0.3 test case %s)"


def errors(h, r, sel):
    sc = np.nanmax(abs(r[sel]))
    d = abs(h[sel] - r[sel]) / sc * 100
    return np.nanmax(d), np.nanmedian(d), sc


def first_order(case, hdir, rdir, tc, out, heading=0.0):
    f = R.hams_first(hdir)
    fig, ax = plt.subplots(3, 4, figsize=(15, 9))
    lines = []
    for r, (i, _, _, name) in enumerate(MODES):
        unit_rot = i > 3
        for c, (kind, fn) in enumerate((("A", "AddedMass_%d%d.dat" % (i, i)), ("B", "Damping_%d%d.dat" % (i, i)),
                                        ("X", "%sf1st.rao" % ({1: "fx", 3: "fz", 5: "my"}[i])),
                                        ("xi", {1: "surge.rao", 3: "heave.rao", 5: "pitch.rao"}[i]))):
            path = os.path.join(rdir, fn)
            if not os.path.exists(path):
                continue
            w, a, p, heads, unit = R.hs_rao(path)
            k = heads.index(heading) if heading in heads else 0
            if kind in ("A", "B"):
                h = f["AB"][(i, i)]
                wh = h[:, 0]
                vh = h[:, 1] * R.RHO_HS if kind == "A" else h[:, 2] * R.RHO_HS * h[:, 0]
            else:
                wh, cplx = f["X" if kind == "X" else "RAO"][(round(heading, 3), i)]
                vh = abs(cplx) * (RG if kind == "X" else 1.0)
                if kind == "xi" and unit_rot:
                    vh = vh * 180.0 / np.pi
            ref = a[:, k]
            ax[r][c].plot(w, ref, "b-", lw=1.6, label=SRC % tc)
            ax[r][c].plot(wh, vh, "r--", lw=1.3, label="HAMS-MREL")
            sym = {"A": "A_%d%d" % (i, i), "B": "B_%d%d" % (i, i), "X": "|X_%d|" % i, "xi": "|ξ_%d|" % i}[kind]
            ax[r][c].set_title("%s [%s]" % (sym, unit), fontsize=9)
            ax[r][c].set_xlabel("ω [rad/s]", fontsize=8)
            ax[r][c].grid(alpha=0.3)
            ri = np.interp(wh, w, ref, left=np.nan, right=np.nan)
            m = np.isfinite(ri)
            d = abs(vh[m] - ri[m]) / np.nanmax(abs(ref)) * 100
            lines.append("%-6s max %6.2f %%  median %6.2f %% of max|HydroStar| (%d frequencies, %.3f-%.3f rad/s)"
                         % (sym, d.max(), np.median(d), m.sum(), wh[m].min(), wh[m].max()))
    ax[0][0].legend(fontsize=7)
    fig.suptitle("%s, first order, heading %g°: HAMS-MREL against %s" % (case, heading, SRC % tc), fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "%s_first_order.png" % case), dpi=100)
    plt.close(fig)
    return lines


def qtf(case, hfile, rdir, tc, out, prefix, kernel, heads, wmax, grid):
    s = "⁻" if kernel == "diff" else "⁺"
    lines = []
    fig, ax = plt.subplots(3, 3, figsize=(13, 11))
    fig2, ax2 = plt.subplots(3, 2, figsize=(13, 10))
    W1, W2 = np.meshgrid(grid, grid, indexing="ij")
    for r, (i, tag, sym, name) in enumerate(MODES):
        q = R.hs_qtf(os.path.join(rdir, "%s%s.rao" % (prefix, tag)))
        ref = R.hs_pair_matrix(q, heads, grid) / RG
        ham = R.hams_qtf(hfile, i, heads[0], heads[1], grid)
        sel = np.isfinite(ref) & np.isfinite(ham) & (W2 >= W1) & (W2 <= wmax + 1e-6)
        sel &= ~((ref == 0) & (W2 > W1))   # HydroStar writes 0 beyond its list
        emax, emed, sc = errors(ham, ref, sel)
        sel6 = sel & (W2 - W1 <= 0.6 + 1e-6)
        e6max, e6med, _ = errors(ham, ref, sel6)
        imax = np.unravel_index(np.nanargmax(np.where(sel, abs(ham - ref), -1)), ham.shape)
        lines.append("%s (%g°, %g°): all pairs to %.2f rad/s: max %.2f %% at (%.3f, %.3f), median %.2f %%; dw <= 0.6: "
                     "max %.2f %%, median %.2f %%; scale max|HydroStar| = %.4g (rho g); %d pairs"
                     % (sym.format(s=s), heads[0], heads[1], wmax, emax, grid[imax[0]], grid[imax[1]], emed, e6max, e6med,
                        sc, sel.sum()))
        # maps over the full square: the upper triangle as given, the lower by symmetry (|T| is symmetric)
        def full(m):
            a = np.where(sel, m, np.nan)
            return np.where(np.isfinite(a), a, a.T)
        rm, hm = abs(full(ref)), abs(full(ham))
        vmax = np.nanmax(rm)
        ext = [grid[0], grid[-1], grid[0], grid[-1]]
        for c, (m, title) in enumerate(((rm, "HydroStar"), (hm, "HAMS-MREL"),
                                        (rm - hm, "HydroStar − HAMS-MREL"))):
            kw = dict(vmin=0, vmax=vmax, cmap="viridis") if c < 2 else dict(vmin=-0.2 * vmax, vmax=0.2 * vmax, cmap="RdBu_r")
            im = ax[r][c].imshow(m, origin="lower", extent=ext, aspect="equal", **kw)
            plt.colorbar(im, ax=ax[r][c], fraction=0.046)
            ax[r][c].set_title("|%s|/ρg, %s" % (sym.format(s=s), title) + (" (±20 %% scale)" if c == 2 else ""), fontsize=9)
            ax[r][c].set_xlabel("ω_2 [rad/s]", fontsize=8)
            ax[r][c].set_ylabel("ω_1 [rad/s]", fontsize=8)
            ax[r][c].set_xlim(grid[0], wmax)
            ax[r][c].set_ylim(grid[0], wmax)
        # off-diagonal lines: first for modes 1 and 3, third for 5 (the paper's choice), and the diagonal
        off = 3 if i == 5 else 1
        for c, k in enumerate((off, 0)):
            idx = [(a, a + k) for a in range(len(grid) - k) if sel[a, a + k]]
            w = [grid[a] for a, _ in idx]
            ax2[r][c].plot(w, [ref[a, b].real for a, b in idx], "b-", label="HydroStar Re")
            ax2[r][c].plot(w, [ref[a, b].imag for a, b in idx], "b-.", label="HydroStar Im")
            ax2[r][c].plot(w, [ham[a, b].real for a, b in idx], "r--", label="HAMS-MREL Re")
            ax2[r][c].plot(w, [ham[a, b].imag for a, b in idx], "r:", label="HAMS-MREL Im")
            what = "diagonal (mean drift)" if k == 0 else "off-diagonal %d (ω_2 = ω_1 + %.4f rad/s)" % (k, k * 0.0628)
            ax2[r][c].set_title("%s/ρg, %s" % (sym.format(s=s), what), fontsize=9)
            ax2[r][c].set_xlabel("ω_1 [rad/s]", fontsize=8)
            ax2[r][c].grid(alpha=0.3)
    ax2[0][0].legend(fontsize=7)
    title = "%s, %s QTF%s, headings (β_1, β_2) = (%g°, %g°): HAMS-MREL (%s) against %s" % (
        case, "difference" if kernel == "diff" else "sum", s, heads[0], heads[1], os.path.basename(hfile), SRC % tc)
    for fg, nm in ((fig, ""), (fig2, "_offdiag")):
        fg.suptitle(title, fontsize=10)
        fg.tight_layout()
        fg.savefig(os.path.join(out, "%s_qtf_%s_%g_%g%s.png" % (case, kernel, heads[0], heads[1], nm)), dpi=95)
        plt.close(fg)
    return lines


def main():
    p = argparse.ArgumentParser()
    p.add_argument("case")
    p.add_argument("hdir")
    p.add_argument("rdir")
    p.add_argument("out")
    p.add_argument("--tc", default="")
    p.add_argument("--qtf", default=None)
    p.add_argument("--kernel", default="diff")
    p.add_argument("--hams-file", default="QTF_Indirect.11d")
    p.add_argument("--heads", default="0,0")
    p.add_argument("--wmax", type=float, default=3.14)
    p.add_argument("--nfreq", type=int, default=50)
    p.add_argument("--no-first", action="store_true")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    grid = [round(0.0628 * (k + 1), 4) for k in range(a.nfreq)]
    lines = []
    if not a.no_first:
        lines += first_order(a.case, a.hdir, a.rdir, a.tc, a.out)
    if a.qtf:
        heads = tuple(float(v) for v in a.heads.split(","))
        lines += qtf(a.case, os.path.join(a.hdir, a.hams_file), a.rdir, a.tc, a.out, a.qtf, a.kernel, heads, a.wmax, grid)
    with open(os.path.join(a.out, "%s_summary_%s_%s.txt" % (a.case, a.kernel, a.heads.replace(",", "_"))), "w",
              encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
