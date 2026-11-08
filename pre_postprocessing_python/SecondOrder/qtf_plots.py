"""
Plots of the quadratic QTFs of HAMS-MREL against the exact sharp-cylinder reference and Nemoh's DUOK part.

Usage:
    python qtf_plots.py <results> --reference <qtf_reference.csv> --out <folder> [--nemoh <folder of nemoh_* runs>]
                        [--dw 0.1]

For each motion (fixed, free) one figure per kernel, with a row per quantity (F_x, F_z, M_y over rho g, M_y in m) and two
columns, real and imaginary part, in the internal e^{-i omega t} convention:
- difference kernel: T-(omega, omega + dw) against omega, for the given dw;
- sum kernel: T+(omega, omega) against omega (the double-frequency load of a regular wave).
Every q_sharp_<motion>_src_m* case under <results> is a line, and each nemoh_sharp_<motion>_m* run under --nemoh.
The readers are those of qtf_compare.py.
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import qtf_compare as qc  # noqa: E402

PANELS = {"m1": 448, "m2": 1792, "m3": 4032, "m4": 7168}


def series(table, omegas, dw, index, kernel, source):
    """(x, values) along the plotted line, from a lookup function source(wi, wj, index, kernel) -> complex or None."""
    xs, vs = [], []
    for w in omegas:
        wj = round(w + (dw if kernel == "diff" else 0.0), 4)
        v = source(round(w, 4), wj, index, kernel)
        if v is not None:
            xs.append(w)
            vs.append(v)
    return xs, vs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("results")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--nemoh")
    parser.add_argument("--out", required=True)
    parser.add_argument("--dw", type=float, default=0.1)
    args = parser.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    ref = qc.read_reference(args.reference)
    written = []
    for motion in ("fixed", "free"):
        omegas = sorted({k[2] for k in ref if k[0] == motion})
        for kernel, ext, nfile in (("diff", "10d", "OUT_QTFM_N.dat"), ("sum", "10s", "OUT_QTFP_N.dat")):
            fig, axes = plt.subplots(3, 2, figsize=(11, 10), sharex=True)
            lines = []

            def ref_source(wi, wj, index, kern, motion=motion):
                q = {1: "Fx", 3: "Fz", 5: "My"}[index]
                if (motion, kern, wi, wj) in ref:
                    return ref[(motion, kern, wi, wj)][q]
                if (motion, kern, wj, wi) in ref:
                    v = ref[(motion, kern, wj, wi)][q]
                    return v.conjugate() if kern == "diff" else v
                return None
            lines.append(("exact", ref_source, dict(color="k", lw=2.2)))
            for path in sorted(glob.glob(os.path.join(args.results, "q_sharp_%s_src_m*" % motion))):
                tag = path.rsplit("_", 1)[-1]
                table = qc.read_hams(os.path.join(path, "output", "Wamit_format", "QTF_Quadratic." + ext))
                if table:
                    lines.append(("HAMS %d panels" % PANELS.get(tag, 0),
                                  lambda wi, wj, i, k, t=table: qc.lookup(t, wi, wj, i, k), dict(lw=1.4)))
            if args.nemoh:
                for path in sorted(glob.glob(os.path.join(args.nemoh, "nemoh_sharp_%s_m*" % motion))):
                    tag = path.rsplit("_", 1)[-1]
                    table = qc.read_nemoh(os.path.join(path, "results", "QTF", nfile))
                    if table:
                        lines.append(("Nemoh %d panels" % PANELS.get(tag, 0),
                                      lambda wi, wj, i, k, t=table: qc.lookup(t, wi, wj, i, k),
                                      dict(ls="--", marker="o", ms=3, lw=1.0)))
            for row, (quantity, index) in enumerate(qc.QUANTITIES):
                for label, source, style in lines:
                    xs, vs = series(None, omegas, args.dw, index, kernel, source)
                    axes[row, 0].plot(xs, [v.real for v in vs], label=label, **style)
                    axes[row, 1].plot(xs, [v.imag for v in vs], label=label, **style)
                axes[row, 0].set_ylabel("%s / (ρg)" % quantity)
                for col, part in enumerate(("Re", "Im")):
                    axes[row, col].set_title("%s %s" % (part, quantity), fontsize=10)
                    axes[row, col].grid(alpha=0.3)
            for col in range(2):
                axes[2, col].set_xlabel("ω (rad/s)" if kernel == "sum" else "ω_i (rad/s), ω_j = ω_i + %.2f" % args.dw)
            axes[0, 0].legend(fontsize=8)
            what = ("difference-frequency T⁻(ω, ω + %.2f)" % args.dw) if kernel == "diff" else "sum-frequency T⁺(ω, ω)"
            fig.suptitle("Sharp truncated cylinder, %s body: quadratic %s, e$^{-iωt}$, per unit amplitudes" % (motion, what))
            fig.tight_layout()
            name = os.path.join(args.out, "qtf_%s_%s.png" % (motion, kernel))
            fig.savefig(name, dpi=110)
            plt.close(fig)
            written.append(name)
    print("\n".join(written))
    return 0


if __name__ == "__main__":
    sys.exit(main())
