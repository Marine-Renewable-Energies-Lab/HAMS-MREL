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
Then one mesh-convergence figure per motion (convergence): the error against the exact reference against panel count.
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
    written += convergence(args, ref)
    print("\n".join(written))
    return 0


def convergence(args, ref, max_dw=0.3):
    """One figure per motion: the error of each mesh against the exact reference, |value - exact| / max |exact| over the
    pairs (qtf_compare.errors; difference pairs with omega_j - omega_i <= max_dw, sum pairs all), median (solid) and
    largest (dotted), against the panel count, for F_x, F_z, M_y and both kernels. HAMS filled, Nemoh open markers."""
    written = []
    for motion in ("fixed", "free"):
        fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
        for row, (kernel, ext, nfile) in enumerate((("diff", "10d", "OUT_QTFM_N.dat"), ("sum", "10s", "OUT_QTFP_N.dat"))):
            for col, (quantity, index) in enumerate(qc.QUANTITIES):
                ax = axes[row, col]
                sources = [("HAMS", os.path.join(args.results, "q_sharp_%s_src_m*" % motion), qc.read_hams,
                            os.path.join("output", "Wamit_format", "QTF_Quadratic." + ext), "o")]
                if args.nemoh:
                    sources.append(("Nemoh", os.path.join(args.nemoh, "nemoh_sharp_%s_m*" % motion), qc.read_nemoh,
                                    os.path.join("results", "QTF", nfile), "s"))
                for code, pattern, reader, fname, marker in sources:
                    pts = []
                    for path in sorted(glob.glob(pattern)):
                        n = PANELS.get(path.rsplit("_", 1)[-1])
                        table = reader(os.path.join(path, fname))
                        if not n or not table:
                            continue
                        values = {p: qc.lookup(table, p[0], p[1], index, kernel)
                                  for p in qc.pairs(ref, motion, kernel, max_dw)}
                        e = qc.errors(values, ref, motion, kernel, quantity, max_dw)
                        if e:
                            pts.append((n, 100 * e[1], 100 * e[0]))
                    if not pts:
                        continue
                    pts.sort()
                    face = None if code == "HAMS" else "none"
                    color = "C0" if code == "HAMS" else "C1"
                    ax.loglog([q[0] for q in pts], [q[1] for q in pts], "-" + marker, color=color, mfc=face,
                              label="%s median" % code)
                    ax.loglog([q[0] for q in pts], [q[2] for q in pts], ":" + marker, color=color, mfc=face,
                              label="%s largest" % code)
                ax.set_title("%s, %s" % (quantity, "difference" if kernel == "diff" else "sum"), fontsize=10)
                ax.grid(alpha=0.3, which="both")
                if col == 0:
                    ax.set_ylabel("error / max |exact| (%)")
                if row == 1:
                    ax.set_xlabel("panels")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle("Sharp truncated cylinder, %s body: mesh convergence of the quadratic QTFs against the exact reference\n"
                     "difference pairs with omega_j - omega_i <= %.1f rad/s, sum pairs all, 0.40 to 1.40 rad/s" % (motion, max_dw))
        fig.tight_layout()
        name = os.path.join(args.out, "qtf_%s_convergence.png" % motion)
        fig.savefig(name, dpi=110)
        plt.close(fig)
        written.append(name)
    return written


if __name__ == "__main__":
    sys.exit(main())
