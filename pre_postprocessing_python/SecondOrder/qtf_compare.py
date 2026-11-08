"""
Compare the quadratic QTFs of HAMS-MREL with the exact reference of the sharp truncated cylinder and with Nemoh's DUOK part.

Usage:
    python qtf_compare.py <results> --reference <qtf_reference.csv> [--nemoh <folder of nemoh_* runs>] [--max-dw 0.3]

<results> holds the q_* case folders of `drift_convergence.py --qtf` (q_*/output/Wamit_format/QTF_Quadratic.10d, .10s).
The reference CSV is from `truncated_cylinder_reference.py qtf`: T/(rho g) for omega_i <= omega_j, e^{-i omega t}, both
kernels, F_x, F_z, M_y (moments about the origin, which is XR of the q_* cases). HAMS files hold conj(T)/(rho g): they are
conjugated back here. The rest of the square follows from T-_ji = conj(T-_ij), T+_ji = T+_ij.

Error measure. For one case, kernel and quantity, over the pairs omega_i <= omega_j of the grid (the difference kernel only
where omega_j - omega_i <= max-dw): |HAMS - exact| / max |exact| on each pair; the table gives its maximum and median over the
pairs, per mesh. Normalising by the largest value of the set keeps pairs where the QTF passes through zero from dominating.

Nemoh. OUT_QTFM_N.dat and OUT_QTFP_N.dat (w1 w2 beta1 beta2 DOF Mod Pha Re/(rho g) Im/(rho g)) list w1 >= w2. Nemoh's
documentation does not fix the conjugation of the off-diagonal values against this convention, so both readings, the value
as T(w1, w2) and as its conjugate, are compared with the exact reference and the closer one is reported and used. Only the
documented output files are read.
"""
import argparse
import csv
import glob
import os
import statistics
import sys

import numpy as np

QUANTITIES = (("Fx", 1), ("Fz", 3), ("My", 5))


def read_reference(path):
    """{(motion, kernel, w_i, w_j): {quantity: complex}} for w_i <= w_j."""
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            key = (r["motion"], r["kernel"], round(float(r["omega_i"]), 4), round(float(r["omega_j"]), 4))
            out[key] = {q: complex(float(r[q + "_re"]), float(r[q + "_im"])) for q, _ in QUANTITIES}
    return out


def read_hams(path):
    """{(w1, w2, I): T/(rho g)} in e^{-i omega t}, heading pair (0, 0), from a QTF_Quadratic file."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) != 9 or float(p[2]) != 0.0 or float(p[3]) != 0.0:
                continue
            out[(round(float(p[0]), 4), round(float(p[1]), 4), int(p[4]))] = complex(float(p[7]), -float(p[8]))
    return out


def read_nemoh(path):
    """{(w1, w2, dof): value/(rho g)} as written, heading pair (0, 0)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            p = line.split()
            try:
                vals = [float(x) for x in p]
            except ValueError:
                continue
            if len(vals) != 9 or vals[2] != 0.0 or vals[3] != 0.0:
                continue
            out[(round(vals[0], 4), round(vals[1], 4), int(vals[4]))] = complex(vals[7], vals[8])
    return out


def lookup(table, wi, wj, index, kernel):
    """T(wi, wj) from a table holding either half of the square."""
    if (wi, wj, index) in table:
        return table[(wi, wj, index)]
    if (wj, wi, index) in table:
        v = table[(wj, wi, index)]
        return v.conjugate() if kernel == "diff" else v
    return None


def pairs(ref, motion, kernel, max_dw):
    keys = sorted(k for k in ref if k[0] == motion and k[1] == kernel)
    return [(k[2], k[3]) for k in keys if kernel == "sum" or k[3] - k[2] <= max_dw + 1e-9]


def errors(values, ref, motion, kernel, quantity, max_dw):
    """(max, median) of |value - exact| / max |exact| over the pairs; values maps (wi, wj) -> complex or None."""
    ps = pairs(ref, motion, kernel, max_dw)
    exact = np.array([ref[(motion, kernel, wi, wj)][quantity] for wi, wj in ps])
    scale = np.max(np.abs(exact))
    errs = [abs(values[p] - e) / scale for p, e in zip(ps, exact) if values.get(p) is not None]
    if not errs:
        return None
    return max(errs), statistics.median(errs), len(errs), len(ps)


def nemoh_values(table, index, kernel, conjugate):
    """(wi, wj) -> T(wi, wj) under the reading of Nemoh's value as T(w1, w2) (conjugate False) or conj T(w1, w2)."""
    out = {}
    for (w1, w2, dof), v in table.items():
        if dof != index:
            continue
        t = v.conjugate() if conjugate else v
        out[(w1, w2)] = t
        out[(w2, w1)] = t.conjugate() if kernel == "diff" else t
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("results")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--nemoh", help="folder holding nemoh_<hull>_<motion>_m<k>/results/QTF")
    parser.add_argument("--max-dw", type=float, default=0.3)
    args = parser.parse_args(argv)
    ref = read_reference(args.reference)
    cases = sorted(os.path.basename(p) for p in glob.glob(os.path.join(args.results, "q_sharp_*")))
    print("Exact sharp cylinder: |HAMS - exact| / max |exact| over the pairs, max and median, percent; "
          "difference kernel for omega_j - omega_i <= %.2f rad/s" % args.max_dw)
    print("%-24s %-5s %-3s %18s %18s %18s" % ("case", "kern", "", "F_x", "F_z", "M_y"))
    for case in cases:
        motion = "fixed" if "_fixed_" in case else "free"
        for kernel, ext in (("diff", "10d"), ("sum", "10s")):
            table = read_hams(os.path.join(args.results, case, "output", "Wamit_format", "QTF_Quadratic." + ext))
            if not table:
                continue
            cells = []
            for quantity, index in QUANTITIES:
                values = {p: lookup(table, p[0], p[1], index, kernel) for p in pairs(ref, motion, kernel, args.max_dw)}
                e = errors(values, ref, motion, kernel, quantity, args.max_dw)
                cells.append("%7.2f %7.2f %s" % (100 * e[0], 100 * e[1], "" if e[2] == e[3] else "(%d/%d)" % e[2:]) if e
                             else "-")
            print("%-24s %-5s %-3s %18s %18s %18s" % (case, kernel, "", *cells))
    if args.nemoh:
        print("\nNemoh DUOK against the exact reference (same measure); reading chosen per kernel by the smaller median")
        for folder in sorted(glob.glob(os.path.join(args.nemoh, "nemoh_sharp_*"))):
            name = os.path.basename(folder)
            motion = "fixed" if "_fixed_" in name else "free"
            for kernel, fname in (("diff", "OUT_QTFM_N.dat"), ("sum", "OUT_QTFP_N.dat")):
                table = read_nemoh(os.path.join(folder, "results", "QTF", fname))
                if not table:
                    continue
                best = None
                for conjugate in (False, True):
                    cells, medians = [], []
                    for quantity, index in QUANTITIES:
                        e = errors(nemoh_values(table, index, kernel, conjugate), ref, motion, kernel, quantity, args.max_dw)
                        cells.append("%7.2f %7.2f" % (100 * e[0], 100 * e[1]) if e else "-")
                        medians.append(e[1] if e else np.inf)
                    score = sum(medians)
                    if best is None or score < best[0]:
                        best = (score, conjugate, cells)
                print("%-24s %-5s %-4s %16s %16s %16s" % (name, kernel, "conj" if best[1] else "as-is", *best[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
