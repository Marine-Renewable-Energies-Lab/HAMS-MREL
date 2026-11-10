"""
Compare the quadratic QTFs of HAMS-MREL with the exact reference of the sharp truncated cylinder and with Nemoh's DUOK part.

Usage:
    python qtf_compare.py check
    python qtf_compare.py <results> --reference <qtf_reference.csv> [--nemoh <folder of nemoh_* runs>] [--max-dw 0.3]
    python qtf_compare.py <results> --potential <qtf_potential.csv> [--nemoh <folder>] [--max-dw 0.3] [--exclude LO HI]

<results> holds the q_* case folders of `drift_convergence.py --qtf` (q_*/output/Wamit_format/QTF_Quadratic.10d, .10s).
The reference CSV is from `truncated_cylinder_reference.py qtf`: T/(rho g) for omega_i <= omega_j, e^{-i omega t}, both
kernels, F_x, F_z, M_y (moments about the origin, which is XR of the q_* cases). HAMS files hold conj(T)/(rho g): they are
conjugated back here. The rest of the square follows from T-_ji = conj(T-_ij), T+_ji = T+_ij.

Error measure. For one case, kernel and quantity, over the pairs omega_i <= omega_j of the grid (the difference kernel only
where omega_j - omega_i <= max-dw): |HAMS - exact| / max |exact| on each pair; the table gives its maximum and median over the
pairs, per mesh. Normalising by the largest value of the set keeps pairs where the QTF passes through zero from dominating.

Nemoh. OUT_QTFM_N.dat and OUT_QTFP_N.dat (w1 w2 beta1 beta2 DOF Mod Pha Re/(rho g) Im/(rho g)) list w1 >= w2. Nemoh's
documentation does not fix the conjugation of the off-diagonal values against this convention, so both readings, the value
as T(w1, w2), its conjugate, and both with the sign changed (a wave of another phase flips T+), are compared with the exact
reference, and the closest is reported and used, with a warning if even that one has a median error above 50%. Only the
documented output files are read.

Potential part (--potential). The CSV is from `truncated_cylinder_reference.py qtf-potential`: the quadratic T-_ij (the
columns of the quadratic reference, so --reference is not needed) and the five body-forcing terms of T_pot,ij
(FK, ID, B1, B2, B3, columns <quantity>_<term>_re/_im), difference kernel, omega_i <= omega_j. The reference total of
QTF_Indirect.11d is the quadratic part plus the terms. For each q_sharp_* case (`drift_convergence.py --qtf-potential`) it
compares QTF_Indirect.11d with that total, .11d - .10d with the sum of the terms, and each term of QTF_PotentialTerms.txt
(16 digits, conjugated back) with its reference; the scale of a term is the larger of its own and the potential part's
largest abs over the pairs, so that a term that is small everywhere is not measured against itself. Nemoh's OUT_QTFM_N.dat
of a DUOK + HASBO run (body forcing, no free-surface integral) is compared with the reference total, readings as above.
"""
import argparse
import csv
import glob
import os
import statistics
import sys

import numpy as np

import qtf_tools

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
    return errors_of(values, {p: ref[(motion, kernel) + p][quantity] for p in ps})


def errors_of(values, exact, scale=None):
    """(max, median, compared, pairs) of |value - exact| / scale over the pairs of exact {(wi, wj): complex}; scale
    defaults to the largest abs(exact). values maps (wi, wj) -> complex or None (not compared)."""
    ps = sorted(exact)
    scale = max(abs(v) for v in exact.values()) if scale is None else scale
    if scale == 0.0:
        return None
    errs = [abs(values[p] - exact[p]) / scale for p in ps if values.get(p) is not None]
    if not errs:
        return None
    return max(errs), statistics.median(errs), len(errs), len(ps)


POTENTIAL_TERMS = qtf_tools.POTENTIAL_TERMS


def read_potential_reference(path, max_dw):
    """{motion: {"total" | "pot" | term: {(w_i, w_j): {quantity: complex}}}} of a qtf-potential CSV, w_i <= w_j,
    w_j - w_i <= max_dw: the reference total of QTF_Indirect.11d (quadratic + terms), the potential part (sum of the
    terms) and each term, internal e^{-i omega t}, over rho g."""
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["kernel"] != "diff":
                continue
            wi, wj = round(float(r["omega_i"]), 4), round(float(r["omega_j"]), 4)
            if wj - wi > max_dw + 1e-9:
                continue
            m = out.setdefault(r["motion"], {key: {} for key in ("total", "pot") + POTENTIAL_TERMS})
            for q, _ in QUANTITIES:
                terms = {t: complex(float(r["%s_%s_re" % (q, t)]), float(r["%s_%s_im" % (q, t)])) for t in POTENTIAL_TERMS}
                pot = sum(terms.values())
                m["total"].setdefault((wi, wj), {})[q] = complex(float(r[q + "_re"]), float(r[q + "_im"])) + pot
                m["pot"].setdefault((wi, wj), {})[q] = pot
                for t in POTENTIAL_TERMS:
                    m[t].setdefault((wi, wj), {})[q] = terms[t]
    return out


def read_hams_terms(path):
    """{(w1, w2, I): [term] * 5} of QTF_PotentialTerms.txt, heading pair (0, 0), conjugated back to e^{-i omega t}."""
    if not os.path.exists(path):
        return {}
    return {(round(k[0], 4), round(k[1], 4), k[4]): [v.conjugate() for v in values]
            for k, values in qtf_tools.read_qtf_terms(path).items() if k[2] == 0.0 and k[3] == 0.0}


def _potential_cell(values, exact, scale=None):
    e = errors_of(values, exact, scale)
    return ("%7.2f %7.2f %s" % (100 * e[0], 100 * e[1], "" if e[2] == e[3] else "(%d/%d)" % e[2:])) if e else "-"


def _difference(a, b):
    """(wi, wj, index) -> a - b of two lookup functions, None where either is missing."""
    def value(wi, wj, index):
        x, y = a(wi, wj, index), b(wi, wj, index)
        return None if x is None or y is None else x - y
    return value


def _term(terms, t):
    """(wi, wj, index) -> term t of a read_hams_terms table; zero on the diagonal, where the file has no row."""
    table = {k: values[t] for k, values in terms.items()}

    def value(wi, wj, index):
        v = lookup(table, wi, wj, index, "diff")
        return 0j if v is None and wi == wj else v
    return value


def compare_potential(results, path, nemoh, max_dw, exclude=None):
    """The --potential tables; see the module docstring. exclude (lo, hi): leave out the pairs with omega_i or omega_j
    in [lo, hi] from the max and median, keeping the scales of all the pairs."""
    ref = read_potential_reference(path, max_dw)
    print("Potential part (qtf_terms = body), exact sharp cylinder: |HAMS - exact| / max |exact| over the pairs, max and "
          "median, percent; omega_j - omega_i <= %.2f rad/s; a term against the larger of its and T_pot's max" % max_dw)
    if exclude:
        print("Pairs with omega_i or omega_j in [%.3f, %.3f] rad/s left out; scales from all the pairs" % tuple(exclude))

    def kept(p):
        return not exclude or not any(exclude[0] - 1e-9 <= w <= exclude[1] + 1e-9 for w in p)
    print("%-24s %-6s %18s %18s %18s" % ("case", "what", "F_x", "F_z", "M_y"))
    for case in sorted(os.path.basename(p) for p in glob.glob(os.path.join(results, "q_sharp_*"))):
        motion = "fixed" if "_fixed_" in case else "free_damped" if "_damp_" in case else "free"
        if motion not in ref:
            print("%-24s no reference for the %s body in %s" % (case, motion, path))
            continue
        folder = os.path.join(results, case, "output", "Wamit_format")
        indirect = read_hams(os.path.join(folder, "QTF_Indirect.11d"))
        quadratic = read_hams(os.path.join(folder, "QTF_Quadratic.10d"))
        terms = read_hams_terms(os.path.join(folder, "QTF_PotentialTerms.txt"))
        if not indirect:
            print("%-24s no QTF_Indirect.11d" % case)
            continue

        def of(table):
            return lambda wi, wj, index: lookup(table, wi, wj, index, "diff")
        rows = [(".11d", of(indirect), "total", None)]
        if quadratic:
            rows.append(("pot", _difference(of(indirect), of(quadratic)), "pot", None))
        if terms:
            rows += [(name, _term(terms, t), name, "pot") for t, name in enumerate(POTENTIAL_TERMS)]
        for label, source, key, against in rows:
            cells = []
            for q, index in QUANTITIES:
                exact = {p: v[q] for p, v in ref[motion][key].items()}
                scale = max(abs(v) for v in exact.values()) if exclude else None
                if against:
                    scale = max(max(abs(v) for v in exact.values()), max(abs(v[q]) for v in ref[motion][against].values()))
                exact = {p: v for p, v in exact.items() if kept(p)}
                cells.append(_potential_cell({p: source(p[0], p[1], index) for p in exact}, exact, scale))
            print("%-24s %-6s %18s %18s %18s" % (case, label, *cells))
    if nemoh:
        print("\nNemoh DUOK + HASBO against the exact total (quadratic + body forcing); reading chosen by the smaller median")
        for folder in sorted(glob.glob(os.path.join(nemoh, "nemoh_sharp_*"))):
            name = os.path.basename(folder)
            motion = "fixed" if "_fixed_" in name else "free"
            table = read_nemoh(os.path.join(folder, "results", "QTF", "OUT_QTFM_N.dat"))
            if not table:
                continue
            if motion not in ref:
                print("%-24s no reference for the %s body in %s" % (name, motion, path))
                continue
            best = None
            for label, conjugate, sign in READINGS:
                cells, medians = [], []
                for q, index in QUANTITIES:
                    exact = {p: v[q] for p, v in ref[motion]["total"].items()}
                    scale = max(abs(v) for v in exact.values()) if exclude else None
                    exact = {p: v for p, v in exact.items() if kept(p)}
                    e = errors_of(nemoh_values(table, index, "diff", conjugate, sign), exact, scale)
                    cells.append("%7.2f %7.2f" % (100 * e[0], 100 * e[1]) if e else "-")
                    medians.append(e[1] if e else np.inf)
                if best is None or sum(medians) < best[0]:
                    best = (sum(medians), label, cells, max(medians))
            print("%-24s %-6s %16s %16s %16s%s" % (name, best[1], *best[2],
                                                    "  WARNING: no reading within 50%" if best[3] > 0.5 else ""))
    return 0


READINGS = (("as-is", False, 1.0), ("conj", True, 1.0), ("-1", False, -1.0), ("-conj", True, -1.0))


def nemoh_values(table, index, kernel, conjugate, sign=1.0):
    """(wi, wj) -> T(wi, wj) under the reading of Nemoh's value v as T(w1, w2) = sign * v, or sign * conj(v)."""
    out = {}
    for (w1, w2, dof), v in table.items():
        if dof != index:
            continue
        t = sign * (v.conjugate() if conjugate else v)
        out[(w1, w2)] = t
        out[(w2, w1)] = t.conjugate() if kernel == "diff" else t
    return out


def check_compare_potential():
    """compare_potential on a temporary fixture: a reference CSV with a quadratic part and five terms (a free body, two
    frequencies), and HAMS files written from it (.11d = quadratic + sum of the terms, .10d, QTF_PotentialTerms.txt,
    conjugated and in the file layouts): every error is 0 (to 1e-6, the 7 printed digits of .11d), so the reader of the
    reference CSV, the total = quadratic + terms assembly, the term order and the index mapping agree; a term moved to the
    wrong column gives a nonzero error; and a fixed case without a fixed reference prints its note."""
    import contextlib
    import io
    import tempfile
    pairs = ((0.7, 0.7), (0.7, 0.95), (0.95, 0.95))
    rng = np.random.default_rng(3)
    quad = {p: {q: complex(*rng.normal(size=2)) for q, _ in QUANTITIES} for p in pairs}
    terms = {p: {q: [complex(*rng.normal(size=2)) if p[0] != p[1] else 0j for _ in POTENTIAL_TERMS]
                 for q, _ in QUANTITIES} for p in pairs}

    def internal(wi, wj, q, what):
        if (wi, wj) in quad:
            return what(wi, wj, q)
        return what(wj, wi, q).conjugate()
    with tempfile.TemporaryDirectory() as d:
        ref_path = os.path.join(d, "ref.csv")
        cols = ["motion", "omega_i", "omega_j", "kernel"] + ["%s_%s" % (q, part) for q, _ in QUANTITIES for part in ("re", "im")]
        cols += ["%s_%s_%s" % (q, t, part) for q, _ in QUANTITIES for t in POTENTIAL_TERMS for part in ("re", "im")]
        with open(ref_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for wi, wj in pairs:
                row = ["free", wi, wj, "diff"]
                for q, _ in QUANTITIES:
                    row += [quad[wi, wj][q].real, quad[wi, wj][q].imag]
                for q, _ in QUANTITIES:
                    for t in range(len(POTENTIAL_TERMS)):
                        row += [terms[wi, wj][q][t].real, terms[wi, wj][q][t].imag]
                w.writerow(row)
        for case, swap in (("q_sharp_free_src_m1", None), ("q_sharp_free_src_m2", (1, 3)), ("q_sharp_fixed_src_m1", None)):
            folder = os.path.join(d, "results", case, "output", "Wamit_format")
            os.makedirs(folder)
            files = {name: open(os.path.join(folder, name), "w") for name in
                     ("QTF_Quadratic.10d", "QTF_Indirect.11d", "QTF_PotentialTerms.txt")}
            files["QTF_PotentialTerms.txt"].write("# header\n")
            for wi in (0.7, 0.95):
                for wj in (0.7, 0.95):
                    for q, index in QUANTITIES:
                        a = internal(wi, wj, q, lambda x, y, z: quad[x, y][z])
                        tt = [internal(wi, wj, q, lambda x, y, z, t=t: terms[x, y][z][t]) for t in range(5)]
                        for name, v in (("QTF_Quadratic.10d", a), ("QTF_Indirect.11d", a + sum(tt))):
                            v = v.conjugate()
                            files[name].write("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E\n" % (
                                wi, wj, 0.0, 0.0, index, abs(v), 0.0, v.real, v.imag))
                        if wi != wj:
                            if swap:
                                tt[swap[0]], tt[swap[1]] = tt[swap[1]], tt[swap[0]]
                            files["QTF_PotentialTerms.txt"].write("%14.6E%14.6E%14.6E%14.6E%6d" % (wi, wj, 0.0, 0.0, index)
                                                                  + "".join("%24.15E%24.15E" % (v.real, -v.imag)
                                                                            for v in tt) + "\n")
            for f in files.values():
                f.close()
        text = io.StringIO()
        with contextlib.redirect_stdout(text):
            compare_potential(os.path.join(d, "results"), ref_path, None, 0.3)
    lines = text.getvalue().splitlines()
    good = [l for l in lines if l.startswith("q_sharp_free_src_m1 ")]
    bad = [l for l in lines if l.startswith("q_sharp_free_src_m2 ") and l.split()[1] in ("ID", "B2")]
    numbers = [float(x) for l in good for x in l.split()[2:]]
    assert len(good) == 7 and max(numbers) <= 1e-4, "errors of the exact fixture: %s" % good
    assert bad and all(max(float(x) for x in l.split()[2:]) > 1.0 for l in bad), "swapped terms not seen: %s" % bad
    assert any("no reference for the fixed body" in l for l in lines), "no note for the fixed case"
    return "exact fixture max %.1e %%, swapped ID/B2 seen, fixed case noted" % max(numbers)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    if argv is None and sys.argv[1:2] == ["check"]:
        print("compare_potential    %s" % check_compare_potential())
        return 0
    parser.add_argument("results")
    parser.add_argument("--reference", help="CSV of truncated_cylinder_reference.py qtf (needed without --potential)")
    parser.add_argument("--potential", help="CSV of truncated_cylinder_reference.py qtf-potential: compare QTF_Indirect.11d "
                                            "and QTF_PotentialTerms.txt")
    parser.add_argument("--nemoh", help="folder holding nemoh_<hull>_<motion>_m<k>/results/QTF")
    parser.add_argument("--max-dw", type=float, default=0.3)
    parser.add_argument("--exclude", type=float, nargs=2, metavar=("LO", "HI"),
                        help="with --potential: leave out the pairs with a frequency in [LO, HI] rad/s (a resonance band)")
    args = parser.parse_args(argv)
    if args.potential:
        return compare_potential(args.results, args.potential, args.nemoh, args.max_dw, args.exclude)
    if not args.reference:
        parser.error("--reference is needed without --potential")
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
                for label, conjugate, sign in READINGS:
                    cells, medians = [], []
                    for quantity, index in QUANTITIES:
                        e = errors(nemoh_values(table, index, kernel, conjugate, sign), ref, motion, kernel, quantity,
                                   args.max_dw)
                        cells.append("%7.2f %7.2f" % (100 * e[0], 100 * e[1]) if e else "-")
                        medians.append(e[1] if e else np.inf)
                    score = sum(medians)
                    if best is None or score < best[0]:
                        best = (score, label, cells, max(medians))
                print("%-24s %-5s %-6s %16s %16s %16s%s" % (name, kernel, best[1], *best[2],
                                                            "  WARNING: no reading within 50%" if best[3] > 0.5 else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
