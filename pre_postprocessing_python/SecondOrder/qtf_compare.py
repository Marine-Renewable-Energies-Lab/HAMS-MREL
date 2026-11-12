"""
Compare the quadratic QTFs of HAMS-MREL with the exact reference of the sharp truncated cylinder and with Nemoh's DUOK part.

Usage:
    python qtf_compare.py check
    python qtf_compare.py <results> --reference <qtf_reference.csv> [--nemoh <folder of nemoh_* runs>] [--max-dw 0.3]
    python qtf_compare.py <results> --potential <qtf_potential.csv> [--nemoh <folder>] [--max-dw 0.3] [--exclude LO HI]
    python qtf_compare.py <results> --against-nemoh <folder> [--hull oc4] [--reading as-is|conj|-1|-conj] [--max-dw 0.3]
                          [--exclude LO HI]
    python qtf_compare.py <results> --direct [<qtf_potential.csv>] [--max-dw 0.3] [--exclude LO HI] [--by-omega]
    python qtf_compare.py <results> --against-wamit <file.12d> [--hull oc4w] [--reading as-is|conj|-1|-conj]
                          [--max-dw 0.3] [--max-omega 1.90] [--exclude LO HI]

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

Direct method (--direct). For each q_*direct* case (`drift_convergence.py --qtf-potential`, qtf_direct = 1), T_dir
(QTF_Direct.12d - .10d) against T_pot (QTF_Indirect.11d - .10d) of the same run: |T_dir - T_pot| / max |T_pot| over the
rows of forces or of moments (all modes, bodies and heading pairs), max and median, the measure of the integration tests.
They are read from QTF_DirectTerms.txt and QTF_PotentialTerms.txt (16 digits) where the run wrote them, otherwise from the
printed files. --by-omega adds the largest error at each difference frequency (an irregular frequency of the Omega problem
shows there as a spike). Given the qtf-potential CSV, the q_sharp_*direct* cases also get .12d and .11d against the exact
total and .12d - .10d and .11d - .10d against the exact potential part, as --potential, heading 0.

WAMIT (--against-wamit). A WAMIT numeric .12d (difference QTF; one header line, then PER_i PER_j BETA_i BETA_j I Mod Pha Re
Im, periods in s, nondimensional by rho g L^k with L = 1 m: the layout HAMS writes for .10d and .11d) is read "as written":
omega = 2 pi / PER (snapped to the HAMS grid, as the periods have 5 digits), and the value taken in the convention of the
HAMS files (WAMIT's e^{+i omega t}, which HAMS writes by conjugating), so both are conjugated back alike. It holds one half
of the square (omega_i >= omega_j); the other half is its conjugate. For each q_<hull>_* case, QTF_Indirect.11d (quadratic
+ body forcing) and QTF_Quadratic.10d (quadratic only) are compared with it, |HAMS - WAMIT| / max |WAMIT|, max and median,
over the pairs of the HAMS file with omega <= max-omega: those off the diagonal, those with 0 < omega_j - omega_i <= max-dw,
and the diagonal (the mean drift, real) alone. The reading is fixed ("as-is" by default) and never chosen by agreement;
another reading is only named in a WARNING when it fits much better off the diagonal, as for --against-nemoh. The diagonal
fixes the sign (a reading "-1" fails there) but not the conjugation (the diagonal is real); the median ratio Re(HAMS) /
Re(WAMIT) on the diagonal is printed per quantity.
"""
import argparse
import csv
import glob
import os
import statistics
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # run with python -I: the sibling modules

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


def read_part_sums(path, parts):
    """{(oufr1, oufr2, beta1, beta2, I): sum of the parts} of a terms file with `parts` complex columns after the five keys
    (QTF_PotentialTerms.txt 5, QTF_DirectTerms.txt 2), in the file convention, keys rounded as qtf_tools.read_qtf; {} if
    the file is missing."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                continue
            p = line.split()
            if len(p) != 5 + 2 * parts:
                continue
            key = (round(float(p[0]), 6), round(float(p[1]), 6), round(float(p[2]), 6), round(float(p[3]), 6), int(p[4]))
            out[key] = sum(complex(float(p[5 + 2 * t]), float(p[6 + 2 * t])) for t in range(parts))
    return out


def _is_force(key):
    return (key[4] - 1) % 6 < 3


def direct_against_indirect(t_dir, t_pot, kept=None):
    """(forces, moments): (max, median) of |T_dir - T_pot| over max |T_pot| of the rows of that kind (forces: I - 1 mod 6
    < 3), over the keys of t_pot (all modes, bodies and heading pairs) for which kept(key) is true; None for a kind without
    rows. The measure of the integration tests (direct_against_indirect of test_integration_drift.f90)."""
    out = []
    for force in (True, False):
        keys = [k for k in t_pot if _is_force(k) == force]
        missing = [k for k in keys if k not in t_dir]
        if missing:
            raise ValueError("T_dir has no row for %d rows of T_pot, e.g. %s" % (len(missing), missing[0]))
        scale = max([abs(t_pot[k]) for k in keys] or [0.0])
        errs = [abs(t_dir[k] - t_pot[k]) / scale for k in keys if kept is None or kept(k)] if scale > 0.0 else []
        out.append((max(errs), statistics.median(errs)) if errs else None)
    return tuple(out)


def compare_direct(results, path, max_dw, exclude=None, by_omega=False):
    """The --direct tables; see the module docstring. path: CSV of truncated_cylinder_reference.py qtf-potential, or None.
    Returns the rows (case, (forces, moments) of T_dir against T_pot)."""
    ref = read_potential_reference(path, max_dw) if path else {}
    print("Direct method (QTF_Direct.12d, qtf_direct = 1) against the indirect one (QTF_Indirect.11d) of the same run: "
          "|T_dir - T_pot| / max |T_pot| over the rows of forces or of moments, all modes and heading pairs, max and median, "
          "percent; omega_j - omega_i <= %.2f rad/s" % max_dw)
    if exclude:
        print("Pairs with omega_i or omega_j in [%.3f, %.3f] rad/s left out; scales from all the pairs" % tuple(exclude))

    def kept(p):
        return not exclude or not any(exclude[0] - 1e-9 <= w <= exclude[1] + 1e-9 for w in p)
    print("%-26s %-6s %18s %18s" % ("case", "from", "forces", "moments"))
    rows, sharp = [], []
    for case in sorted(os.path.basename(p) for p in glob.glob(os.path.join(results, "q_*direct*"))):
        folder = os.path.join(results, case, "output", "Wamit_format")
        # The terms files carry 16 digits; without them (write_terms = 0) the printed .12d - .10d and .11d - .10d
        t_dir = read_part_sums(os.path.join(folder, "QTF_DirectTerms.txt"), 2)
        t_pot = read_part_sums(os.path.join(folder, "QTF_PotentialTerms.txt"), len(POTENTIAL_TERMS))
        source = "terms"
        if not t_dir or not t_pot:
            files = [os.path.join(folder, name) for name in ("QTF_Quadratic.10d", "QTF_Direct.12d", "QTF_Indirect.11d")]
            if not all(os.path.exists(f) for f in files):
                print("%-26s no QTF_Direct.12d, QTF_Indirect.11d and QTF_Quadratic.10d" % case)
                continue
            quadratic = qtf_tools.read_qtf(files[0])
            t_dir = qtf_tools.potential_part(quadratic, qtf_tools.read_qtf(files[1]))
            t_pot = qtf_tools.potential_part(quadratic, qtf_tools.read_qtf(files[2]))
            source = "files"
        t_pot = {k: v for k, v in t_pot.items() if k[0] != k[1] and abs(k[1] - k[0]) <= max_dw + 1e-9}
        if not t_pot:
            print("%-26s no pairs within %.2f rad/s" % (case, max_dw))
            continue
        cells = direct_against_indirect(t_dir, t_pot, lambda k: kept(k[:2]))
        rows.append((case, cells))
        print("%-26s %-6s %18s %18s" % (case, source,
                                        *["%7.2f %7.2f" % (100 * c[0], 100 * c[1]) if c else "-" for c in cells]))
        if by_omega:
            # The largest error at each difference frequency, over the scales of all the pairs: a spike at an irregular
            # frequency of the Omega problem shows here
            scales = [max([abs(v) for k, v in t_pot.items() if _is_force(k) == force] or [0.0]) for force in (True, False)]
            for omega in sorted({round(abs(k[1] - k[0]), 6) for k in t_pot}):
                at = [k for k in t_pot if abs(round(abs(k[1] - k[0]), 6) - omega) < 1e-9 and kept(k[:2])]
                worst = [max([abs(t_dir[k] - t_pot[k]) for k in at if _is_force(k) == force] or [0.0]) / scales[i]
                         if scales[i] > 0.0 else 0.0 for i, force in enumerate((True, False))]
                print("    Omega %8.4f  %4d rows  forces %7.2f  moments %7.2f" % (omega, len(at), 100 * worst[0],
                                                                                  100 * worst[1]))
        if case.startswith("q_sharp_"):
            sharp.append((case, folder))
    if not ref:
        return rows
    print("\nAgainst the exact sharp cylinder (heading 0): |HAMS - exact| / max |exact| over the pairs, max and median, "
          "percent; .12d and .11d against the exact total, dir (.12d - .10d) and pot (.11d - .10d) against the exact "
          "potential part")
    print("%-26s %-5s %18s %18s %18s" % ("case", "what", "F_x", "F_z", "M_y"))
    for case, folder in sharp:
        motion = "fixed" if "_fixed_" in case else "free_damped" if "_damp_" in case else "free"
        if motion not in ref:
            print("%-26s no reference for the %s body in %s" % (case, motion, path))
            continue
        direct = read_hams(os.path.join(folder, "QTF_Direct.12d"))
        indirect = read_hams(os.path.join(folder, "QTF_Indirect.11d"))
        quadratic = read_hams(os.path.join(folder, "QTF_Quadratic.10d"))

        def of(table):
            return lambda wi, wj, index: lookup(table, wi, wj, index, "diff")
        sources = [(".12d", of(direct), "total"), (".11d", of(indirect), "total")]
        if quadratic:
            sources += [("dir", _difference(of(direct), of(quadratic)), "pot"),
                        ("pot", _difference(of(indirect), of(quadratic)), "pot")]
        for label, source, key in sources:
            cells = []
            for q, index in QUANTITIES:
                exact = {p: v[q] for p, v in ref[motion][key].items()}
                scale = max(abs(v) for v in exact.values())
                exact = {p: v for p, v in exact.items() if kept(p)}
                cells.append(_potential_cell({p: source(p[0], p[1], index) for p in exact}, exact, scale))
            print("%-26s %-5s %18s %18s %18s" % (case, label, *cells))
    return rows


READINGS = (("as-is", False, 1.0), ("conj", True, 1.0), ("-1", False, -1.0), ("-conj", True, -1.0))


def compare_against_nemoh(results, nemoh, max_dw, exclude=None, hull="oc4", reading="as-is"):
    """The --against-nemoh table, for a hull without an exact reference: QTF_Indirect.11d of each q_<hull>_<motion>_m<k>
    case against OUT_QTFM_N.dat of nemoh_<hull>_<motion>_m<k> (DUOK + HASBO, the same mesh), |HAMS - Nemoh| / max |Nemoh|
    over the pairs omega_i <= omega_j of the HAMS file with omega_j - omega_i <= max_dw, max and median, percent. Nemoh's
    values are taken under one fixed reading, "as-is" by default: the one the sharp-cylinder runs established against the
    exact reference. With no reference here, the reading is never chosen by agreement with HAMS, which would hide a sign or
    conjugation error of either code; the other readings are only computed to print a WARNING when one of them has a sum of
    medians under WARN_READING_RATIO of the fixed reading's. exclude (lo, hi) as in compare_potential.
    Returns the rows (case, [(max, median) or None per quantity], warning or None)."""
    fixed = [r for r in READINGS if r[0] == reading]
    if not fixed:
        raise ValueError("unknown Nemoh reading " + str(reading))
    print("HAMS QTF_Indirect.11d against Nemoh DUOK + HASBO, %s: |HAMS - Nemoh| / max |Nemoh| over the pairs, max and median, "
          "percent; omega_j - omega_i <= %.2f rad/s; Nemoh reading %s" % (hull, max_dw, reading))
    if exclude:
        print("Pairs with omega_i or omega_j in [%.3f, %.3f] rad/s left out; scales from all the pairs" % tuple(exclude))
    print("%-24s %16s %16s %16s" % ("case", "F_x", "F_z", "M_y"))
    rows = []
    for folder in sorted(glob.glob(os.path.join(results, "q_%s_*" % hull))):
        case = os.path.basename(folder)
        partner = os.path.join(nemoh, "nemoh_" + case[2:])
        hams = read_hams(os.path.join(folder, "output", "Wamit_format", "QTF_Indirect.11d"))
        table = read_nemoh(os.path.join(partner, "results", "QTF", "OUT_QTFM_N.dat"))
        if not hams or not table:
            print("%-24s %s" % (case, "no QTF_Indirect.11d" if not hams else "no %s" % os.path.basename(partner)))
            continue
        pairs_all = sorted(set((min(w1, w2), max(w1, w2)) for w1, w2, _ in hams if abs(w2 - w1) <= max_dw + 1e-9))
        kept = [p for p in pairs_all if not exclude or not any(exclude[0] - 1e-9 <= w <= exclude[1] + 1e-9 for w in p)]

        def cells_of(conjugate, sign):
            cells = []
            for _, index in QUANTITIES:
                values = nemoh_values(table, index, "diff", conjugate, sign)
                if not all(p in values for p in pairs_all):
                    cells.append(None)
                    continue
                scale = max(abs(values[p]) for p in pairs_all)
                e = errors_of({p: lookup(hams, p[0], p[1], index, "diff") for p in kept}, {p: values[p] for p in kept}, scale)
                cells.append(e[:2] if e else None)
            return cells, sum(c[1] if c else np.inf for c in cells)
        cells, score = cells_of(*fixed[0][1:])
        warning = None
        for label, conjugate, sign in READINGS:
            if label == reading:
                continue
            other = cells_of(conjugate, sign)[1]
            if other < WARN_READING_RATIO * score:
                warning = "WARNING: reading %s fits much better (sum of medians %.2f%% against %.2f%%)" % (
                    label, 100 * other, 100 * score)
        rows.append((case, cells, warning))
        print("%-24s %16s %16s %16s%s" % (case, *["%7.2f %7.2f" % (100 * c[0], 100 * c[1]) if c else "-" for c in cells],
                                          "  " + warning if warning else ""))
    return rows


# A reading other than the fixed one is reported when its sum of medians is below this fraction of the fixed reading's
WARN_READING_RATIO = 0.2


def check_against_nemoh():
    """compare_against_nemoh on a temporary fixture: a HAMS .11d (conj(T), both halves of the square) and a Nemoh file (T as
    written, w1 >= w2) of the same T: errors 0 under the default "as-is" reading. Nemoh written conjugated: a large error under
    "as-is", with a warning naming "conj", and errors 0 with reading = "conj". A Nemoh run missing is noted."""
    import contextlib
    import io
    import tempfile
    grid = (0.4, 0.5, 0.6)
    rng = np.random.default_rng(5)
    t = {(wi, wj): {index: complex(*rng.normal(size=2)) for _, index in QUANTITIES} for wi in grid for wj in grid if wi <= wj}
    for p in t:
        if p[0] == p[1]:
            t[p] = {index: complex(v.real, 0.0) for index, v in t[p].items()}

    def value(wi, wj, index):
        return t[(wi, wj)][index] if wi <= wj else t[(wj, wi)][index].conjugate()
    with tempfile.TemporaryDirectory() as d:
        for motion, conjugate in (("fixed", False), ("free", True)):
            folder = os.path.join(d, "results", "q_oc4_%s_m1" % motion, "output", "Wamit_format")
            os.makedirs(folder)
            with open(os.path.join(folder, "QTF_Indirect.11d"), "w") as f:
                for wi in grid:
                    for wj in grid:
                        for _, index in QUANTITIES:
                            v = value(wi, wj, index).conjugate()
                            f.write("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E\n" % (
                                wi, wj, 0.0, 0.0, index, abs(v), 0.0, v.real, v.imag))
            folder = os.path.join(d, "nemoh", "nemoh_oc4_%s_m1" % motion, "results", "QTF")
            os.makedirs(folder)
            with open(os.path.join(folder, "OUT_QTFM_N.dat"), "w") as f:
                f.write("# header\n")
                for w1 in grid:
                    for w2 in grid:
                        if w1 < w2:
                            continue
                        for _, index in QUANTITIES:
                            v = value(w1, w2, index)
                            v = v.conjugate() if conjugate else v
                            f.write("%.6E %.6E 0. 0. %d %.6E 0. %.9E %.9E\n" % (w1, w2, index, abs(v), v.real, v.imag))
        os.makedirs(os.path.join(d, "results", "q_oc4_free_m2", "output", "Wamit_format"))
        with open(os.path.join(d, "results", "q_oc4_free_m2", "output", "Wamit_format", "QTF_Indirect.11d"), "w") as f:
            f.write("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E\n" % (0.4, 0.4, 0.0, 0.0, 1, 1.0, 0.0, 1.0, 0.0))
        text = io.StringIO()
        with contextlib.redirect_stdout(text):
            rows = compare_against_nemoh(os.path.join(d, "results"), os.path.join(d, "nemoh"), 0.3)
            conj_rows = compare_against_nemoh(os.path.join(d, "results"), os.path.join(d, "nemoh"), 0.3, reading="conj")
    found = {case: (cells, warning) for case, cells, warning in rows}
    as_is = max(c[0] for c in found["q_oc4_fixed_m1"][0])
    wrong = min(c[1] for c in found["q_oc4_free_m1"][0])
    conj = max(c[0] for case, cells, _ in conj_rows if case == "q_oc4_free_m1" for c in cells)
    assert as_is <= 1e-5 and found["q_oc4_fixed_m1"][1] is None, "as-is fixture: %s" % (found["q_oc4_fixed_m1"],)
    assert wrong > 0.1 and found["q_oc4_free_m1"][1] and "conj" in found["q_oc4_free_m1"][1], \
        "the conjugated fixture under as-is: %s" % (found["q_oc4_free_m1"],)
    assert conj <= 1e-5, "the conjugated fixture under conj: %g" % conj
    assert "Nemoh reading as-is" in text.getvalue() and "no nemoh_oc4_free_m2" in text.getvalue(), "header or missing note"
    return ("as-is fixture max %.1e %%, conjugated fixture %.0f %% median under as-is with a warning and %.1e %% under conj, "
            "missing run noted" % (100 * as_is, 100 * wrong, 100 * conj))


def read_wamit(path):
    """{(w1, w2, I): T/(rho g)} in e^{-i omega t}, heading pair (0, 0), from a WAMIT numeric .12d: omega = 2 pi / PER
    (unrounded; snap_wamit puts it on a grid), and the value as written, conjugated back as read_hams does with a HAMS file."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) != 9:
                continue
            try:
                vals = [float(x) for x in p]
            except ValueError:
                continue
            if vals[2] != 0.0 or vals[3] != 0.0:
                continue
            out[(2.0 * np.pi / vals[0], 2.0 * np.pi / vals[1], int(vals[4]))] = complex(vals[7], -vals[8])
    return out


# A WAMIT omega goes to the nearest omega of the HAMS grid within this (rad/s): the 5 digits of PER leave about 1e-4 at 3 rad/s
WAMIT_SNAP = 1e-3


def snap_wamit(table, grid):
    """table of read_wamit with each omega replaced by the nearest one of grid (within WAMIT_SNAP); rows with an omega off the
    grid are dropped. Returns (snapped table, number of grid omega found in the file)."""
    grid = sorted(set(grid))

    def snapped(w):
        best = min(grid, key=lambda g: abs(g - w)) if grid else None
        return best if best is not None and abs(best - w) <= WAMIT_SNAP else None
    out, found = {}, set()
    for (w1, w2, index), v in table.items():
        a, b = snapped(w1), snapped(w2)
        if a is None or b is None:
            continue
        found.update((a, b))
        out[(a, b, index)] = v
    return out, len(found)


def compare_against_wamit(results, path, max_dw, max_omega=1.90, exclude=None, hull="oc4w", reading="as-is"):
    """The --against-wamit table: QTF_Indirect.11d and QTF_Quadratic.10d of each q_<hull>_* case against the WAMIT .12d at
    path; see the module docstring. exclude (lo, hi) as in compare_potential. Returns the rows (case, file, set, [(max,
    median) or None per quantity], warning or None, [median Re(HAMS) / Re(WAMIT) on the diagonal, or None, per quantity])."""
    fixed = [r for r in READINGS if r[0] == reading]
    if not fixed:
        raise ValueError("unknown WAMIT reading " + str(reading))
    raw = read_wamit(path)
    if not raw:
        print("no WAMIT rows in %s" % path)
        return []
    print("HAMS against WAMIT %s, %s: |HAMS - WAMIT| / max |WAMIT| over the pairs of the set, max and median, percent; "
          "omega <= %.2f rad/s; WAMIT reading %s (as written, conjugated back like the HAMS files)" % (
              os.path.basename(path), hull, max_omega, reading))
    if exclude:
        print("Pairs with omega_i or omega_j in [%.3f, %.3f] rad/s left out; scales from all the pairs of the set"
              % tuple(exclude))
    print("%-20s %-5s %-9s %16s %16s %16s   diagonal Re(HAMS) / Re(WAMIT), median" % ("case", "file", "set", "F_x", "F_z",
                                                                                     "M_y"))
    rows = []
    for folder in sorted(glob.glob(os.path.join(results, "q_%s_*" % hull))):
        case = os.path.basename(folder)
        files = [(label, read_hams(os.path.join(folder, "output", "Wamit_format", name)))
                 for label, name in ((".11d", "QTF_Indirect.11d"), (".10d", "QTF_Quadratic.10d"))]
        files = [(label, table) for label, table in files if table]
        if not files:
            print("%-20s no QTF_Indirect.11d or QTF_Quadratic.10d" % case)
            continue
        grid = sorted({w for _, table in files for key in table for w in key[:2] if w <= max_omega + 1e-9})
        wamit, on_grid = snap_wamit(raw, grid)
        if on_grid < len(grid):
            print("%-20s only %d of the %d HAMS omega are in the WAMIT file" % (case, on_grid, len(grid)))
        for label, hams in files:
            pairs_all = sorted(set((min(w1, w2), max(w1, w2)) for w1, w2, _ in hams if max(w1, w2) <= max_omega + 1e-9))
            sets = (("all", [p for p in pairs_all if p[0] != p[1]]),
                    ("dw<=%.2f" % max_dw, [p for p in pairs_all if p[0] != p[1] and p[1] - p[0] <= max_dw + 1e-9]),
                    ("diagonal", [p for p in pairs_all if p[0] == p[1]]))

            def cells_of(chosen, conjugate, sign):
                cells = []
                for _, index in QUANTITIES:
                    values = nemoh_values(wamit, index, "diff", conjugate, sign)
                    if not chosen or not all(p in values for p in chosen):
                        cells.append(None)
                        continue
                    scale = max(abs(values[p]) for p in chosen)
                    kept = [p for p in chosen
                            if not exclude or not any(exclude[0] - 1e-9 <= w <= exclude[1] + 1e-9 for w in p)]
                    e = errors_of({p: lookup(hams, p[0], p[1], index, "diff") for p in kept},
                                  {p: values[p] for p in kept}, scale)
                    cells.append(e[:2] if e else None)
                return cells, sum(c[1] if c else np.inf for c in cells)
            ratios = []
            for _, index in QUANTITIES:
                values = nemoh_values(wamit, index, "diff", *fixed[0][1:])
                r = []
                for p in sets[2][1]:
                    h = lookup(hams, p[0], p[1], index, "diff")
                    if h is not None and p in values and values[p].real != 0.0:
                        r.append(h.real / values[p].real)
                ratios.append(statistics.median(r) if r else None)
            for name, chosen in sets:
                cells, score = cells_of(chosen, *fixed[0][1:])
                warning = None
                if name == "all":
                    for other_label, conjugate, sign in READINGS:
                        if other_label == reading:
                            continue
                        other = cells_of(chosen, conjugate, sign)[1]
                        if other < WARN_READING_RATIO * score:
                            warning = "WARNING: reading %s fits much better (sum of medians %.2f%% against %.2f%%)" % (
                                other_label, 100 * other, 100 * score)
                rows.append((case, label, name, cells, warning, ratios))
                diagonal = "   " + " ".join("%.3f" % r if r is not None else "-" for r in ratios) if name == "diagonal" else ""
                print("%-20s %-5s %-9s %16s %16s %16s%s%s" % (
                    case, label, name, *["%7.2f %7.2f" % (100 * c[0], 100 * c[1]) if c else "-" for c in cells], diagonal,
                    "  " + warning if warning else ""))
    return rows


def _wamit_period(period):
    """A period as WAMIT writes it: 0.ddddd E+nn, 5 digits."""
    exponent = int(np.floor(np.log10(period))) + 1
    mantissa = round(period / 10.0 ** exponent, 5)
    if mantissa >= 1.0:
        mantissa, exponent = mantissa / 10.0, exponent + 1
    return "%.5fE%+03d" % (mantissa, exponent)


def check_against_wamit():
    """compare_against_wamit on a temporary fixture: a HAMS .11d and .10d (conj(T), the whole square; .10d differs from .11d
    off the diagonal only) and a WAMIT .12d of the .11d (one header line, periods with 5 digits, the half omega_i >= omega_j,
    modes in WAMIT's order 1 3 5 2 4 6, values written as HAMS writes them, and an omega that HAMS does not have): .11d errors
    0 under "as-is" in every set, diagonal ratio 1; .10d errors 0 on the diagonal only; a WAMIT file written conjugated gives
    a warning naming "conj" and errors 0 under reading "conj"; one with the sign flipped fails on the diagonal (ratio -1); a
    case without files is noted."""
    import contextlib
    import io
    import tempfile
    grid = (0.25, 0.30, 0.35, 1.85)
    rng = np.random.default_rng(7)
    t = {(wi, wj): {index: complex(*rng.normal(size=2)) for index in range(1, 7)} for wi in grid for wj in grid if wi <= wj}
    pot = {p: {index: (0j if p[0] == p[1] else complex(*rng.normal(size=2))) for index in range(1, 7)} for p in t}
    for p in t:
        if p[0] == p[1]:
            t[p] = {index: complex(v.real, 0.0) for index, v in t[p].items()}

    def value(wi, wj, index, table):
        return table[(wi, wj)][index] if wi <= wj else table[(wj, wi)][index].conjugate()

    def row(a, b, index, v):
        v = v.conjugate()   # written as HAMS writes: conj(T)
        return "%15s%15s%15.5E%15.5E%5d%15.5E%15.5E%15.5E%15.5E\n" % (a, b, 0.0, 0.0, index, abs(v), np.degrees(np.angle(v)),
                                                                     v.real, v.imag)
    with tempfile.TemporaryDirectory() as d:
        folder = os.path.join(d, "results", "q_oc4w_free_m1", "output", "Wamit_format")
        os.makedirs(folder)
        with open(os.path.join(folder, "QTF_Indirect.11d"), "w") as f11, \
                open(os.path.join(folder, "QTF_Quadratic.10d"), "w") as f10:
            for wi in grid:
                for wj in grid:
                    for index in range(1, 7):
                        total = value(wi, wj, index, t)
                        f11.write(row("%.6E" % wi, "%.6E" % wj, index, total))
                        f10.write(row("%.6E" % wi, "%.6E" % wj, index, total - value(wi, wj, index, pot)))
        os.makedirs(os.path.join(d, "results", "q_oc4w_free_m2"))
        paths = {}
        for name, conjugate, sign in (("asis", False, 1.0), ("conj", True, 1.0), ("minus", False, -1.0)):
            paths[name] = os.path.join(d, "%s.12d" % name)
            with open(paths[name], "w") as f:
                f.write(" WAMIT Numeric Output -- Filename  fixture.12d\n")
                for w1 in grid + (0.40,):
                    for w2 in grid + (0.40,):
                        if w1 < w2:
                            continue
                        for index in (1, 3, 5, 2, 4, 6):
                            v = sign * (complex(1.0, 0.0) if 0.40 in (w1, w2) else value(w1, w2, index, t))
                            f.write(row(_wamit_period(2 * np.pi / w1), _wamit_period(2 * np.pi / w2), index,
                                        v.conjugate() if conjugate else v))
        text = io.StringIO()
        with contextlib.redirect_stdout(text):
            rows = compare_against_wamit(os.path.join(d, "results"), paths["asis"], 0.3)
            conj_rows = compare_against_wamit(os.path.join(d, "results"), paths["conj"], 0.3)
            conj_fixed = compare_against_wamit(os.path.join(d, "results"), paths["conj"], 0.3, reading="conj")
            minus_rows = compare_against_wamit(os.path.join(d, "results"), paths["minus"], 0.3)
    found = {(r[1], r[2]): r for r in rows}
    worst = max(c[0] for key in ((".11d", "all"), (".11d", "dw<=0.30"), (".11d", "diagonal")) for c in found[key][3])
    ratio = max(abs(r - 1.0) for r in found[(".11d", "diagonal")][5])
    diag10 = max(c[0] for c in found[(".10d", "diagonal")][3])
    off10 = min(c[1] for c in found[(".10d", "all")][3])
    conj_warning = [r[4] for r in conj_rows if r[1] == ".11d" and r[2] == "all"][0]
    conj_worst = max(c[0] for r in conj_fixed if r[1] == ".11d" for c in r[3])
    minus = [r for r in minus_rows if r[1] == ".11d" and r[2] == "diagonal"][0]
    assert worst <= 1e-4 and ratio <= 1e-4 and found[(".11d", "all")][4] is None, "as-is fixture: %s" % (rows,)
    assert diag10 <= 1e-4 and off10 > 0.1, ".10d fixture: diagonal %g, off the diagonal %g" % (diag10, off10)
    assert conj_warning and "conj" in conj_warning and conj_worst <= 1e-4, "conjugated fixture: %s" % (conj_rows,)
    assert min(c[0] for c in minus[3]) > 1.9 and max(abs(r + 1.0) for r in minus[5]) <= 1e-4, "sign fixture: %s" % (minus,)
    assert "no QTF_Indirect.11d or QTF_Quadratic.10d" in text.getvalue(), "missing case not noted"
    return (".11d fixture max %.1e %%, diagonal ratio 1 to %.0e, .10d %.0f %% median off the diagonal, conjugated file "
            "warned and %.1e %% under conj, flipped sign seen on the diagonal, missing case noted"
            % (100 * worst, ratio, 100 * off10, 100 * conj_worst))


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


def check_compare_direct():
    """compare_direct on a temporary fixture: a reference CSV (free body, quadratic part and five terms, two frequencies) and
    three cases. q_sharp_direct_free_m1 has T_dir = 1.01 T_pot in its terms files and its .12d: the direct-indirect error is
    1% at the largest |T_pot| (to the 16 digits of the terms files), .11d against the exact total is 0 and .12d is not.
    q_sharp_direct_free_m2 has T_dir = -T_pot: 200%. q_round_direct_free_m1 has no terms files, two heading pairs and
    .12d = .11d: read from the printed files, 0 to their 7 digits. --by-omega prints one line per difference frequency."""
    import contextlib
    import io
    import tempfile
    pairs = ((0.7, 0.7), (0.7, 0.95), (0.95, 0.95))
    rng = np.random.default_rng(5)
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
        for case, factor, headings, with_terms in (("q_sharp_direct_free_m1", 1.01, (0.0,), True),
                                                   ("q_sharp_direct_free_m2", -1.0, (0.0,), True),
                                                   ("q_round_direct_free_m1", 1.0, (0.0, 30.0), False)):
            folder = os.path.join(d, "results", case, "output", "Wamit_format")
            os.makedirs(folder)
            names = ["QTF_Quadratic.10d", "QTF_Indirect.11d", "QTF_Direct.12d"]
            if with_terms:
                names += ["QTF_PotentialTerms.txt", "QTF_DirectTerms.txt"]
            files = {name: open(os.path.join(folder, name), "w") for name in names}
            for wi in (0.7, 0.95):
                for wj in (0.7, 0.95):
                    for beta in headings:
                        for q, index in QUANTITIES:
                            a = internal(wi, wj, q, lambda x, y, z: quad[x, y][z])
                            tt = [internal(wi, wj, q, lambda x, y, z, t=t: terms[x, y][z][t]) for t in range(5)]
                            if beta:   # another heading pair: other values, the same structure
                                a, tt = 2.0 * a, [0.5 * v for v in tt]
                            pot = sum(tt)
                            for name, v in (("QTF_Quadratic.10d", a), ("QTF_Indirect.11d", a + pot),
                                            ("QTF_Direct.12d", a + factor * pot)):
                                v = v.conjugate()
                                files[name].write("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E\n" % (
                                    wi, wj, beta, beta, index, abs(v), 0.0, v.real, v.imag))
                            if with_terms and wi != wj:
                                key = "%14.6E%14.6E%14.6E%14.6E%6d" % (wi, wj, beta, beta, index)
                                files["QTF_PotentialTerms.txt"].write(key + "".join(
                                    "%24.15E%24.15E" % (v.real, -v.imag) for v in tt) + "\n")
                                parts = (factor * tt[0], factor * (pot - tt[0]))
                                files["QTF_DirectTerms.txt"].write(key + "".join(
                                    "%24.15E%24.15E" % (v.real, -v.imag) for v in parts) + "\n")
            for f in files.values():
                f.close()
        text = io.StringIO()
        with contextlib.redirect_stdout(text):
            rows = dict(compare_direct(os.path.join(d, "results"), ref_path, 0.3, by_omega=True))
    lines = text.getvalue().splitlines()
    one, flipped, round_ = rows["q_sharp_direct_free_m1"], rows["q_sharp_direct_free_m2"], rows["q_round_direct_free_m1"]
    assert all(abs(c[0] - 0.01) < 1e-12 for c in one), "T_dir = 1.01 T_pot: %s" % (one,)
    assert all(abs(c[0] - 2.0) < 1e-12 for c in flipped), "T_dir = -T_pot: %s" % (flipped,)
    assert all(c[0] < 1e-5 for c in round_), "T_dir = T_pot, printed files: %s" % (round_,)
    assert any(l.startswith("q_round_direct_free_m1") and " files " in l for l in lines), "round case not from the files"
    assert sum(1 for l in lines if l.strip().startswith("Omega")) == 3, "one --by-omega line per case and Omega"

    def cells(case, what):
        found = [l for l in lines if l.startswith(case + " ") and l.split()[1] == what]
        assert len(found) == 1, "no %s row of %s" % (what, case)
        return [float(x) for x in found[0].split()[2:]]
    assert max(cells("q_sharp_direct_free_m1", ".11d")) <= 1e-4, "the .11d of the fixture is the exact total"
    assert max(cells("q_sharp_direct_free_m1", "pot")) <= 1e-4, "the .11d - .10d of the fixture is the exact potential part"
    assert max(cells("q_sharp_direct_free_m1", "dir")[0::2]) > 0.5, "1% of T_pot off the exact potential part"
    assert not any(l.startswith("q_round_direct") for l in lines[lines.index(next(l for l in lines if "exact sharp" in l)):]), \
        "the rounded hull has no exact reference"
    return "1.01 T_pot gives %.4f%%, -T_pot %.1f%%, printed files %.1e%%, exact rows of the sharp cases only" % (
        100 * one[0][0], 100 * flipped[0][0], 100 * max(c[0] for c in round_))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    if argv is None and sys.argv[1:2] == ["check"]:
        print("compare_potential    %s" % check_compare_potential())
        print("against_nemoh        %s" % check_against_nemoh())
        print("against_wamit        %s" % check_against_wamit())
        print("compare_direct       %s" % check_compare_direct())
        return 0
    parser.add_argument("results")
    parser.add_argument("--reference", help="CSV of truncated_cylinder_reference.py qtf (needed without --potential)")
    parser.add_argument("--potential", help="CSV of truncated_cylinder_reference.py qtf-potential: compare QTF_Indirect.11d "
                                            "and QTF_PotentialTerms.txt")
    parser.add_argument("--nemoh", help="folder holding nemoh_<hull>_<motion>_m<k>/results/QTF")
    parser.add_argument("--max-dw", type=float, default=0.3)
    parser.add_argument("--exclude", type=float, nargs=2, metavar=("LO", "HI"),
                        help="with --potential or --against-nemoh: leave out the pairs with a frequency in [LO, HI] rad/s "
                             "(a resonance band)")
    parser.add_argument("--against-nemoh", metavar="FOLDER",
                        help="compare QTF_Indirect.11d of q_<hull>_* with Nemoh's nemoh_<hull>_* under FOLDER (no exact "
                             "reference)")
    parser.add_argument("--direct", nargs="?", const="", metavar="CSV",
                        help="compare QTF_Direct.12d of the q_*direct* cases with QTF_Indirect.11d of the same run, and, "
                             "given the CSV of truncated_cylinder_reference.py qtf-potential, both with the exact reference")
    parser.add_argument("--by-omega", action="store_true",
                        help="with --direct: also the largest direct-indirect error at each difference frequency")
    parser.add_argument("--against-wamit", metavar="FILE",
                        help="compare QTF_Indirect.11d and QTF_Quadratic.10d of q_<hull>_* with a WAMIT .12d file")
    parser.add_argument("--hull", help="with --against-nemoh or --against-wamit: the hull of the case names (default oc4, "
                                       "or oc4w with --against-wamit)")
    parser.add_argument("--max-omega", type=float, default=1.90,
                        help="with --against-wamit: the largest omega of the pairs compared (default 1.90 rad/s)")
    parser.add_argument("--reading", choices=[label for label, _, _ in READINGS], default="as-is",
                        help="with --against-nemoh: Nemoh's reading (default as-is, established on the sharp cylinder "
                             "against the exact reference; never chosen by agreement with HAMS); with --against-wamit: "
                             "WAMIT's (default as-is, as written)")
    args = parser.parse_args(argv)
    if args.direct is not None:
        rows = compare_direct(args.results, args.direct or None, args.max_dw, args.exclude, args.by_omega)
        return 0 if rows else 1
    if args.against_wamit:
        rows = compare_against_wamit(args.results, args.against_wamit, args.max_dw, args.max_omega, args.exclude,
                                     args.hull or "oc4w", args.reading)
        return 0 if rows else 1
    if args.against_nemoh:
        rows = compare_against_nemoh(args.results, args.against_nemoh, args.max_dw, args.exclude, args.hull or "oc4",
                                     args.reading)
        return 0 if rows else 1
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
