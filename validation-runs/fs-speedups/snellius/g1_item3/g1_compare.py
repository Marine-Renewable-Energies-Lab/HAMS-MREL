"""G1 item 3 after the Snellius runs: qtf_fs_circle = panels against samples, each against the exact references.

Usage: python g1_compare.py <run folder> <repo> --sum-ref <qtfdirectsum_reference.csv> --diff-ref <qtffs_reference.csv>
                            [--out <folder>] [--routes samples panels]

<run folder> holds g1_<route>/<case> (run_one.sh). The expected cases are the lines of sharp.txt and round.txt beside this
script, not what is on disk. Writes g1_compare.txt (and g1_compare.png with matplotlib) to --out (default: the run folder).

Like with like (qtf_compare.py is imported for its readers, never edited; its --direct mode compares the qtf_terms = full
.12d/.11d of the directsum cases with the reference parts that leave out the free surface, review C2 of Task 7). Both
references hold, per motion (fixed, free) and pair omega_i <= omega_j, the quadratic part, the five body-forcing terms and
the free-surface part FS (columns <q>_FS_re/_im). qtffs_reference.csv: kernel diff, 0.40-1.40 by 0.01, omega_j - omega_i
<= 0.3. qtfdirectsum_reference.csv: kernel sum, 0.40-1.40 by 0.05. read_potential_reference gives "total" (quadratic +
terms), "pot" (terms), "FS", "total_fs" = total + FS, "pot_fs" = pot + FS. HAMS side (sharp m2 cases, heading 0,
e^{-i omega t}):
- difference: FS_ind = FS_inner + FS_outer (QTF_PotentialTerms.txt, terms 6 + 7) and D_FS (QTF_DirectTerms.txt, part 3)
  against FS; .11d - FS_ind and .12d - D_FS against total; P_BODY (the five terms) and D_BODY (D_I2 + D_D) against pot;
  .11d and .12d against total_fs.
- sum: P_FS and D_FS (QTF_DirectTermsSum.txt, parts 5 and 3) against FS; .12s - D_FS against total; P_BODY and D_BODY
  against pot; .12s against total_fs.
Each: |HAMS - exact| / max |exact| over the pairs the run has (F_x, F_z, M_y), max and median, percent; the free body also
with 0.75-0.80 rad/s left out (the surge-pitch resonance, as compare_stage5_excl080.txt of 7c-2). Then the change
|new - today| / max |exact|, and the G1 rule per quantity on the max and on the median: the new route passes when its error
is at most today's + max(0.5 pp, 10 % of today's). Round m2 (no exact reference): the in-run direct-indirect check of both
kernels per route. Then every output file of the new route against today's (qtf_files_compare.py) and the core-hours.

Any problem (a case or route missing, an output file missing or without rows, no pair in common with a reference, a missing
result.json or a nonzero code) is listed under PROBLEMS at the end and the script exits 1.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import qtf_files_compare as qfc  # noqa: E402

QUANTITIES = (("Fx", 1), ("Fz", 3), ("My", 5))
EXCLUDE = (0.75, 0.80)


def expected_cases():
    cases = []
    for name in ("sharp.txt", "round.txt"):
        for line in open(os.path.join(HERE, name)):
            if line.strip():
                case = line.split()[1]
                if case not in cases:
                    cases.append(case)
    return cases


def read_terms(path):
    """{(w1, w2, I): [complex]} of a terms file (heading pair 0, 0), keys rounded to 4 digits, conjugated back to
    e^{-i omega t} (the files hold conj(T)/(rho g), as .10d/.10s)."""
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.lstrip().startswith("#"):
            continue
        p = line.split()
        if len(p) < 7 or (len(p) - 5) % 2:
            continue
        try:
            v = [float(x) for x in p]
        except ValueError:
            continue
        if v[2] != 0.0 or v[3] != 0.0:
            continue
        out[(round(v[0], 4), round(v[1], 4), int(v[4]))] = [complex(v[5 + 2 * t], -v[6 + 2 * t])
                                                            for t in range((len(v) - 5) // 2)]
    return out


class Sources:
    """The HAMS quantities of one case and kernel: [(name, function (wi, wj, index) -> complex or None, reference key)]."""

    FILES = {"diff": ("QTF_Quadratic.10d", "QTF_Indirect.11d", "QTF_Direct.12d", "QTF_PotentialTerms.txt",
                      "QTF_DirectTerms.txt"),
             "sum": ("QTF_Quadratic.10s", "QTF_Direct.12s", "QTF_DirectTermsSum.txt")}

    def __init__(self, qc, folder, kernel, problems, tag):
        self.qc, self.kernel = qc, kernel
        self.t = {}
        for name in self.FILES[kernel]:
            path = os.path.join(folder, name)
            table = read_terms(path) if name.endswith(".txt") else qc.read_hams(path)
            if not table:
                problems.append("%s: %s missing or without rows" % (tag, name))
            self.t[name] = table
        self.ok = all(self.t.values())

    def of(self, name, parts=None):
        table = self.t[name]
        if parts is not None:
            table = {k: sum(v[i] for i in parts) for k, v in table.items()}

        def value(wi, wj, index):
            v = self.qc.lookup(table, wi, wj, index, self.kernel)
            # the difference terms files have no diagonal row: the potential part is zero there (Omega = 0)
            if v is None and parts is not None and self.kernel == "diff" and wi == wj:
                return 0j
            return v
        return value

    def items(self):
        d = self.qc._difference
        if self.kernel == "diff":
            fs_ind = self.of("QTF_PotentialTerms.txt", (5, 6))
            d_fs = self.of("QTF_DirectTerms.txt", (2,))
            return [("FS_ind", fs_ind, "FS"), ("D_FS", d_fs, "FS"),
                    (".11d-FS", d(self.of("QTF_Indirect.11d"), fs_ind), "total"),
                    (".12d-D_FS", d(self.of("QTF_Direct.12d"), d_fs), "total"),
                    ("P_BODY", self.of("QTF_PotentialTerms.txt", (0, 1, 2, 3, 4)), "pot"),
                    ("D_BODY", self.of("QTF_DirectTerms.txt", (0, 1)), "pot"),
                    (".11d", self.of("QTF_Indirect.11d"), "total_fs"),
                    (".12d", self.of("QTF_Direct.12d"), "total_fs")]
        d_fs = self.of("QTF_DirectTermsSum.txt", (2,))
        return [("P_FS", self.of("QTF_DirectTermsSum.txt", (4,)), "FS"), ("D_FS", d_fs, "FS"),
                (".12s-D_FS", d(self.of("QTF_Direct.12s"), d_fs), "total"),
                ("P_BODY", self.of("QTF_DirectTermsSum.txt", (3,)), "pot"),
                ("D_BODY", self.of("QTF_DirectTermsSum.txt", (0, 1)), "pot"),
                (".12s", self.of("QTF_Direct.12s"), "total_fs")]


def stats(errs):
    errs = sorted(errs)
    n = len(errs)
    if not n:
        return None
    med = errs[n // 2] if n % 2 else 0.5 * (errs[n // 2 - 1] + errs[n // 2])
    return 100 * errs[-1], 100 * med


def margin(today):
    return max(0.5, 0.1 * today)


def against_exact(f, exact, kept, index):
    """(max, median) percent of |f - exact| / max |exact| over the pairs where f has a value (scale over all of them,
    statistics over the kept ones); None if no pair."""
    vals = {p: f(p[0], p[1], index) for p in exact}
    have = [p for p in exact if vals[p] is not None]
    if not have:
        return None
    scale = max(abs(exact[p]) for p in have)
    return stats([abs(vals[p] - exact[p]) / scale for p in have if kept(p)]) if scale > 0 else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("repo")
    parser.add_argument("--sum-ref", required=True)
    parser.add_argument("--diff-ref", required=True)
    parser.add_argument("--out")
    parser.add_argument("--routes", nargs=2, default=("samples", "panels"), metavar=("TODAY", "NEW"),
                        help="the two route folders g1_<TODAY>, g1_<NEW> (default samples panels)")
    args = parser.parse_args()
    sys.path.insert(0, os.path.join(args.repo, "pre_postprocessing_python", "SecondOrder"))
    import qtf_compare as qc
    out = args.out or args.run
    today, new = args.routes
    problems, lines, verdict = [], [], []
    refs = {"diff": qc.read_potential_reference(args.diff_ref, 0.3, "diff"),
            "sum": qc.read_potential_reference(args.sum_ref, 0.0, "sum")}
    for kernel, ref in refs.items():
        for motion in ("fixed", "free"):
            if motion not in ref or not ref[motion]["FS"]:
                problems.append("reference %s has no %s rows with the FS part" % (kernel, motion))
    cases = expected_cases()
    folders = {}
    for route in (today, new):
        for case in cases:
            folders[(route, case)] = os.path.join(args.run, "g1_" + route, case, "output", "Wamit_format")
            if not os.path.isdir(folders[(route, case)]):
                problems.append("%s / %s: no output/Wamit_format" % (route, case))

    # 1. sharp cases against the exact references, and the G1 rule
    lines.append("1. Against the exact sharp cylinder, |HAMS - exact| / max |exact|, max / median, percent; change = "
                 "|%s - %s| / max |exact|; G1: %s <= %s + max(0.5 pp, 10 %% of %s), on the max and on the median"
                 % (new, today, new, today, today))
    lines.append("%-27s %-4s %-4s %-10s %-3s %17s %17s %17s  %s" % ("case", "kern", "set", "quantity", "q", today, new,
                                                                     "change", "G1 max / median"))
    for case in [c for c in cases if c.startswith("q_sharp_")]:
        motion = "fixed" if "_fixed_" in case else "free"
        for kernel in ("diff", "sum"):
            ref = refs[kernel].get(motion)
            src = {}
            for route in (today, new):
                if os.path.isdir(folders[(route, case)]):
                    s = Sources(qc, folders[(route, case)], kernel, problems, "%s / %s" % (route, case))
                    if s.ok:
                        src[route] = s.items()
            if not ref or len(src) < 2:
                continue
            selections = [("all", lambda p: True)]
            if motion == "free":
                selections.append(("excl", lambda p: not any(EXCLUDE[0] - 1e-9 <= w <= EXCLUDE[1] + 1e-9 for w in p)))
            for (label, f_old, key), (_, f_new, _) in zip(src[today], src[new]):
                for sel, kept in selections:
                    for q, index in QUANTITIES:
                        exact = {p: v[q] for p, v in ref[key].items()}
                        a, b = against_exact(f_old, exact, kept, index), against_exact(f_new, exact, kept, index)
                        if a is None or b is None:
                            problems.append("%s %s %s %s: no pair in common with the reference" % (case, kernel, label, q))
                            continue
                        both = [p for p in exact if f_old(p[0], p[1], index) is not None
                                and f_new(p[0], p[1], index) is not None]
                        scale = max(abs(exact[p]) for p in both)
                        chg = stats([abs(f_new(p[0], p[1], index) - f_old(p[0], p[1], index)) / scale
                                     for p in both if kept(p)]) or (float("nan"),) * 2
                        ok_max = b[0] <= a[0] + margin(a[0])
                        ok_med = b[1] <= a[1] + margin(a[1])
                        verdict.append(ok_max and ok_med)
                        lines.append("%-27s %-4s %-4s %-10s %-3s %8.2f %8.2f %8.2f %8.2f %8.3f %8.3f  %s / %s" % (
                            case, kernel, sel, label, q, a[0], a[1], b[0], b[1], chg[0], chg[1],
                            "pass" if ok_max else "FAIL", "pass" if ok_med else "FAIL"))
    lines.append("G1 over the rows above: %s (%d rows, %d fail)" % ("PASS" if verdict and all(verdict) else "FAIL",
                                                                   len(verdict), verdict.count(False)))

    # 2. in-run direct against indirect, every case (the round hull has no exact reference)
    lines += ["", "2. In-run check |T_dir - T_pot| / max |T_pot| over the rows of forces / of moments, all modes and "
              "heading pairs, max / median, percent"]
    for case in cases:
        for route in (today, new):
            f = folders[(route, case)]
            if not os.path.isdir(f):
                continue
            for kernel in ("diff", "sum"):
                if kernel == "diff":
                    t_dir = qc.read_part_sums(os.path.join(f, "QTF_DirectTerms.txt"), 3)
                    t_pot = qc.read_part_sums(os.path.join(f, "QTF_PotentialTerms.txt"), 7)
                else:
                    terms = qc.read_direct_sum_terms(os.path.join(f, "QTF_DirectTermsSum.txt"))
                    t_dir = {k: sum(v[:3]) for k, v in terms.items()}
                    t_pot = {k: sum(v[3:]) for k, v in terms.items()}
                if not t_dir or not t_pot:
                    problems.append("%s / %s: no %s terms for the in-run check" % (route, case, kernel))
                    continue
                try:
                    cells = qc.direct_against_indirect(t_dir, t_pot)
                except ValueError as err:
                    problems.append("%s / %s %s in-run check: %s" % (route, case, kernel, err))
                    continue
                lines.append("%-27s %-8s %-4s forces %s  moments %s" % ((case, route, kernel) + tuple(
                    "%7.2f %7.2f" % (100 * c[0], 100 * c[1]) if c else "      -       -" for c in cells)))

    # 3. every output file, new against today
    summaries = {}
    for case in cases:
        a, b = folders[(new, case)], folders[(today, case)]
        if os.path.isdir(a) and os.path.isdir(b):
            text, summary, file_problems = qfc.compare_folders(a, b, "3. %s, %s against %s" % (case, new, today))
            summaries[case] = summary
            problems += ["%s: %s" % (case, p) for p in file_problems]
            lines += ["", text]

    # 4. time
    lines += ["", "4. Wall time and core-hours (24 CPUs billed)",
              "%-34s %-8s %10s %8s %5s" % ("case", "route", "wall s", "core-h", "code")]
    total, complete = 0.0, True
    for route in (today, new):
        for case in cases:
            path = os.path.join(args.run, "g1_" + route, case, "result.json")
            if not os.path.exists(path):
                problems.append("%s / %s: no result.json" % (route, case))
                complete = False
                continue
            r = json.load(open(path))
            if r.get("code", 1) != 0:
                problems.append("%s / %s: result.json code %s" % (route, case, r.get("code")))
            ch = r["wall_s"] * 24 / 3600.0
            total += ch
            lines.append("%-34s %-8s %10d %8.1f %5d" % (case, route, r["wall_s"], ch, r.get("code", -1)))
    lines.append("total %.1f core-h%s" % (total, "" if complete else " (INCOMPLETE: result.json missing)"))
    lines += ["", "PROBLEMS: %d" % len(problems)] + ["  " + p for p in problems]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "g1_compare.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote", os.path.join(out, "g1_compare.txt"))
    figure(out, summaries)
    sys.exit(1 if problems else 0)


def figure(out, summaries):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    files = sorted({name for s in summaries.values() for name in s if s[name]})
    if not files:
        return
    colours = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    width = 0.8 / max(1, len(summaries))
    for j, (case, s) in enumerate(sorted(summaries.items())):
        vals = [max(s[name][0], 1e-16) if s.get(name) else float("nan") for name in files]
        ax.bar([i + j * width for i in range(len(files))], vals, width, label=case, color=colours[j % len(colours)])
    ax.set_yscale("log")
    ax.set_xticks([i + 0.4 - width / 2 for i in range(len(files))])
    ax.set_xticklabels(files, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("max |new - today| / max |today| (group of modes)")
    ax.set_title("G1 item 3: qtf_fs_circle = panels against samples, per output file", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "g1_compare.png"), dpi=150)
    print("wrote", os.path.join(out, "g1_compare.png"))


if __name__ == "__main__":
    main()
