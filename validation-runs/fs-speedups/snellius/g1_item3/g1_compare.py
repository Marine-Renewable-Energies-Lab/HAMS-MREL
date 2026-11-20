"""G1 item 3 after the Snellius runs: qtf_fs_circle = panels against samples, and both against the exact references.

Usage: python g1_compare.py <run folder> <repo> --sum-ref <qtfdirectsum_reference.csv> --diff-ref <qtffs_reference.csv>
                            [--out <folder>]

<run folder> holds g1_samples/<case> and g1_panels/<case> (run_one.sh). Writes to --out (default: the run folder):
- g1_compare.txt:
  1. per route, qtf_compare.py --direct-sum (the sum QTF .12s, .12s - .10s, P and D_FS against the exact sharp cylinder,
     and the in-run direct-indirect check, round hull included), with all pairs and with 0.75-0.80 rad/s left out (the
     free-body resonance, as compare_stage5_excl080.txt of 7c-2);
  2. per route, qtf_compare.py --direct (the difference QTF .12d and .11d against the exact sharp cylinder);
  3. per case, every output file of panels against samples (qtf_files_compare.py);
  4. wall time and core-hours of each run (result.json; 24 CPUs billed).
- g1_compare.png (if matplotlib is there): the largest panels-against-samples difference per QTF file and case, log axis.
qtf_compare.py is only called, not edited. The G1 rule (spec G1): no quantity worse than samples by more than the larger of
0.5 percentage points and 10 % of the samples error, against the exact reference.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import qtf_files_compare as qfc  # noqa: E402

ROUTES = ("samples", "panels")


def run(cmd):
    res = subprocess.run(cmd, capture_output=True, text=True)
    return "$ qtf_compare.py " + " ".join(os.path.basename(c) for c in cmd[2:]) + "\n" + res.stdout + res.stderr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("repo")
    parser.add_argument("--sum-ref", required=True)
    parser.add_argument("--diff-ref", required=True)
    parser.add_argument("--out")
    args = parser.parse_args()
    out = args.out or args.run
    os.makedirs(out, exist_ok=True)
    qc = os.path.join(args.repo, "pre_postprocessing_python", "SecondOrder", "qtf_compare.py")
    py = sys.executable
    parts = []
    for route in ROUTES:
        root = os.path.join(args.run, "g1_" + route)
        parts.append("=== %s: sum, all pairs\n%s" % (route, run([py, qc, root, "--direct-sum", args.sum_ref])))
        parts.append("=== %s: sum, 0.75-0.80 left out\n%s" % (route, run([py, qc, root, "--direct-sum", args.sum_ref,
                                                                          "--exclude", "0.75", "0.80"])))
        parts.append("=== %s: difference\n%s" % (route, run([py, qc, root, "--direct", args.diff_ref])))
    cases = sorted(set(os.listdir(os.path.join(args.run, "g1_panels"))) & set(os.listdir(os.path.join(args.run, "g1_samples"))))
    summaries = {}
    for case in cases:
        new = os.path.join(args.run, "g1_panels", case, "output", "Wamit_format")
        ref = os.path.join(args.run, "g1_samples", case, "output", "Wamit_format")
        if not (os.path.isdir(new) and os.path.isdir(ref)):
            parts.append("=== %s: output missing" % case)
            continue
        text, summary = qfc.compare_folders(new, ref, "%s, panels against samples" % case)
        summaries[case] = summary
        parts.append("=== " + text)
    lines = ["=== wall time and core-hours (24 CPUs billed)", "%-34s %-8s %10s %8s %5s" % ("case", "route", "wall s", "core-h",
                                                                                        "code")]
    total = 0.0
    for route in ROUTES:
        for case in cases:
            path = os.path.join(args.run, "g1_" + route, case, "result.json")
            if os.path.exists(path):
                r = json.load(open(path))
                ch = r["wall_s"] * 24 / 3600.0
                total += ch
                lines.append("%-34s %-8s %10d %8.1f %5d" % (case, route, r["wall_s"], ch, r["code"]))
    lines.append("total %.1f core-h" % total)
    parts.append("\n".join(lines))
    with open(os.path.join(out, "g1_compare.txt"), "w") as f:
        f.write("\n\n".join(parts) + "\n")
    print("wrote", os.path.join(out, "g1_compare.txt"))
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
    ax.set_ylabel("max |panels - samples| / max |samples| (group of modes)")
    ax.set_title("G1 item 3: qtf_fs_circle = panels against samples, per output file", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "g1_compare.png"), dpi=150)
    print("wrote", os.path.join(out, "g1_compare.png"))


if __name__ == "__main__":
    main()
