"""The OC4 WAMIT-model sum case with the fs-speedups items on, against job 29453291 (standard library; the figure needs
matplotlib and is skipped without it).

Usage: python3 compare_oc4.py <new case folder> <29453291 case folder> [--out <folder>]

A case folder holds log.txt, result.json, time.txt (if /usr/bin/time ran) and output/Wamit_format. Writes to --out (default:
the new case folder) compare_oc4.txt and, with matplotlib, compare_oc4.png:
1. every QTF file of the new run against 29453291 (qtf_files_compare.py: per mode, term or tag, max and median of
   |new - ref| / max |ref| over the rows of the force or moment group), and the other output files (byte-identical or not);
2. wall time and core-hours of both runs (result.json wall_s; 24 CPUs billed on fat_genoa), and the speed-up;
3. the time of each logged pass (the "QTF ...: omega ..., <t> s for ..." lines of log.txt), summed per kind, both runs.
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import qtf_files_compare as qfc  # noqa: E402

CPUS = 24
PASS = re.compile(r"^\s*(QTF [^:]+?):\s+omega\s+[-\d.E+]+,\s+([\d.]+) s for")


def wall(case):
    path = os.path.join(case, "result.json")
    if not os.path.exists(path):
        return None
    return float(json.load(open(path))["wall_s"])


def passes(case):
    out = {}
    path = os.path.join(case, "log.txt")
    if not os.path.exists(path):
        return out
    for line in open(path, errors="replace"):
        m = PASS.match(line)
        if m:
            n, t = out.get(m.group(1), (0, 0.0))
            out[m.group(1)] = (n + 1, t + float(m.group(2)))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("new")
    parser.add_argument("ref")
    parser.add_argument("--out")
    args = parser.parse_args()
    out = args.out or args.new
    text, summary = qfc.compare_folders(os.path.join(args.new, "output", "Wamit_format"),
                                        os.path.join(args.ref, "output", "Wamit_format"), "OC4 fs-speedups against 29453291")
    lines = [text, "", "Time (%d CPUs billed)" % CPUS]
    wn, wr = wall(args.new), wall(args.ref)
    for name, w in (("new", wn), ("29453291", wr)):
        if w is not None:
            lines.append("  %-9s wall %8.0f s = %6.2f h, %7.1f core-h" % (name, w, w / 3600.0, w * CPUS / 3600.0))
    if wn and wr:
        lines.append("  speed-up %.2f x, saved %.1f core-h" % (wr / wn, (wr - wn) * CPUS / 3600.0))
    pn, pr = passes(args.new), passes(args.ref)
    lines.append("")
    lines.append("Logged passes, count and total s: %-44s %18s %18s" % ("", "new", "29453291"))
    for kind in sorted(set(pn) | set(pr)):
        a, b = pn.get(kind, (0, 0.0)), pr.get(kind, (0, 0.0))
        lines.append("  %-76s %4d %10.0f s %4d %10.0f s" % (kind, a[0], a[1], b[0], b[1]))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "compare_oc4.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    files = [name for name in sorted(summary) if summary[name]]
    if not files:
        return
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = range(len(files))
    ax.bar([i - 0.2 for i in x], [max(summary[n][0], 1e-16) for n in files], 0.4, label="max", color="#0072B2")
    ax.bar([i + 0.2 for i in x], [max(summary[n][1], 1e-16) for n in files], 0.4, label="median", color="#E69F00")
    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(files, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("|new - 29453291| / max |29453291| (group of modes)")
    title = "OC4 WAMIT model, sum QTF: fs-speedups defaults + qtf_fs_circle = panels against job 29453291"
    if wn and wr:
        title += "\nwall %.2f h against %.2f h (%.1f against %.1f core-h)" % (wn / 3600, wr / 3600, wn * CPUS / 3600,
                                                                          wr * CPUS / 3600)
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "compare_oc4.png"), dpi=150)


if __name__ == "__main__":
    main()
