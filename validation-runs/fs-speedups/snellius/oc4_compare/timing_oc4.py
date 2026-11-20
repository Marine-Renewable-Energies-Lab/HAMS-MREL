"""Wall time and core-hours of the OC4 WAMIT-model runs side by side: old HAMS, new HAMS, Nemoh.

Usage:
    python -I timing_oc4.py --hams-old DIR JOBID --hams-new DIR JOBID --nemoh DIR JOBID --sacct SACCT.txt --out DIR

SACCT.txt is the output of
    sacct -j <ids> -X -P --format=JobID,ElapsedRaw,NCPUS,CPUTimeRAW
(pipe-separated). The billed core-hours are CPUTimeRAW / 3600, i.e. 24 CPUs times the wall time.

Split of the time:
- HAMS: from the per-pass lines of log.txt.
  - free-surface part: "QTF free surface fields" + "QTF free surface psi";
  - body part of the second order: "QTF B2 control surface fields" + "QTF B2 control surface grad psi" (the fields on the
    control surfaces C and F), the "QTF sum indirect check ... in all" line and the "Mean drift source solve" line;
  - first order and the rest: the wall time minus the two (first-order radiation-diffraction solves, assembly, output).
- Nemoh: from the stage times of result.json (run_nemoh.sbatch). preProc + solver + postProc is the first order;
  QTFpreProc + QTFsolver + QTFpostProc is the second order. Nemoh's logs do not split the second order into body and
  free-surface parts, so it is one block.
Writes timing_oc4.txt and timing_oc4.png in --out. Exits 1 on any missing input.
"""
import argparse
import json
import re
import sys
from pathlib import Path

PASS = re.compile(r"^\s*(QTF free surface fields|QTF free surface psi|QTF B2 control surface fields|"
                  r"QTF B2 control surface grad psi):.*?,\s+([0-9.]+) s for")
CHECK_ALL = re.compile(r"QTF sum indirect check:\s+([0-9.]+) s in all")
SOURCE = re.compile(r"Mean drift source solve:\s+([0-9.]+) s to assemble.*?,\s+([0-9.]+) s to solve")


def read_sacct(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        f = line.strip().split("|")
        if len(f) >= 4 and f[0].isdigit():
            rows[f[0]] = {"wall_s": int(f[1]), "ncpus": int(f[2]), "core_h": int(f[3]) / 3600.0}
    return rows


def hams_split(run):
    log = Path(run, "log.txt")
    if not log.is_file():
        raise FileNotFoundError(log)
    fs = body = 0.0
    counts = {}
    for line in log.read_text(errors="replace").splitlines():
        m = PASS.match(line)
        if m:
            t = float(m.group(2))
            counts[m.group(1)] = counts.get(m.group(1), 0) + 1
            if m.group(1).startswith("QTF free surface"):
                fs += t
            else:
                body += t
            continue
        m = CHECK_ALL.search(line)
        if m:
            body += float(m.group(1))
            continue
        m = SOURCE.search(line)
        if m:
            body += float(m.group(1)) + float(m.group(2))
    if not counts:
        raise ValueError(f"no pass lines in {log}")
    if "completes successfully" not in log.read_text(errors="replace"):
        raise ValueError(f"{log} has no 'completes successfully' line")
    return {"fs": fs, "body2": body, "counts": counts}


def nemoh_split(run):
    st = json.loads(Path(run, "result.json").read_text())["stages"]
    first = st["preProc"] + st["solver"] + st["postProc"]
    second = st["QTFpreProc"] + st["QTFsolver"] + st["QTFpostProc"]
    return {"first": first, "second": second, "stages": st}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hams-old", nargs=2, required=True, metavar=("DIR", "JOBID"))
    ap.add_argument("--hams-new", nargs=2, required=True, metavar=("DIR", "JOBID"))
    ap.add_argument("--nemoh", nargs=2, required=True, metavar=("DIR", "JOBID"))
    ap.add_argument("--sacct", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    problems = []
    sacct = read_sacct(a.sacct)
    runs = []
    for label, (d, jid) in (("HAMS old (stencil + samples)", a.hams_old),
                            ("HAMS new (all speed-ups, circle = panels)", a.hams_new),
                            ("Nemoh v3.0.3", a.nemoh)):
        if jid not in sacct:
            problems.append(f"job {jid} not in {a.sacct}")
            continue
        r = {"label": label, "jobid": jid, **sacct[jid]}
        try:
            if label.startswith("Nemoh"):
                s = nemoh_split(d)
                r["parts"] = [("first order (HAMS: and the rest)", s["first"]), ("second order (body + FS, not split in the logs)",
                                                             s["second"])]
                r["note"] = "stages " + ", ".join(f"{k} {v} s" for k, v in s["stages"].items())
            else:
                s = hams_split(d)
                rest = r["wall_s"] - s["fs"] - s["body2"]
                r["parts"] = [("first order (HAMS: and the rest)", rest), ("second order, body part", s["body2"]),
                              ("second order, free-surface part", s["fs"])]
                r["note"] = "passes " + ", ".join(f"{k} x{v}" for k, v in s["counts"].items())
        except (OSError, ValueError, KeyError) as e:
            problems.append(f"{label}: {e}")
            continue
        runs.append(r)
    if problems:
        print("PROBLEMS:\n  " + "\n  ".join(problems))
        sys.exit(1)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    lines = ["OC4 WAMIT model (5022 hull panels), fat_genoa, 24 CPUs billed per job", "",
             f"{'run':44s} {'job':>9s} {'wall h':>7s} {'core-h':>7s}  split (h)"]
    for r in runs:
        split = "; ".join(f"{n} {t / 3600:.2f}" for n, t in r["parts"])
        lines.append(f"{r['label']:44s} {r['jobid']:>9s} {r['wall_s'] / 3600:7.2f} {r['core_h']:7.1f}  {split}")
        lines.append(f"{'':44s} {'':>9s} {'':>7s} {'':>7s}  ({r['note']})")
    old, new = runs[0], runs[1]
    lines += ["", f"new / old HAMS: wall {new['wall_s'] / old['wall_s']:.3f}, core-h {new['core_h'] / old['core_h']:.3f}",
              f"new HAMS / Nemoh: wall {new['wall_s'] / runs[2]['wall_s']:.3f}, core-h "
              f"{new['core_h'] / runs[2]['core_h']:.3f}",
              "Nemoh ran its stages on 1 thread; the job was billed for 24 CPUs (fat_genoa minimum)."]
    txt = "\n".join(lines)
    (out / "timing_oc4.txt").write_text(txt + "\n")
    print(txt)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"first order (HAMS: and the rest)": "#9aa5b1",
              "second order, body part": "#2a6fdb", "second order, free-surface part": "#e0892b",
              "second order (body + FS, not split in the logs)": "#7a5bc7"}
    fig, ax = plt.subplots(figsize=(9, 4.4))
    names = [f"{r['label']}\njob {r['jobid']}" for r in runs][::-1]
    seen = set()
    for i, r in enumerate(runs[::-1]):
        left = 0.0
        for n, t in r["parts"]:
            h = t / 3600
            ax.barh(i, h, left=left, color=colors[n], edgecolor="white", linewidth=2,
                    label=None if n in seen else n)
            seen.add(n)
            left += h
        ax.text(left + 0.08, i, f"{r['wall_s'] / 3600:.2f} h wall, {r['core_h']:.1f} core-h", va="center", fontsize=9)
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("wall time (h)")
    ax.set_xlim(0, max(r["wall_s"] for r in runs) / 3600 * 1.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#ddd", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, frameon=False)
    ax.set_title("OC4 WAMIT model: wall time and core-hours (24 CPUs billed per job)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "timing_oc4.png", dpi=150)


if __name__ == "__main__":
    main()
