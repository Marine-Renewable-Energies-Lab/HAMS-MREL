"""Compare the output files of two HAMS-MREL runs of the same case, file by file (standard library only).

Usage: python3 qtf_files_compare.py <new Wamit_format folder> <reference Wamit_format folder> [--label NAME]

For every file present in both folders:
- WAMIT-format QTF files (QTF_*.10d/.10s/.11d/.11s/.12d/.12s): rows keyed by (omega_1, omega_2, beta_1, beta_2, mode),
  value Re + i Im (the last two columns). Per mode: |new - ref| / max |ref| over the rows of the forces (modes 1-3) or of
  the moments (modes 4-6) the mode belongs to, as qtf_compare.py --direct (so that a mode near zero for the heading, sway
  at heading 0, is not measured against itself), its max and median over the rows; the line "all" takes the worst mode.
- Term files (QTF_PotentialTerms*.txt, QTF_DirectTerms*.txt; 16 digits): rows keyed as above, then Re Im per term. Per term:
  |new - ref| / (max over the rows, terms and modes of the force or moment group of |ref|), max and median over the rows
  and modes.
- QTF_PsiCheck.txt: rows keyed by every column but the last two (tag, omega, ...), value Re + i Im. Per tag (PSI, EXC):
  |new - ref| / max |ref| of the tag, max and median.
- Every other file: whether it is byte-identical; if not, the number of differing lines and the largest relative
  difference of the numbers on them (relative to the largest number of the file).
Rows of one file that the other lacks are counted and reported. Measure as qtf_compare.py (normalise by the largest value of
the set so that rows passing through zero do not dominate).
"""
import argparse
import os
import re
import statistics

QTF_FILE = re.compile(r"^QTF_.*\.1[012][ds]$")
TERMS_FILE = re.compile(r"^QTF_(Potential|Direct)Terms.*\.txt$")


def numbers(line):
    out = []
    for tok in line.split():
        try:
            out.append(float(tok))
        except ValueError:
            return None
    return out


def key_of(vals, n=5):
    return tuple(round(v, 6) for v in vals[:n])


def read_rows(path, skip_tag=False):
    rows = {}
    dup = 0
    for line in open(path, errors="replace"):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if skip_tag:
            parts = line.split()
            try:
                vals = [float(t) for t in parts[1:]]
            except ValueError:
                continue
            key = (parts[0],) + tuple(round(v, 9) for v in vals[:-2])
            val = complex(vals[-2], vals[-1])
        else:
            vals = numbers(line)
            if vals is None or len(vals) < 7:
                continue
            key = key_of(vals)
            val = vals[5:]
        if key in rows:
            dup += 1
        rows[key] = val
    return rows, dup


def group(mode):
    """0 for a force (modes 1-3 of each body), 1 for a moment (4-6)."""
    return ((int(mode) - 1) % 6) // 3


def stats(errs):
    if not errs:
        return float("nan"), float("nan")
    return max(errs), statistics.median(errs)


def compare_qtf(new, ref):
    a, _ = read_rows(new)
    b, _ = read_rows(ref)
    common = sorted(set(a) & set(b))
    out = {}
    modes = sorted({k[4] for k in common})
    gscale = {}
    for k in common:
        g = group(k[4])
        gscale[g] = max(gscale.get(g, 0.0), abs(complex(b[k][-2], b[k][-1])))
    for m in modes:
        ks = [k for k in common if k[4] == m]
        za = [complex(a[k][-2], a[k][-1]) for k in ks]
        zb = [complex(b[k][-2], b[k][-1]) for k in ks]
        scale = gscale[group(m)] or 1.0
        out["mode %d" % m] = stats([abs(x - y) / scale for x, y in zip(za, zb)])
    if out:
        out["all"] = (max(v[0] for v in out.values()), max(v[1] for v in out.values()))
    return out, len(set(a) - set(b)), len(set(b) - set(a)), len(common)


def compare_terms(new, ref):
    a, _ = read_rows(new)
    b, _ = read_rows(ref)
    common = sorted(set(a) & set(b))
    out = {}
    if not common:
        return out, len(a), len(b), 0
    nterm = min(len(a[common[0]]), len(b[common[0]])) // 2
    gscale = {}
    for k in common:
        g = group(k[4])
        gscale[g] = max([gscale.get(g, 0.0)] + [abs(complex(b[k][2 * t], b[k][2 * t + 1])) for t in range(nterm)])
    for m in sorted({k[4] for k in common}):
        ks = [k for k in common if k[4] == m]
        scale = gscale[group(m)] or 1.0
        for t in range(nterm):
            errs = [abs(complex(a[k][2 * t], a[k][2 * t + 1]) - complex(b[k][2 * t], b[k][2 * t + 1])) / scale for k in ks]
            name = "term %d" % (t + 1)
            mx, md = stats(errs)
            old = out.get(name, (0.0, 0.0))
            out[name] = (max(old[0], mx), max(old[1], md))
    out["all"] = (max(v[0] for v in out.values()), max(v[1] for v in out.values()))
    return out, len(set(a) - set(b)), len(set(b) - set(a)), len(common)


def compare_psi(new, ref):
    a, _ = read_rows(new, skip_tag=True)
    b, _ = read_rows(ref, skip_tag=True)
    common = sorted(set(a) & set(b))
    out = {}
    for tag in sorted({k[0] for k in common}):
        ks = [k for k in common if k[0] == tag]
        scale = max(abs(b[k]) for k in ks) or 1.0
        out[tag] = stats([abs(a[k] - b[k]) / scale for k in ks])
    if out:
        out["all"] = (max(v[0] for v in out.values()), max(v[1] for v in out.values()))
    return out, len(set(a) - set(b)), len(set(b) - set(a)), len(common)


def compare_other(new, ref):
    ta = open(new, "rb").read()
    tb = open(ref, "rb").read()
    if ta == tb:
        return "byte-identical"
    la = ta.decode(errors="replace").splitlines()
    lb = tb.decode(errors="replace").splitlines()
    ndiff = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
    big = 0.0
    for line in lb:
        v = numbers(line) or []
        big = max([big] + [abs(x) for x in v])
    worst = 0.0
    for x, y in zip(la, lb):
        if x == y:
            continue
        va, vb = numbers(x), numbers(y)
        if va and vb and len(va) == len(vb):
            worst = max([worst] + [abs(p - q) for p, q in zip(va, vb)])
    rel = worst / big if big else float("nan")
    return "%d lines differ, largest difference %.3e of the largest number" % (ndiff, rel)


def compare_folders(new, ref, label=""):
    lines = []
    names = sorted(set(os.listdir(new)) & set(os.listdir(ref)))
    only = sorted(set(os.listdir(new)) ^ set(os.listdir(ref)))
    head = "%s: %s against %s" % (label, new, ref) if label else "%s against %s" % (new, ref)
    lines.append(head)
    lines.append("|new - ref| / max |ref| (per mode, term or tag; 'all' the worst), max and median over the rows")
    summary = {}
    for name in names:
        p, q = os.path.join(new, name), os.path.join(ref, name)
        if QTF_FILE.match(name):
            res, extra, missing, n = compare_qtf(p, q)
        elif TERMS_FILE.match(name):
            res, extra, missing, n = compare_terms(p, q)
        elif name == "QTF_PsiCheck.txt":
            res, extra, missing, n = compare_psi(p, q)
        else:
            lines.append("  %-28s %s" % (name, compare_other(p, q)))
            continue
        summary[name] = res.get("all")
        lines.append("  %-28s %d rows in common, %d only in new, %d only in ref" % (name, n, extra, missing))
        for key, (mx, md) in res.items():
            lines.append("      %-10s max %.3e  median %.3e" % (key, mx, md))
    if only:
        lines.append("  in one folder only: %s" % ", ".join(only))
    return "\n".join(lines), summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("new")
    parser.add_argument("ref")
    parser.add_argument("--label", default="")
    args = parser.parse_args()
    text, _ = compare_folders(args.new, args.ref, args.label)
    print(text)


if __name__ == "__main__":
    main()
