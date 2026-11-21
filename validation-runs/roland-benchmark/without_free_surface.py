"""QTF_Indirect.11d of a qtf_terms = full run without its free-surface terms, written in the same WAMIT format, so that
compare_diagonals.py and compare_maps.py can show HAMS with and without the free-surface integral.

.11d = .10d + the seven terms of QTF_PotentialTerms.txt (output_files.md); the last two are FS_inner and FS_outer. This
writes .11d - FS_inner - FS_outer row by row (the terms are scaled and conjugated as .11d, and keyed by omega_i, omega_j,
beta_1, beta_2 and the index). Diagonal rows have no terms and are copied. The row format is that of .11d
(4ES14.6, I6, 4ES14.6); the phase column is in degrees.

With --direct it does the same for QTF_Direct.12d: .12d - D_FS, the third part of QTF_DirectTerms.txt (phi_D of the
free-surface forcing; output_files.md).

Usage: python -I without_free_surface.py WAMIT_FORMAT_DIR OUT_FILE [--direct]
"""
import cmath
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "pre_postprocessing_python", "SecondOrder"))
import qtf_tools  # noqa: E402


def read_direct_terms(path):
    """{key: [D_I2, D_D, D_FS]} of QTF_DirectTerms.txt, keyed as qtf_tools.read_qtf_terms."""
    out = {}
    for line in open(path):
        if line.lstrip().startswith("#"):
            continue
        p = line.split()
        if len(p) != 11:
            continue
        key = (round(float(p[0]), 6), round(float(p[1]), 6), round(float(p[2]), 6), round(float(p[3]), 6), int(p[4]))
        out[key] = [complex(float(p[5 + 2 * t]), float(p[6 + 2 * t])) for t in range(3)]
    return out


def main(folder, out, direct=False):
    if direct:
        raw = read_direct_terms(os.path.join(folder, "QTF_DirectTerms.txt"))
        terms = {k: [0j] * 5 + [v[2], 0j] for k, v in raw.items()}   # D_FS in the FS_inner slot
        source = "QTF_Direct.12d"
    else:
        terms = qtf_tools.read_qtf_terms(os.path.join(folder, "QTF_PotentialTerms.txt"))
        source = "QTF_Indirect.11d"
    if not terms or len(next(iter(terms.values()))) != 7:
        raise ValueError("no free-surface terms in %s (qtf_terms = full needed)" % folder)
    rows, removed = [], 0
    for line in open(os.path.join(folder, source)):
        p = line.split()
        if len(p) != 9:
            continue
        key = (round(float(p[0]), 6), round(float(p[1]), 6), round(float(p[2]), 6), round(float(p[3]), 6), int(p[4]))
        v = complex(float(p[7]), float(p[8]))
        if key in terms:
            v -= terms[key][5] + terms[key][6]
            removed += 1
        elif float(p[0]) != float(p[1]):
            raise ValueError("no terms for the off-diagonal row %s" % (key,))
        rows.append("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E" % (
            float(p[0]), float(p[1]), float(p[2]), float(p[3]), int(p[4]), abs(v), math.degrees(cmath.phase(v)),
            v.real, v.imag))
    open(out, "w").write("\n".join(rows) + "\n")
    print("wrote %d rows (%d with the free-surface terms removed) to %s" % (len(rows), removed, out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--direct" in sys.argv[3:])
