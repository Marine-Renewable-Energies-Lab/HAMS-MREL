"""
Writes test/data/qtf_direct_field_six/expected_direct.txt from a QTF_Field.txt of the run of qtf_direct_field_six.

Usage:
    python -I extract_direct_field_expected.py QTF_Field.txt [expected_direct.txt]

The file holds the rows the direct solve covers (has_direct = 1) with the columns PER1 PER2 BETA1 BETA2 K IPT Re_dir Im_dir,
the tokens copied as the run wrote them. The columns are read by the names of the "# Columns:" line. The integration test
qtf_direct_field_six_values compares a new run with this file to a relative 1e-6. Regenerate it only when the direct part
is meant to change: run the integration tests, then this script on
<build>/output_qtf_direct_field_six/Wamit_format/QTF_Field.txt.
"""

import os
import sys

KEYS = ("PER1", "PER2", "BETA1", "BETA2", "K", "IPT", "Re_dir", "Im_dir")


def main():
    source = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    target = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "..", "test", "data", "qtf_direct_field_six",
                                                                  "expected_direct.txt")
    names = None
    rows = []
    with open(source) as f:
        for line in f:
            if line.startswith("# Columns:"):
                names = line.split(":", 1)[1].split()
            elif names and line.strip() and not line.startswith("#"):
                tokens = line.split()
                if int(tokens[names.index("has_direct")]) == 1:
                    rows.append([tokens[names.index(key)] for key in KEYS])
    assert rows, "no row with has_direct = 1"
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w") as f:
        f.write("! Direct part of the rows of qtf_direct_field_six that the direct solve covers, from QTF_Field.txt of a run of\n")
        f.write("! this build (extract_direct_field_expected.py). Columns: %s\n" % " ".join(KEYS))
        for row in rows:
            f.write(" ".join(row) + "\n")
    print("%d rows written to %s" % (len(rows), target))


if __name__ == "__main__":
    main()
