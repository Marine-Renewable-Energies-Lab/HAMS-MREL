"""The 200 m cut-off check of the VolturnUS-S off-diagonal gap (docs/second_order_validation.md, section 8).

The case of make_volturnus_case.py (--fmax 2.50 --lid) with the inner free surface to R = 200 m, WAMIT's RINNER, in 10 m
panels, indirect only (no qtf_direct) and with qtf_fs_circle = panels. Only FS_outer of this run is used: the free-surface
integral beyond 200 m, whose partial waves come from the hull integral and do not depend on the inner mesh.
QTF_Indirect.11d of q_volturnus_free_full_ext minus it is HAMS with the free surface cut at 200 m (volturnus_cut200.py).
Usage: python -I make_volturnus_r200.py IEA_HYDRODATA_DIR OUT_CASE_DIR
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_volturnus_case as mvc  # noqa: E402

if len(sys.argv) != 3:
    sys.exit(__doc__)
mvc.FS_RADIUS = 200.0
mvc.main(["make_volturnus_case.py", "--fmax", "2.50", "--fs-panel", "10", "--lid", sys.argv[1], sys.argv[2]])
inp = os.path.join(sys.argv[2], "input", "SecondOrder.in")
lines = [line for line in open(inp).read().splitlines() if line.strip() != "qtf_direct = 1"]
lines += ["qtf_fs_circle = panels"]
open(inp, "w", newline="\n").write("\n".join(lines) + "\n")
