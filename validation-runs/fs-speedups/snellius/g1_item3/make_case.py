"""Write one Phase 7c-2 case of drift_convergence.py for the G1 item-3 runs, with the circle route set.

Usage: python3 make_case.py <repo> <root> <case> <samples|panels> [--threads 24]

The case is written as `drift_convergence.py run` writes it (write_hams_case: hull, lid, ControlFile.in with the thread
count, SecondOrder.in with the case's options), then `qtf_fs_circle = <route>` is appended to SecondOrder.in. Every other
key stays at the branch default (field_deep_kh = 10, qtf_fs_inner = stencil). Prints the input folder. drift_convergence.py
is only imported, not edited (another branch edits it).
"""
import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", help="checkout of HAMS-MREL")
    parser.add_argument("root", help="folder of the runs of one route")
    parser.add_argument("case", help="a case of drift_convergence.py --qtf-potential, e.g. q_sharp_directsum_fixed_m2")
    parser.add_argument("route", choices=("samples", "panels"))
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    sys.path.insert(0, os.path.join(args.repo, "pre_postprocessing_python", "SecondOrder"))
    import drift_convergence as dc
    found = [case for case in dc.qtf_potential_cases() if case.name == args.case]
    if len(found) != 1:
        raise SystemExit("no case named %s in drift_convergence.py --qtf-potential" % args.case)
    path = dc.write_hams_case(found[0], args.root, args.threads)
    second = os.path.join(path, "SecondOrder.in")
    text = open(second).read()
    if "qtf_fs_circle" in text:
        raise SystemExit("%s already sets qtf_fs_circle" % second)
    with open(second, "a") as f:
        f.write("qtf_fs_circle = %s\n" % args.route)
    print(path)


if __name__ == "__main__":
    main()
