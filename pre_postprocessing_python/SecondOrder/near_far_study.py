"""
Near-field against far-field mean drift on the hull of drift_cyl as the mesh is refined, held fixed, free, and free on
external springs, against the implementation plan's 5% limit on F_x. Not part of CI: a run takes a minute at 1792 panels
and many minutes beyond. See test/data/drift_cyl/README.md.

Usage:
    python near_far_study.py <path to hamsmrel> <work directory> [--bilge-radius R] [--multiples 1,2,3]

With no --bilge-radius, the hull has the sharp bilge of drift_cyl. With one, it is the hull of drift_round, with the wall
rows of each fixture: cosine when fixed and free, uniform in z on the springs. --rows cosine or uniform uses those on all
three. Each multiple scales all the panel counts. The exit code is 1 if F_x misses the limit at the finest mesh, at frequencies away from
resonance for the moving hull. Prints the terms of F_x, for finding where a difference comes from.
"""

import argparse
import os
import subprocess
import sys

import make_drift_fixtures as fixtures
import make_mesh

LIMIT = 0.05
FREQUENCIES = fixtures.DRIFT_CYL_FREQUENCIES
# Case name, motion_state, external restoring, frequencies below this are printed but not checked, as in
# test_integration_drift.f90: near the pitch resonance, and near the heave resonance on the springs
CASES = (("fixed", "fixed", None, 0.0), ("free", "free", None, 0.8), ("springs", "free", fixtures.EXTERNAL_RESTORING, 1.0))


def write_hull(path, multiple, bilge_radius, external, uniform_rows):
    if bilge_radius is None:
        mesh = dict(fixtures.DRIFT_CYL_MESH, external_restoring=external)
        for key in ("ntheta", "nz", "nr"):
            mesh[key] *= multiple
        return make_mesh.write_cylinder(os.path.join(path, "HullMesh.pnl"), os.path.join(path, "Hydrostatic.in"), **mesh)
    mesh = dict(fixtures.ROUND_MESH, bilge_radius=bilge_radius, external_restoring=external, uniform_rows=uniform_rows)
    for key in ("ntheta", "nz", "narc", "nr"):
        mesh[key] *= multiple
    return make_mesh.write_rounded_cylinder(os.path.join(path, "HullMesh.pnl"), os.path.join(path, "Hydrostatic.in"), **mesh)


def read_fx(path):
    """F_x per frequency from a mean drift file."""
    values = {}
    for line in open(path):
        columns = line.split()
        if int(columns[3]) == 1:
            values[round(float(columns[0]), 4)] = float(columns[6])
    return values


def read_terms(path):
    """The terms of F_x per frequency."""
    values = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        columns = line.split()
        if int(columns[4]) == 1:
            values.setdefault(round(float(columns[0]), 4), []).append(float(columns[6]))
    return values


def run_case(hams, work, name, motion, external, multiple, bilge_radius, rows):
    """Run one case. Returns the panel count and the rows (omega, far, near, terms)."""
    base = os.path.join(work, "%s_x%d" % (name, multiple))
    path = os.path.join(base, "input")
    os.makedirs(path, exist_ok=True)
    uniform_rows = rows == "uniform" or (rows == "fixture" and external is not None)
    mesh = write_hull(path, multiple, bilge_radius, external, uniform_rows)
    fixtures.write_control(os.path.join(path, "ControlFile.in"), frequencies=FREQUENCIES, threads=os.cpu_count())
    fixtures.write_second_order(os.path.join(path, "SecondOrder.in"), mean_drift=1, drift_methods="far near",
                                motion_state=motion, write_terms=1)
    output = os.path.join(base, "output")
    with open(os.path.join(base, "log.txt"), "w") as log:
        subprocess.run([hams, path + os.sep, output + os.sep], stdout=log, stderr=subprocess.STDOUT, check=True)
    folder = os.path.join(output, "Wamit_format")
    far = read_fx(os.path.join(folder, "MeanDriftMomentum.8"))
    near = read_fx(os.path.join(folder, "MeanDriftPressure.9"))
    terms = read_terms(os.path.join(folder, "MeanDriftTerms.txt"))
    return len(mesh.panels), [(omega, far[omega], near[omega], terms[omega]) for omega in sorted(far)]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("hams")
    parser.add_argument("work")
    parser.add_argument("--bilge-radius", type=float)
    parser.add_argument("--multiples", default="1,2,3")
    parser.add_argument("--rows", choices=("fixture", "cosine", "uniform"), default="fixture")
    args = parser.parse_args()
    multiples = [int(text) for text in args.multiples.split(",")]

    failed = False
    for name, motion, external, unchecked_below in CASES:
        for multiple in multiples:
            panels, rows = run_case(os.path.abspath(args.hams), args.work, name, motion, external, multiple, args.bilge_radius,
                                    args.rows)
            for omega, far, near, terms in rows:
                error = abs(near - far) / abs(far)
                checked = multiple == multiples[-1] and omega >= unchecked_below
                failed = failed or (checked and error >= LIMIT)
                print("%-7s %6d panels  omega %.3f  far %10.5f  near %10.5f  error %7.2f%%  terms %s" % (
                    name, panels, omega, far, near, 100.0 * error, " ".join("%.4f" % t for t in terms)), flush=True)
    print("5%% limit at the finest mesh: %s" % ("missed" if failed else "met"))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
