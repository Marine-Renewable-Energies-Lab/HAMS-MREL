"""Lid-mesh convergence study of the dissipative damping lid against the eigenfunction reference.

The case is that of the mm_cyl_* fixtures (make_lid_fixtures.py): a cylinder of radius 2 m standing on the sea bed in 10 m of
water, k = 0.15 (kh = 1.5), heading 0, the lid r 2 to 4 on z = 0. The hull is the 48 x 20 wall, not the 24 x 10 one of the
fixtures, since the 24 x 10 hull sets the error floor of the lid effect (test/data/damping_lid/README.md). The lid has 48 panels
around and NR rings. Field points: (0, 8, 0), (8, 0, 0), then the lid centroids, for the lid power.

The cases are written to WORK (not to test/data), run with the given hamsmrel executable, and compared with the reference
values of annular_lid_reference.py (README, "Reference values"). Prints a Markdown table. Not run by the test suite.

    python -I lid_convergence_study.py --exe ../../install/bin/hamsmrel.exe --work C:/scratch/study [--rings 8 16 32]
        [--eps 0.2 1] [--around 48] [--threads 4] [--time-only]
"""
import argparse
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_lid_fixtures as fx  # noqa: E402

OMEGA = 1.1538937638
RHO = 1000.0
# annular_lid_reference.py, lid projection, extrapolated in 1/N from N = 160, 320 (README, "Reference values"): X1 over
# rho g A, the elevation at (0, 8, 0) and at (8, 0, 0) over A, and the lid power in W per m**2 of incident amplitude.
MACCAMY_FUCHS = [complex(1.65550518, 23.42095485), complex(1.01387576, 0.04392707), complex(0.29600397, -0.97000684)]
REFERENCE = {
    0.2: ([complex(1.73472971, 23.39194867), complex(1.00952678, 0.04447564), complex(0.29053829, -0.96832198)], 5268.414),
    1.0: ([complex(2.02520747, 23.23109427), complex(0.99339102, 0.04910955), complex(0.27074161, -0.95859350)], 24729.64),
}
FAR = [(0.0, 8.0, 0.0), (8.0, 0.0, 0.0)]


def lid_areas(lid_text):
    """Panel areas of a quad lid, half the cross product of the diagonals, as read_lid_areas in the test."""
    lines = lid_text.split(fx.EOL)
    k = lines.index("    # Start Definition of Node Coordinates     ! node_number   x   y   z") + 1
    nodes = {}
    while not lines[k].startswith("    # End"):
        f = lines[k].split()
        nodes[int(f[0])] = [float(v) for v in f[1:4]]
        k += 1
    k += 3
    areas = []
    while not lines[k].startswith("    # End"):
        f = [int(v) for v in lines[k].split()]
        p = [nodes[i] for i in f[2:2 + f[1]]]
        if len(p) == 3:
            u = [p[1][c] - p[0][c] for c in range(3)]
            v = [p[2][c] - p[0][c] for c in range(3)]
        else:
            u = [p[2][c] - p[0][c] for c in range(3)]
            v = [p[3][c] - p[1][c] for c in range(3)]
        areas.append(0.5 * abs(u[0] * v[1] - u[1] * v[0]))
        k += 1
    return areas


def write(work, name, lid_text, hull, points, threads):
    fx.TARGET = work
    fx.write_case(name, lid_text, lid_flag=0 if lid_text is None else 1, points=points, depth=10.0, wave_number=0.15,
                  hull_files=hull)
    path = os.path.join(work, name, "input", "ControlFile.in")
    with open(path, "rb") as f:
        text = f.read().decode("ascii")
    text = text.replace("Number of threads       1", "Number of threads       %d" % threads)
    with open(path, "wb") as f:
        f.write(text.encode("ascii"))


def run(exe, work, name):
    out = os.path.join(work, name, "output")
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(work, name + ".log"), "w") as log:
        code = subprocess.call([exe, os.path.join(work, name, "input") + os.sep, out + os.sep], stdout=log,
                               stderr=subprocess.STDOUT)
    if code != 0:
        raise RuntimeError("%s: exit code %d" % (name, code))
    return time.time() - t0


def results(work, name):
    """X1 and the diffraction elevations at all field points, in the order of the field points."""
    out = os.path.join(work, name, "output", "Wamit_format")
    x1 = None
    with open(os.path.join(out, "ExcForce.3")) as f:
        for line in f:
            v = line.split()
            if len(v) == 7 and int(v[2]) == 1 and abs(float(v[1])) < 1e-6:
                x1 = complex(float(v[5]), float(v[6]))
    eta = {}
    with open(os.path.join(out, "PressureElevation.6p")) as f:
        for line in f:
            v = line.split()
            if len(v) == 7:
                eta[int(float(v[2]))] = complex(float(v[5]), float(v[6]))
    return x1, [eta[i] for i in sorted(eta)]


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exe", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--rings", type=int, nargs="+", default=[8, 16, 32])
    p.add_argument("--eps", type=float, nargs="+", default=[0.2, 1.0])
    p.add_argument("--around", type=int, default=48)
    p.add_argument("--hull", type=int, nargs=2, default=[48, 20], help="panels around and rows of the wall")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--time-only", action="store_true", help="run the finest lid at the first epsilon only, and stop")
    args = p.parse_args()
    exe = os.path.abspath(args.exe)
    work = os.path.abspath(args.work)
    hull = (fx.wall(2.0, 10.0, args.hull[0], args.hull[1]), fx.wall_hydrostatics(2.0, 10.0))
    tag = "h%dx%d" % tuple(args.hull)

    if args.time_only:
        nr, eps = max(args.rings), args.eps[0]
        lid = fx.annulus(2.0, 4.0, args.around, nr, dampers=[(1, args.around * nr, eps)])
        name = "study_%s_%dx%d_eps%g" % (tag, args.around, nr, eps)
        write(work, name, lid, hull, FAR + fx.centroids(lid, args.around * nr), args.threads)
        print("%s: %d hull + %d lid panels, %.1f s" % (name, args.hull[0] * args.hull[1], args.around * nr,
                                                     run(exe, work, name)))
        return

    bare_name = "study_%s_nolid" % tag
    write(work, bare_name, None, hull, FAR, args.threads)
    t = run(exe, work, bare_name)
    bx, beta = results(work, bare_name)
    bare = [bx, beta[0], beta[1]]
    print("no lid, %s: %.1f s; against MacCamy-Fuchs %s" % (tag, t, ", ".join("%.2e" % rel(a, b)
                                                                            for a, b in zip(bare, MACCAMY_FUCHS))))
    print()
    print("| epsilon (m/s) | Lid | Lid panels | X1 | eta(0, 8, 0) | eta(8, 0, 0) | Lid effect: X1, eta(0, 8, 0), eta(8, 0, 0) "
          "| Lid power | Run (s) |")
    print("| ---: | :--- | ---: | ---: | ---: | ---: | :--- | ---: | ---: |")
    for eps in args.eps:
        ref, ref_power = REFERENCE[eps]
        for nr in args.rings:
            n = args.around * nr
            lid = fx.annulus(2.0, 4.0, args.around, nr, dampers=[(1, n, eps)])
            name = "study_%s_%dx%d_eps%g" % (tag, args.around, nr, eps)
            write(work, name, lid, hull, FAR + fx.centroids(lid, n), args.threads)
            t = run(exe, work, name)
            x1, eta = results(work, name)
            v = [x1, eta[0], eta[1]]
            raw = [rel(a, b) for a, b in zip(v, ref)]
            effect = [abs((a - c) - (b - m)) / abs(b - m) for a, b, c, m in zip(v, ref, bare, MACCAMY_FUCHS)]
            area = lid_areas(lid)
            power = 0.5 * RHO * OMEGA ** 2 * sum(eps * abs(e) ** 2 * s for e, s in zip(eta[2:2 + n], area))
            print("| %g | %d x %d | %d | %.2e | %.2e | %.2e | %s | %+.2e | %.0f |" % (
                eps, args.around, nr, n, raw[0], raw[1], raw[2], ", ".join("%.2e" % e for e in effect),
                (power - ref_power) / ref_power, t))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
