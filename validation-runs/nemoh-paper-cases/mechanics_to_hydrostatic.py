"""Write HAMS-MREL's Hydrostatic.in from the Mechanics folder that Nemoh v3.0.3's hydrosCal writes (Inertia.dat, Kh.dat,
and Km.dat if present) and the centre of gravity of Mesh.cal (its fourth line).

Used on the cluster for the Nemoh test cases 8b (cylinder) and 11 (OC3 Hywind), which ship without Mechanics: their mass
and inertia are Nemoh's own estimate, and the HydroStar reference files of those cases are compared on that basis.
Nemoh uses rho = 1025 kg/m^3 and g = 9.81 m/s^2, HAMS rho = 1000 and g = 9.80665: the mass and inertia are scaled by
1000/1025, the stiffness by (1000 x 9.80665)/(1025 x 9.81), as make_cases.py does. Plain Python 3, no packages.

Usage: python3 mechanics_to_hydrostatic.py NEMOH_CASE_DIR OUT_HYDROSTATIC_IN
"""
import os
import sys

MASS_SCALE = 1000.0 / 1025.0
STIFF_SCALE = 1000.0 * 9.80665 / (1025.0 * 9.81)


def matrix(path):
    with open(path) as f:
        rows = [[float(v) for v in ln.split()] for ln in f if ln.strip()]
    if len(rows) != 6 or any(len(r) != 6 for r in rows):
        raise SystemExit("%s is not a 6 x 6 matrix" % path)
    return rows


def main(case, out):
    mech = os.path.join(case, "Mechanics")
    with open(os.path.join(case, "Mesh.cal")) as f:
        cg = [float(v) for v in f.read().split("\n")[3].split()[:3]]
    zero = [[0.0] * 6 for _ in range(6)]
    mass = [[v * MASS_SCALE for v in r] for r in matrix(os.path.join(mech, "Inertia.dat"))]
    kh = [[v * STIFF_SCALE for v in r] for r in matrix(os.path.join(mech, "Kh.dat"))]
    km_path = os.path.join(mech, "Km.dat")
    km = [[v * STIFF_SCALE for v in r] for r in matrix(km_path)] if os.path.exists(km_path) else zero
    if mass[0][0] <= 0.0 or kh[2][2] <= 0.0:
        raise SystemExit("hydrosCal gave a non-positive mass or heave stiffness")
    with open(out, "w") as f:
        f.write(" Center of Gravity:\n  %.12E  %.12E  %.12E\n" % tuple(cg))
        for title, m in (("Body Mass Matrix", mass), ("External Linear Damping Matrix", zero),
                         ("External Quadratic Damping Matrix", zero), ("Hydrostatic Restoring Matrix", kh),
                         ("External Restoring Matrix", km)):
            f.write(" %s:\n" % title)
            for row in m:
                f.write("".join("  %12.5E" % v for v in row) + "\n")
    print("Hydrostatic.in from %s: CG %s, mass %.6E kg (HAMS density), C33 %.6E, C55 %.6E N m"
          % (mech, cg, mass[0][0], kh[2][2], kh[4][4]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
