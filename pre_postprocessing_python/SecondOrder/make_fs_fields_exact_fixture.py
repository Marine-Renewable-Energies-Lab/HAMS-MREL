"""The fixture of test_qtf_free_surface's fs_fields_exact_* tests: HAMS's evaluation of the first-order fields on the
free surface (fs_inner_fields: phi, grad phi and phi_zz by the 5-point stencil; fs_circle_modes: the partial waves from
the circle samples) with the exact potential of the sharp truncated cylinder (truncated_cylinder_reference.py: radius
10, draft 10, depth 100 m, fixed body, heading 0) on the hull panels, no BEM solve, against the exact series.

    python make_fs_fields_exact_fixture.py [hull.pnl] [out]

Hull: the study's sharp m1 HullMesh.pnl (448 panels), copied to test/data/fs_sum_exact/HullMesh_sharp_m1.pnl; the
centroids are vertex averages, as read_pnl_hull (second_order_fixtures) makes them. On each panel: phi_P = phi - phi_I,
the exact scattered potential at the centroid; the test sets d phi_P/dn = -d phi_I/dn (fixed body), as
maccamy_fuchs_state. Every value is times C = i, HAMS's incident wave over the reference's (-i g/omega against -g/omega,
make_fs_sum_exact_fixture.py).

File: one block per frequency. 'FREQ omega k g h', 'PANELS n' then n lines 'Re Im' of phi_P at the panels in file order,
'POINTS n' then n lines 'x y  Re Im of phi, gx, gy, gz, phizz' (total fields on z = 0, near the waterline and out to
13.5 m), 'CIRCLE n' then n lines 'x y  Re Im of phi, gx, gy, gz, phizz' of the scattered part only at r = 15.5 and 20 m
(beyond the R = 15 m circle, where the partial waves hold), 'END'."""
import math, os, sys
import numpy as np
import truncated_cylinder_reference as t
import make_mesh

OMEGAS = (1.0, 1.4)
N_MODES = 16000             # converged to 1e-11 against 32000 at 0.1 m from the waterline (phi, grad, phi_zz)
RINGS = (10.1, 10.2, 10.5, 11.0, 12.0, 13.5)
ANGLES = (0.0, 45.0, 90.0, 135.0, 180.0)
CIRCLE = ((15.5, 20.0), (0.0, 60.0, 120.0, 180.0))
C = 1j


def read_pnl(path):
    lines = open(path).read().splitlines()
    i = next(j for j, l in enumerate(lines) if "Number of Panels" in l)
    npan, nnode = (int(v) for v in lines[i + 1].split()[:2])
    i = next(j for j, l in enumerate(lines) if "Start Definition of Node Coordinates" in l)
    nodes = np.array([[float(v) for v in lines[i + 1 + n].split()[1:4]] for n in range(nnode)])
    i = next(j for j, l in enumerate(lines) if "Start Definition of Node Relations" in l)
    cents = []
    for n in range(npan):
        p = [int(v) for v in lines[i + 1 + n].split()]
        nv = p[1]
        cents.append(nodes[np.array(p[2:2 + nv]) - 1].mean(axis=0))
    return np.array(cents)


def fields(sol, x, y):
    """(phi, gx, gy, gz, phizz) of sol's radiating part on z = 0 at (x, y) arrays (_z0_modal, _theta_fields)."""
    r, th = np.hypot(x, y), np.arctan2(y, x)
    out = [np.zeros(x.size, complex) for _ in range(5)]
    for i in range(x.size):
        modal = t._z0_modal(sol, np.array([r[i]]))
        phi, ur, ut, uz, zz = (v[0, 0] for v in t._theta_fields(modal, np.array([r[i]]), sol.omega ** 2 / make_mesh.G,
                                                                  np.array([th[i]]), sorted(sol.orders)))
        c, s = math.cos(th[i]), math.sin(th[i])
        for j, v in enumerate((phi, c * ur - s * ut, s * ur + c * ut, uz, zz)):
            out[j][i] = v
    return out


def incident(w, k, x, y):
    phi = -make_mesh.G / w * np.exp(1j * k * x)
    return phi, 1j * k * phi, 0 * phi, w ** 2 / make_mesh.G * phi, k ** 2 * phi


# The grid of fs_fields_exact_integral: the inner free-surface integral of the sum pair (1.4, 1.0) rad/s, Omega+ = 2.4,
# over a < r < a + 5 m (where the waterline matters) with the exact psi_k, from HAMS's fields against the exact ones.
# Radial: Gauss-Legendre (GRID_GAUSS nodes) on panels graded towards the waterline; angular: theta = pi j / 16,
# j = 0..16, trapezoidal over [0, pi] doubled (the integrand is even in theta for heading 0, and the products have
# orders below 32, so the rule is exact in theta). Weight r dr dtheta.
GRID_EDGES = (0.0, 0.25, 0.5, 1.0, 2.0, 3.5, 5.0)
GRID_GAUSS = 4
GRID_NTH = 16
GRID_PAIR = (1.4, 1.0)


def grid(out, kernel="sum"):
    cyl = t.Cylinder()
    res = dict(t.FS_CHECK_RESOLUTION, n_modes=N_MODES)
    g = make_mesh.G
    w1, w2 = GRID_PAIR
    plus = kernel == "sum"
    big = w1 + w2 if plus else w1 - w2
    sols = {w: t._first_order_case(cyl, w, "fixed", res["n_modes"], res["p_max"], res["orders"])[0] for w in GRID_PAIR}
    rads = t._psi_solutions(cyl, big, res["n_modes"], t._spec_key(res["p_max"]))
    xg, wg = np.polynomial.legendre.leggauss(GRID_GAUSS)
    rr, wr = [], []
    for lo, hi in zip(GRID_EDGES[:-1], GRID_EDGES[1:]):
        rr += list(cyl.radius + 0.5 * (hi - lo) * xg + 0.5 * (hi + lo))
        wr += list(0.5 * (hi - lo) * wg)
    th = math.pi * np.arange(GRID_NTH + 1) / GRID_NTH
    wt = np.full(th.size, 2.0 * math.pi / GRID_NTH)          # doubled half-trapezoid: 2 x (pi/16), ends halved
    wt[0] = wt[-1] = math.pi / GRID_NTH
    x = np.array([r * math.cos(a) for r in rr for a in th])
    y = np.array([r * math.sin(a) for r in rr for a in th])
    wts = np.array([wri * ri * wa for ri, wri in zip(rr, wr) for wa in wt])
    tot = {}
    for w in GRID_PAIR:
        sc, inc = fields(sols[w], x, y), incident(w, sols[w].k, x, y)
        tot[w] = ([sc[j] + inc[j] for j in range(5)], inc)
    psi = []
    for k in (1, 3, 5):
        (m,) = rads[k].orders
        psi.append(np.array([t._z0_modal(rads[k], np.array([np.hypot(x[i], y[i])]))[0][0, 0] for i in range(x.size)])
                   * np.cos(m * np.arctan2(y, x)))
    # sum: Q+ with state j as it is, HAMS's T = C**2 T_ref; difference: Q- with the conjugate of state j, C conj(C) = 1
    cj = (lambda v: tuple(v)) if plus else (lambda v: tuple(np.conj(x) for x in v))
    q = t._fs_forcing(w1, tuple(tot[w1][0]), w2, cj(tot[w2][0]), plus=plus) - \
        t._fs_forcing(w1, tuple(tot[w1][1]), w2, cj(tot[w2][1]), plus=plus)
    factor = C ** 2 if plus else C * np.conj(C)
    expect = [factor * t.FS_SIGN * 1j * big / g ** 2 * np.sum(wts * p * q) for p in psi]
    lines = ["# fs_fields_exact_grid: written by make_fs_fields_exact_fixture.py --grid%s (GRID_* there), kernel %s"
             % ("" if plus else "-diff", kernel),
             "GRIDPAIR %.10f %.10f %.16e %.16e %.16e %.16e" % (w1, w2, sols[w1].k, sols[w2].k, g, cyl.depth),
             "GRID %d" % x.size]
    for i in range(x.size):
        lines.append("%.16e %.16e %.16e " % (x[i], y[i], wts[i]) + " ".join("%.16e %.16e" % (p[i].real, p[i].imag) for p in psi))
    lines.append("EXPECT_GRID " + " ".join("%.16e %.16e" % (v.real, v.imag) for v in expect))
    lines.append("END")
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote %s (%d points); expected T_FS inner (a, a + 5 m) %s" % (out, x.size, np.round(expect, 4)))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data = os.path.join(here, "..", "..", "test", "data", "fs_sum_exact")
    if "--grid" in sys.argv:
        return grid(os.path.join(data, "fs_fields_exact_grid.txt"))
    if "--grid-diff" in sys.argv:
        return grid(os.path.join(data, "fs_fields_exact_grid_diff.txt"), kernel="diff")
    hull = sys.argv[1] if len(sys.argv) > 1 else os.path.join(data, "HullMesh_sharp_m1.pnl")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(data, "fs_fields_exact.txt")
    cents = read_pnl(hull)
    cyl = t.Cylinder()
    res = dict(t.FS_CHECK_RESOLUTION, n_modes=N_MODES)
    lines = ["# fs_fields_exact: written by pre_postprocessing_python/SecondOrder/make_fs_fields_exact_fixture.py "
             "(header there); hull %s" % os.path.basename(hull)]
    for w in OMEGAS:
        sol = t._first_order_case(cyl, w, "fixed", res["n_modes"], res["p_max"], res["orders"])[0]
        lines.append("FREQ %.10f %.16e %.16e %.16e" % (w, sol.k, make_mesh.G, cyl.depth))
        lines.append("PANELS %d" % len(cents))
        for xc in cents:
            r, th, z = math.hypot(xc[0], xc[1]), math.atan2(xc[1], xc[0]), xc[2]
            # the centroid of a wall panel is a chord midpoint inside the cylinder (r = a cos(pi/n)), where the exterior
            # series diverges: take the exact value on the hull at its azimuth and depth (as maccamy_fuchs_state); a
            # bottom panel (z = -d) keeps its r < a, inside the gap region
            if abs(z + cyl.draft) > 1e-6:
                r = cyl.radius
            tot = complex(np.asarray(sol.potential(np.array([r]), np.array([th]), np.array([z]))).ravel()[0])
            kk = sol.k
            inc = -make_mesh.G / w * math.cosh(kk * (z + cyl.depth)) / math.cosh(kk * cyl.depth) * np.exp(1j * kk * r * math.cos(th))
            v = C * (tot - inc)
            lines.append("%.16e %.16e" % (v.real, v.imag))
        x = np.array([r * math.cos(math.radians(a)) for r in RINGS for a in ANGLES])
        y = np.array([r * math.sin(math.radians(a)) for r in RINGS for a in ANGLES])
        sc, inc = fields(sol, x, y), incident(w, sol.k, x, y)
        tot = [sc[j] + inc[j] for j in range(5)]
        lines.append("POINTS %d" % x.size)
        for i in range(x.size):
            lines.append("%.16e %.16e " % (x[i], y[i]) + " ".join("%.16e %.16e" % ((C * tot[j][i]).real, (C * tot[j][i]).imag)
                                                                  for j in range(5)))
        x = np.array([r * math.cos(math.radians(a)) for r in CIRCLE[0] for a in CIRCLE[1]])
        y = np.array([r * math.sin(math.radians(a)) for r in CIRCLE[0] for a in CIRCLE[1]])
        sc = fields(sol, x, y)
        lines.append("CIRCLE %d" % x.size)
        for i in range(x.size):
            lines.append("%.16e %.16e " % (x[i], y[i]) + " ".join("%.16e %.16e" % ((C * sc[j][i]).real, (C * sc[j][i]).imag)
                                                                  for j in range(5)))
        lines.append("END")
        print("omega %.2f done" % w, flush=True)
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote %s (%d lines)" % (out, len(lines)))


if __name__ == "__main__":
    main()
