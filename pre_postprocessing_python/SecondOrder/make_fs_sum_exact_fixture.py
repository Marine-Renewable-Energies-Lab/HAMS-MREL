"""The fixture of test_qtf_free_surface's fs_sum_exact_* tests: the sum-frequency free-surface term T_FS+ of HAMS's
assembly routines (fs_inner_pair, fs_outer_rfar + fs_outer_pair) with the exact first-order fields and psi_k of the
sharp truncated cylinder (truncated_cylinder_reference.py: radius 10, draft 10, depth 100 m, fixed body, heading 0)
prescribed, against the reference's own T_FS+ of the same fields. No BEM solve: it tests the assembly alone.

    python make_fs_sum_exact_fixture.py [out]      (default ../../test/data/fs_sum_exact/fs_sum_exact.txt)

Conventions. HAMS's incident wave is phi_I = -i (g/omega) e^{ikx} (FirstOrderState.incident_phi), the reference's
-(g/omega) e^{ikx}: HAMS's state is c = i times the reference's. The scattered fields are written times c, so that with
HAMS's own incident wave they are HAMS's states; Q is bilinear, so HAMS's T = c**2 T_ref = -T_ref, and the expected
values below are written already times c**2, over rho g (HAMS's t / (rho g)), components psi_1, psi_3, psi_5.
psi_k (unit velocity) is the same in both codes (QTF_PsiCheck against the exact A, B: 2 %, findings section 6).

Inner: rings r = RINGS, NTH points in theta each, weight 1: expected (i Omega/g**2) sum psi_k Q_D(x) with the
reference's _fs_forcing (plus), Q_D = Q[tot, tot] - Q[I, I]. Outer: the partial waves of the scattered states and of
psi_k on R = OUTER_R (pw_set: a(m, 0) of H_m(k r) e^{i m theta}, a(m, n) of K_m(k_n r) e^{i m theta}, from the cos(m theta)
series: a(+-m) = c_m / 2 with (-1)**m on -m for H, since H_-m = (-1)**m H_m and K_-m = K_m), NEV evanescent modes;
expected: the reference's fs_integral over its outer nodes and tail contours on a grid of radius OUTER_R.

File: one block per pair. Line 'PAIR w1 w2 Omega k1 k2 K g h' (w1 >= w2), 'INNER n' then n lines
'x y  Re Im of (phi, gx, gy, gz, phizz) of scattered 1, of scattered 2, then psi_1 psi_3 psi_5', 'EXPECT_INNER' 3 Re Im,
'OUTER R mmax nev', 'KEV set' lines (the nev k_n of each of the 5 sets: state 1, state 2, psi_1, psi_3, psi_5),
'A set m n Re Im' for every nonzero coefficient, 'EXPECT_OUTER' 3 Re Im, 'END'."""
import math, os, sys
import numpy as np
from scipy.special import hankel1, kve
import truncated_cylinder_reference as t
import make_mesh

PAIRS = ((1.2, 1.2), (1.4, 1.0))      # (w1, w2), w1 >= w2: the worst sum pairs of the fixed m2 run (Omega+ = 2.4)
N_MODES = 4000                         # FS_SUM_N, the resolution of fs_sum_bvp
RINGS = (10.3, 10.6, 11.0, 12.0, 13.5)
NTH = 72
OUTER_R = 15.0                         # HAMS's default split radius (qtf_fs_radius)
NEV = 80                               # e^{-k_80 (R - a)} = e^{-12.6}
C = 1j                                 # HAMS's state / the reference's


def _xy_fields(modal, r, nu, theta, orders):
    phi, ur, ut, uz, zz = t._theta_fields(modal, np.array([r]), nu, theta, orders)
    phi, ur, ut, uz, zz = (v[0] for v in (phi, ur, ut, uz, zz))
    return phi, np.cos(theta) * ur - np.sin(theta) * ut, np.sin(theta) * ur + np.cos(theta) * ut, uz, zz


def _incident(state, x, y):
    """The reference's incident fields (phi, gx, gy, gz, phizz) on z = 0 at (x, y), heading 0."""
    k, amp = state.k, -make_mesh.G / state.omega
    phi = amp * np.exp(1j * k * x)
    return phi, 1j * k * phi, 0 * phi, state.omega ** 2 / make_mesh.G * phi, k ** 2 * phi


def _pw_coefficients(sol, nev):
    """{(m, n): a} of sol's radiating part on z = 0 in pw_set form (header)."""
    a_rad = sol.cyl.radius
    out = {}
    for m, o in sol.orders.items():
        c0 = o.scattered / hankel1(m, sol.k * a_rad)
        # K_m(x) = kve(m, x) e^{-x}, so b K_m(k_n r)/K_m(k_n a) = [b / K_m(k_n a)] K_m(k_n r)
        cn = [o.evanescent[n] / (kve(m, sol.kn[n] * a_rad) * math.exp(-sol.kn[n] * a_rad)) for n in range(min(nev, len(sol.kn)))]
        for n, v in enumerate([c0] + cn):
            if v == 0:
                continue
            if m == 0:
                out[(0, n)] = v
            else:
                out[(m, n)] = v / 2.0
                out[(-m, n)] = v / 2.0 * ((-1) ** m if n == 0 else 1.0)
    return out


def block(w1, w2):
    cyl = t.Cylinder()
    res = dict(t.FS_CHECK_RESOLUTION, n_modes=N_MODES)
    g, h = make_mesh.G, cyl.depth
    big = round(w1 + w2, 9)
    tot = {w: t._first_order_case(cyl, w, "fixed", res["n_modes"], res["p_max"], res["orders"])[0] for w in (w1, w2)}
    rads = t._psi_solutions(cyl, big, res["n_modes"], t._spec_key(res["p_max"]))
    grid = t.fs_grid(cyl, OUTER_R, (w2, w1), (big,), res["n_modes"])
    states = {w: t.fs_state(grid, w, tot[w]) for w in (w1, w2)}
    lines = ["PAIR %.10f %.10f %.10f %.16e %.16e %.16e %.16e %.16e" % (w1, w2, big, tot[w1].k, tot[w2].k, rads[1].k, g, h)]
    # inner rings
    theta = 2 * math.pi * np.arange(NTH) / NTH
    rows, expect = [], np.zeros(3, complex)
    for r in RINGS:
        x, y = r * np.cos(theta), r * np.sin(theta)
        sc = {w: _xy_fields(t._z0_modal(tot[w], np.array([r])), r, w ** 2 / g, theta, sorted(tot[w].orders)) for w in (w1, w2)}
        inc = {w: _incident(states[w], x, y) for w in (w1, w2)}
        totf = {w: tuple(sc[w][i] + inc[w][i] for i in range(5)) for w in (w1, w2)}
        psi = []
        for k in (1, 3, 5):
            (m,) = rads[k].orders
            psi.append(t._z0_modal(rads[k], np.array([r]))[0][0, 0] * np.cos(m * theta))
        # _fs_forcing takes (phi, u_r, u_theta, u_z, phi_zz); Q is a sum over the three gradient components, which is
        # invariant under the rotation (r, theta) -> (x, y), so the x, y fields go in the same slots
        q = (t._fs_forcing(w1, totf[w1], w2, totf[w2], plus=True) - t._fs_forcing(w1, inc[w1], w2, inc[w2], plus=True))
        for i, k in enumerate((1, 3, 5)):
            expect[i] += 1j * big / g ** 2 * np.sum(psi[i] * q) * t.FS_SIGN
        for j in range(NTH):
            vals = [C * sc[w1][i][j] for i in range(5)] + [C * sc[w2][i][j] for i in range(5)] + [p[j] for p in psi]
            rows.append("%.16e %.16e " % (x[j], y[j]) + " ".join("%.16e %.16e" % (v.real, v.imag) for v in vals))
    lines.append("INNER %d" % len(rows))
    lines += rows
    lines.append("EXPECT_INNER " + " ".join("%.16e %.16e" % (v.real, v.imag) for v in C ** 2 * expect))
    # outer
    sets = [tot[w1], tot[w2], rads[1], rads[3], rads[5]]
    coefs = [_pw_coefficients(s, NEV) for s in sets]
    # self-check: the pw_set form reproduces the exact radiating part on z = 0 at R and at 1.5 R
    from scipy.special import kv
    for s, c in zip(sets, coefs):
        for r in (OUTER_R, 1.5 * OUTER_R):
            th = np.array([0.0, 0.7, 2.0])
            ref = t._theta_fields(t._z0_modal(s, np.array([r])), np.array([r]), 0.0, th, sorted(s.orders))[0][0]
            mine = sum(v * (hankel1(m, s.k * r) if n == 0 else kv(m, s.kn[n - 1] * r)) * np.exp(1j * m * th)
                       for (m, n), v in c.items())
            err = np.max(np.abs(mine - ref)) / np.max(np.abs(ref))
            if err > 1e-6:
                raise SystemExit("pw coefficients miss the exact field by %.1e at r = %g" % (err, r))
    scale = [C, C, 1.0, 1.0, 1.0]
    mmax = max(abs(m) for c in coefs for (m, n) in c)
    lines.append("OUTER %.10f %d %d" % (OUTER_R, mmax, NEV))
    for i, s in enumerate(sets):
        lines.append("KEV %d %.16e " % (i + 1, s.k) + " ".join("%.16e" % v for v in s.kn[:NEV]))
    for i, c in enumerate(coefs):
        for (m, n), v in sorted(c.items()):
            v = scale[i] * v
            lines.append("A %d %d %d %.16e %.16e" % (i + 1, m, n, v.real, v.imag))
    pair = t.fs_pair(grid, states[w1], states[w2], kernel="sum")
    expect = np.zeros(3, complex)
    for i, k in enumerate((1, 3, 5)):
        test = t.fs_psi_test(grid, rads[k])
        m = test.m
        total = np.sum(grid.outer[1] * grid.outer[0] * test.outer * pair.outer[m])
        for sigma, kappa, r, w, values in pair.tail[m]:
            total += np.sum(1j * sigma * w * test.tau * t.hankel1e(m, pair.nu * r) * values * r * np.exp(1j * kappa * r))
        expect[i] = t.FS_SIGN * 1j * big / g ** 2 * total
    lines.append("EXPECT_OUTER " + " ".join("%.16e %.16e" % (v.real, v.imag) for v in C ** 2 * expect))
    lines.append("END")
    return lines


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "test", "data", "fs_sum_exact", "fs_sum_exact.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    lines = ["# fs_sum_exact: written by pre_postprocessing_python/SecondOrder/make_fs_sum_exact_fixture.py (header there)"]
    for w1, w2 in PAIRS:
        lines += block(w1, w2)
        print("pair (%.2f, %.2f) done" % (w1, w2), flush=True)
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote %s (%d lines)" % (out, len(lines)))


if __name__ == "__main__":
    main()
