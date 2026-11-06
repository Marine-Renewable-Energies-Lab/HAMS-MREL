"""Semi-analytical reference for the damping lid: a bottom-mounted vertical cylinder with a phi = 0 ring around it.

Linear potential flow in water of depth h. A cylinder of radius a stands on the sea bed and pierces the free surface. The total
potential is zero on the annulus a < r < b of z = 0 (the lid). Outside it, r > b, the free surface condition dphi/dz = nu phi
holds, nu = omega**2/g. A plane wave of unit amplitude comes along +x. The problem is solved by eigenfunction matching, one
azimuthal order m at a time, with no boundary element and no Green function:

- lid region a < r < b: vertical modes L_j = cos(lam_j (z + h)), lam_j = (j - 1/2) pi/h, j = 1..N, which vanish on z = 0 and
  have dL/dz = 0 on z = -h. Radial functions R_j = I_m(lam r) K_m'(lam a) - K_m(lam r) I_m'(lam a), so dphi/dr = 0 at r = a.
- outer region r > b: the propagating mode Z_0 = cosh k(z + h)/cosh kh with J_m(kr) (incident) and H_m^(1)(kr) (scattered),
  and the evanescent modes Z_n = cos k_n(z + h)/cos k_n h with K_m(k_n r), n = 1..N-1, k_n tan(k_n h) = -nu.
- matching at r = b: phi and dphi/dr continuous on -h < z < 0. Both equations are projected on the N lid modes L_j ('lid',
  the reference values). The script also solves with both projected on the outer modes ('outer'), and with the two mixed
  choices, and prints the spread between them. 'lid' is used because it is the only one whose pointwise mismatch of both phi
  and dphi/dr at r = b falls with N (check 3): with 'outer' the lid coefficients c_j gamma_j do not decay with j, and the
  pointwise dphi/dr of the lid side stays wrong (even in sign) although every projected quantity converges. All four give the
  same extrapolated force and elevations to 2e-5.

All inner products are closed form: int_{-h}^0 Z_n L_j dz = lam_j (-1)**(j+1)/(lam_j**2 + kap_n**2), kap_0 = k, kap_n = i k_n.

Time factor and output normalisation, from src/PatclVelct.f90 (VINP), src/PotentWavForce.f90 (EFORCE) and
src/PressureElevation.f90 (CalElevation, WamitNondimens): HAMS solves with exp(-i omega t) and an incident potential whose
elevation at the origin is -i A. ExcForce.3 writes (-Im X, -Re X)/(rho g A) and PressureElevation.6p writes (-Im eta, -Re eta)/A,
which is conj(X/(rho g eta_I(0))) and conj(eta/eta_I(0)): the WAMIT convention, exp(+i omega t) and phases relative to the
incident elevation at the origin. Here psi is the potential scaled so that the incident wave is exp(ikx) Z_0(z), so that
eta/eta_I(0) = psi(z = 0), and the force on the body over rho g eta_I(0) is -int psi n_x dS (n out of the body). The script prints
conj of both.

Self-checks; the script exits with 1 if one fails. What each covers, and what it does not:

1. b = a: X1 and the elevations equal MacCamy-Fuchs to 1e-8, for every projection. With b = a, gamma = 0 and the system gives
   a_0 = -dinc_0/beta_0 exactly, and the force takes the outer-mode branch. So it covers the outer modes, beta, the incident
   term, the force and elevation sums and the output convention, but not M, the lid radial functions, gamma or the lid force.
2. Thin lid b = a (1 + eps), eps = 1e-1 .. 1e-6, at N = 80, 320, 1280: X1 through the lid force branch approaches MacCamy-Fuchs,
   monotonically as eps falls, and the remaining difference at eps = 1e-6 falls with N. It runs M, the lid radial functions,
   gamma and the lid force sum together, but only in the limit where the lid does nothing, so it does not test how strongly
   the lid acts.
3. Pointwise matching at r = b, at the depths ZS, orders 0, 1, 2, 5: the mismatch of phi and of dphi/dr between the outer and
   the lid expansions, each evaluated from its own modes, falls with N (N = 20, 80, 320). It checks M, gamma, beta and the
   incident terms against the conditions the solution must meet, at a real lid. It does not check the lid force sum or the
   wall condition dphi/dr = 0 at r = a (which the radial functions meet by construction).
4. Energy: phi = 0 on the lid absorbs nothing, so every order keeps abs(1 + 2 a_0/(eps_m i**m H_m(kb))) = 1. This holds for any
   real M and gamma, so it does not detect a wrong sign or value in them; it checks only that the system is assembled with real
   lid quantities and the right outer radiation condition.
5. Orders: 0..29 against 0..59 at the two points, without and with the lid, below 1e-10.

The spread between the four projections and between extrapolations from different N (printed) bounds the truncation error of
the reference values; it is not a check of the formulation.

Run: python -I annular_lid_reference.py. It takes a few seconds.
"""
import math
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import h1vp, hankel1, ive, jv, jvp, kve

A_RADIUS = 2.0   # cylinder radius a
B_RADIUS = 4.0   # outer radius of the lid b
DEPTH = 10.0     # h
KH = 1.5         # k h; HAMS is given the wave number k = 0.15 directly (Input_frequency_type 2)
G = 9.80665      # g of src/WavDynMods.f90, used only to print omega
ORDERS = 30      # azimuthal orders m = 0..ORDERS-1 for the elevation (the force needs only m = 1)
POINTS = [(0.0, 8.0), (8.0, 0.0)]
ZS = [-1.0, -2.5, -5.0, -7.5, -9.0]   # depths of the pointwise matching check
PROJECTIONS = ("lid", "outer", "phi_lid", "phi_outer")
DEFAULT = PROJECTIONS[0]   # the projection of the reference values; see the module docstring


def evanescent_roots(nu, h, count):
    """The roots k_n, n = 1..count, of k tan(k h) = -nu, k_n h in ((n - 1/2) pi, n pi). Solved as y = n pi - atan(nu h/y),
    y = k_n h, which has no pole in the bracket."""
    out = []
    for n in range(1, count + 1):
        y = brentq(lambda y: y - n * math.pi + math.atan(nu * h / y), (n - 0.5) * math.pi, n * math.pi, xtol=1e-15, rtol=1e-15,
                   maxiter=200)
        out.append(y / h)
    return np.array(out)


class Problem:
    """One geometry and frequency. b == a gives the bare cylinder."""

    def __init__(self, a, b, h, k, n):
        self.a, self.b, self.h, self.k, self.n = a, b, h, k, n
        self.nu = k * math.tanh(k * h)
        self.kn = evanescent_roots(self.nu, h, n - 1)
        self.lam = (np.arange(1, n + 1) - 0.5) * math.pi / h
        kap2 = np.concatenate(([k * k], -self.kn ** 2))                       # kap_n**2
        sign = (-1.0) ** np.arange(n)                                          # (-1)**(j+1), j = 1..n
        self.M = self.lam[None, :] * sign[None, :] / (self.lam[None, :] ** 2 + kap2[:, None])   # M[n, j]
        nrm = np.empty(n)
        nrm[0] = (h / 2 + math.sinh(2 * k * h) / (4 * k)) / math.cosh(k * h) ** 2
        nrm[1:] = (h / 2 + np.sin(2 * self.kn * h) / (4 * self.kn)) / np.cos(self.kn * h) ** 2
        self.N = nrm                                                           # int Z_n**2
        self.P = h / 2                                                         # int L_j**2
        self.zint_lid = sign / self.lam                                        # int L_j dz
        self.zint_out = np.concatenate(([math.tanh(k * h) / k], np.tan(self.kn * h) / self.kn))   # int Z_n dz

    # radial functions, each normalised to 1 at r = b
    def lid_radial(self, m, r):
        """R_j(r)/R_j(b) and R_j'(b)/R_j(b) for every j, from exponentially scaled Bessel functions."""
        lam, a, b = self.lam, self.a, self.b
        ipe_a = 0.5 * (ive(m - 1, lam * a) + ive(m + 1, lam * a))              # I_m'(lam a) exp(-lam a)
        kpe_a = -0.5 * (kve(m - 1, lam * a) + kve(m + 1, lam * a))             # K_m'(lam a) exp(lam a)

        def s(x):                                                              # R(x) exp(-lam (x - a))
            return ive(m, lam * x) * kpe_a - np.exp(-2 * lam * (x - a)) * kve(m, lam * x) * ipe_a

        def sp(x):                                                             # R'(x) exp(-lam (x - a))/lam
            ip = 0.5 * (ive(m - 1, lam * x) + ive(m + 1, lam * x))
            kp = -0.5 * (kve(m - 1, lam * x) + kve(m + 1, lam * x))
            return ip * kpe_a - np.exp(-2 * lam * (x - a)) * kp * ipe_a

        ratio = np.exp(lam * (r - b)) * s(r) / s(b)
        gamma = lam * sp(b) / s(b)
        return ratio, gamma

    def outer_beta(self, m):
        """d/dr of the outer radial functions at r = b over their value there, n = 0..N-1."""
        k, b = self.k, self.b
        beta = np.empty(self.n, dtype=complex)
        beta[0] = k * h1vp(m, k * b) / hankel1(m, k * b)
        kn = self.kn
        beta[1:] = -0.5 * kn * (kve(m - 1, kn * b) + kve(m + 1, kn * b)) / kve(m, kn * b)
        return beta

    def solve(self, m, projection=DEFAULT):
        """Coefficients of order m: a[n] of the outer modes (radial functions 1 at r = b, scattered part only) and c[j] of the
        lid modes (radial functions 1 at r = b)."""
        n, M, N, P = self.n, self.M, self.N, self.P
        s = (1.0 if m == 0 else 2.0) * 1j ** m
        inc = np.zeros(n, dtype=complex)
        dinc = np.zeros(n, dtype=complex)
        inc[0] = s * jv(m, self.k * self.b)
        dinc[0] = s * self.k * jvp(m, self.k * self.b)
        beta = self.outer_beta(m)
        if self.b == self.a:
            gamma = np.zeros(n)
        else:
            gamma = self.lid_radial(m, self.b)[1]
        if projection == "outer":       # phi and dphi/dr on Z_p
            lhs = M * (beta[:, None] - gamma[None, :])
            rhs = N * (beta * inc - dinc)
            c = np.linalg.solve(lhs, rhs)
            a = M @ c / N - inc
        elif projection == "lid":       # phi and dphi/dr on L_q
            lhs = (M * (beta[:, None] - gamma[None, :])).T
            rhs = -(dinc[0] - gamma * inc[0]) * M[0, :]
            a = np.linalg.solve(lhs, rhs)
            c = M.T @ (inc + a) / P
        elif projection == "phi_lid":   # phi on L_q, dphi/dr on Z_p
            W = (M * gamma[None, :]) @ M.T / P                                 # W[p, n] = sum_q gamma_q M_pq M_nq / P
            lhs = np.diag(beta * N) - W
            rhs = W[:, 0] * inc[0] - dinc * N
            a = np.linalg.solve(lhs, rhs)
            c = M.T @ (inc + a) / P
        elif projection == "phi_outer":  # phi on Z_p, dphi/dr on L_q
            lhs = (M.T * beta[None, :]) @ (M / N[:, None]) - np.diag(gamma * P)
            rhs = -(dinc[0] - beta[0] * inc[0]) * M[0, :]
            c = np.linalg.solve(lhs, rhs)
            a = M @ c / N - inc
        else:
            raise ValueError(projection)
        return a, c, inc

    def surge(self, projection=DEFAULT):
        """Surge force over rho g eta_I(0), HAMS output convention (conjugated)."""
        a, c, inc = self.solve(1, projection)
        if self.b == self.a:
            wall = (inc[0] + a[0]) * self.zint_out[0] + np.sum(a[1:] * self.zint_out[1:])
        else:
            ratio = self.lid_radial(1, self.a)[0]
            wall = np.sum(c * ratio * self.zint_lid)
        return np.conj(-math.pi * self.a * wall)

    def elevation(self, x, y, projection=DEFAULT, orders=ORDERS):
        """Total elevation over eta_I(0) at (x, y, 0), r > b, HAMS output convention (conjugated)."""
        r, th = math.hypot(x, y), math.atan2(y, x)
        k, b, kn = self.k, self.b, self.kn
        total = 0.0
        for m in range(orders):
            a, c, inc = self.solve(m, projection)
            s = (1.0 if m == 0 else 2.0) * 1j ** m
            val = s * jv(m, k * r) + a[0] * hankel1(m, k * r) / hankel1(m, k * b)
            val += np.sum(a[1:] * kve(m, kn * r) / kve(m, kn * b) * np.exp(-kn * (r - b)))
            total += val * math.cos(m * th)
        return np.conj(total)

    def matching_residual(self, m, zs, projection=DEFAULT):
        """Pointwise mismatch at r = b of phi and of dphi/dr between the outer and the lid expansions of order m, at the depths
        zs, each over the largest outer value there. Uses the lid modes, gamma and the outer modes, independently of the
        projection."""
        a, c, inc = self.solve(m, projection)
        beta = self.outer_beta(m)
        gamma = self.lid_radial(m, self.b)[1]
        z = np.asarray(zs)[:, None]
        zo = np.concatenate([np.cosh(self.k * (z + self.h)) / math.cosh(self.k * self.h),
                             np.cos(self.kn[None, :] * (z + self.h)) / np.cos(self.kn * self.h)[None, :]], axis=1)
        lo = np.cos(self.lam[None, :] * (z + self.h))
        phi_out = zo @ (inc + a)
        phi_lid = lo @ c
        dr_out = zo[:, 0] * self.k * jvp(m, self.k * self.b) * (1.0 if m == 0 else 2.0) * 1j ** m + zo @ (a * beta)
        dr_lid = lo @ (c * gamma)
        return (np.max(np.abs(phi_out - phi_lid)) / np.max(np.abs(phi_out)),
                np.max(np.abs(dr_out - dr_lid)) / np.max(np.abs(dr_out)))

    def energy_defect(self, m, projection=DEFAULT):
        """abs(1 + 2 a_0/(s H_m(kb))) - 1, zero when no energy is lost or made."""
        a, _, _ = self.solve(m, projection)
        s = (1.0 if m == 0 else 2.0) * 1j ** m
        return abs(1 + 2 * a[0] / (s * hankel1(m, self.k * self.b))) - 1


def maccamy_fuchs_surge(a, h, k):
    return np.conj(4 * math.tanh(k * h) / (k * k * h1vp(1, k * a)))


def maccamy_fuchs_elevation(a, k, x, y, orders=ORDERS):
    r, th = math.hypot(x, y), math.atan2(y, x)
    total = 0.0
    for m in range(orders):
        s = (1.0 if m == 0 else 2.0) * 1j ** m
        total += s * (jv(m, k * r) - jvp(m, k * a) / h1vp(m, k * a) * hankel1(m, k * r)) * math.cos(m * th)
    return np.conj(total)


def values(prob, projection=DEFAULT, orders=ORDERS):
    return [prob.surge(projection)] + [prob.elevation(x, y, projection, orders) for x, y in POINTS]


def fmt(z):
    return "%.6f (%+.6f %+.6fi, phase %+.3f deg)" % (abs(z), z.real, z.imag, math.degrees(math.atan2(z.imag, z.real)))


def main():
    a, b, h = A_RADIUS, B_RADIUS, DEPTH
    k = KH / h
    omega = math.sqrt(G * k * math.tanh(k * h))
    print("a = %g, b = %g, h = %g, k = %.6f (kh = %g), omega = %.10f rad/s, nu = %.10f" % (a, b, h, k, KH, omega, k * math.tanh(k * h)))
    labels = ["surge X1", "eta(0, 8, 0)", "eta(8, 0, 0)"]
    failed = []

    def require(cond, what):
        if not cond:
            failed.append(what)

    print("\nCheck 1, b = a against MacCamy-Fuchs (outer modes only; the lid modes, gamma and the lid force are not used):")
    ref = [maccamy_fuchs_surge(a, h, k)] + [maccamy_fuchs_elevation(a, k, x, y) for x, y in POINTS]
    for proj in PROJECTIONS:
        got = values(Problem(a, a, h, k, 20), proj)
        err = max(abs(g - r) / abs(r) for g, r in zip(got, ref))
        require(err < 1e-8, "b = a, " + proj)
        print("  %-9s worst relative difference %.2e" % (proj, err))
    for lab, r in zip(labels, ref):
        print("  MacCamy-Fuchs %-13s %s  re %.8f  im %.8f" % (lab, fmt(r), r.real, r.imag))
    tail = max(abs(maccamy_fuchs_elevation(a, k, x, y, ORDERS) - maccamy_fuchs_elevation(a, k, x, y, 2 * ORDERS)) for x, y in POINTS)
    print("  orders 0..%d against 0..%d, worst of the two points: %.1e" % (ORDERS - 1, 2 * ORDERS - 1, tail))

    print("\nCheck 2, a thin lid b = a (1 + eps) through the lid path (lid modes, R_j, gamma, lid force sum), projection '%s'."
          "\n  Relative difference of X1 to MacCamy-Fuchs; it must fall as eps -> 0. The approach is slow, as a phi = 0 strip acts"
          "\n  like a conductor, and it stops at the truncation error of the lid modes at the wall, which do not reach z = 0:" % DEFAULT)
    epss = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
    floor = []
    for n in (80, 320, 1280):
        seq = [abs(Problem(a, a * (1 + e), h, k, n).surge() - ref[0]) / abs(ref[0]) for e in epss]
        print("  N = %4d  " % n + "  ".join("eps %.0e: %.1e" % (e, d) for e, d in zip(epss, seq)))
        require(all(x > y for x, y in zip(seq, seq[1:])), "thin lid falls with eps, N = %d" % n)
        floor.append(seq[-1])
    require(floor[0] > floor[1] > floor[2], "thin lid floor falls with N")

    print("\nCheck 3, pointwise matching at r = b, z = %s, b = %g: largest mismatch of phi and of dphi/dr over m = 0, 1, 2, 5,"
          "\n  each relative to the largest outer value. Covers M, gamma, beta and the incident terms pointwise:" % (ZS, b))
    for proj in PROJECTIONS:
        rows = []
        for n in (20, 80, 320):
            res = [Problem(a, b, h, k, n).matching_residual(m, ZS, proj) for m in (0, 1, 2, 5)]
            rows.append((max(r[0] for r in res), max(r[1] for r in res)))
        print("  %-9s " % proj + "  ".join("N = %d: phi %.1e, dphi/dr %.1e" % (n, r[0], r[1]) for n, r in zip((20, 80, 320), rows)))
        if proj == DEFAULT:
            require(rows[0][0] > rows[1][0] > rows[2][0] and rows[0][1] > rows[1][1] > rows[2][1], "matching, " + proj)

    print("\nLid a < r < b, projection '%s', convergence in N:" % DEFAULT)
    table = {}
    for n in (20, 40, 80, 160, 320):
        prob = Problem(a, b, h, k, n)
        table[n] = values(prob)
        defect = max(abs(prob.energy_defect(m)) for m in range(ORDERS))
        require(defect < 1e-8, "energy, N = %d" % n)
        print("  N = %3d  energy defect %.1e" % (n, defect))
        for lab, v in zip(labels, table[n]):
            print("      %-13s %s" % (lab, fmt(v)))
    tail = max(abs(x - y) for x, y in zip(values(Problem(a, b, h, k, 80))[1:], values(Problem(a, b, h, k, 80), orders=2 * ORDERS)[1:]))
    print("  orders 0..%d against 0..%d with the lid, N = 80, worst elevation: %.1e" % (ORDERS - 1, 2 * ORDERS - 1, tail))
    require(tail < 1e-10, "orders with the lid")

    def extrapolate(t, n):            # linear in 1/N through n and 2 n, componentwise
        return [2 * t[2 * n][i] - t[n][i] for i in range(3)]

    print("\nExtrapolated in 1/N, linear through N and 2N, '%s':" % DEFAULT)
    ext = {n: extrapolate(table, n) for n in (20, 40, 160)}
    for i, lab in enumerate(labels):
        print("  %-13s from 20, 40: %.6f   from 40, 80: %.6f   from 160, 320: %.6f" %
              (lab, abs(ext[20][i]), abs(ext[40][i]), abs(ext[160][i])))
    best = ext[160]

    print("\nOther projections, extrapolated from N = 160, 320, relative difference to '%s':" % DEFAULT)
    spread = [abs(ext[40][i] - best[i]) / abs(best[i]) for i in range(3)]
    for proj in PROJECTIONS[1:]:
        t = {n: values(Problem(a, b, h, k, n), proj) for n in (160, 320)}
        e = extrapolate(t, 160)
        d = [abs(e[i] - best[i]) / abs(best[i]) for i in range(3)]
        spread = [max(x, y) for x, y in zip(spread, d)]
        print("  %-9s " % proj + "  ".join("%s %.6f (%.1e)" % (lab, abs(e[i]), d[i]) for i, lab in enumerate(labels)))

    print("\nReference (projection '%s', extrapolated from N = 160, 320). Spread: the largest relative difference to the"
          " extrapolation from N = 40, 80 and to the other projections." % DEFAULT)
    for i, lab in enumerate(labels):
        print("  %-13s %s  re %.8f  im %.8f  spread %.1e" % (lab, fmt(best[i]), best[i].real, best[i].imag, spread[i]))
    print("\nLid over no lid (MacCamy-Fuchs): " +
          ", ".join("%s %.4f" % (lab, abs(x) / abs(y)) for lab, x, y in zip(labels, best, ref)))
    print("\nSelf-checks " + ("passed" if not failed else "FAILED: " + ", ".join(failed)))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
