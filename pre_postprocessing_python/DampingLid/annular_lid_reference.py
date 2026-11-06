"""Semi-analytical reference for the dissipative damping lid: a bottom-mounted vertical cylinder with an annular lid around it.

Linear potential flow in water of depth h, time factor exp(-i omega t), nu = omega**2/g. A cylinder of radius a stands on the sea
bed and pierces the free surface. On the annulus a < r < b of z = 0 (the lid) the pressure p/rho = E dphi/dz acts (the WAMIT
damping lid, E in m/s), which with the kinematic condition gives

    dphi/dz = K phi,   K = nu g/(g - i omega E)

E = 0 is open water, E -> infinity a rigid lid (dphi/dz = 0). Outside it, r > b, dphi/dz = nu phi. A plane wave of unit amplitude
comes along +x. The problem is solved by eigenfunction matching, one azimuthal order m at a time, with no boundary element and
no Green function:

- lid region a < r < b: vertical modes L_j = cos(p_j (z + h))/cos(p_j h), j = 0..N-1, with p tan(p h) = -K, so that
  dL/dz = K L on z = 0 and dL/dz = 0 on z = -h. In the wave number form k = i p these are the roots of k tanh(kh) = K. They are
  complex for E > 0. Each is found by Newton iteration on k sinh(kh) - K cosh(kh) (the same roots as k tanh(kh) - K, without
  the poles), started from the open-water root (k_0, or i k_n) and continued in small steps from E = 0 to E, so that each
  root follows its open-water branch (lid_roots). The residual of k tanh(kh) - K is printed. Radial functions
  R_j = I_m(p r) K_m'(p a) - K_m(p r) I_m'(p a), so dphi/dr = 0 at r = a (Re p >= 0 is chosen; R does not depend on the
  branch). For the rigid lid (E = infinity) p_0 = 0 and R = r**m + a**(2m) r**-m.
- outer region r > b: the propagating mode Z_0 = cosh k(z + h)/cosh kh with J_m(kr) (incident) and H_m^(1)(kr) (scattered),
  and the evanescent modes Z_n = cos k_n(z + h)/cos k_n h with K_m(k_n r), n = 1..N-1, k_n tan(k_n h) = -nu.
- matching at r = b: phi and dphi/dr continuous on -h < z < 0. Both equations are projected on the N lid modes L_j ('lid',
  the reference values). The lid modes are orthogonal in the bilinear product int L_i L_j dz (no complex conjugate), as they
  solve the same Sturm-Liouville problem with the same complex K. The script also solves with both projected on the outer
  modes ('outer') and with the two mixed choices, and prints the spread between them.

All inner products are closed form: int_0^h cos(p s) cos(q s) ds = (S(p - q) + S(p + q))/2, S(x) = sin(x h)/x (S(0) = h).

Time factor and output normalisation, from src/PatclVelct.f90 (VINP), src/PotentWavForce.f90 (EFORCE) and
src/PressureElevation.f90 (CalElevation, WamitNondimens): HAMS solves with exp(-i omega t) and an incident potential whose
elevation at the origin is -i A. ExcForce.3 writes (-Im X, -Re X)/(rho g A) and PressureElevation.6p writes (-Im eta, -Re eta)/A,
which is conj(X/(rho g eta_I(0))) and conj(eta/eta_I(0)): the WAMIT convention, exp(+i omega t) and phases relative to the
incident elevation at the origin. Here psi is the potential scaled so that the incident wave is exp(ikx) Z_0(z), so that
eta/eta_I(0) = psi(z = 0) off the lid, and the force on the body over rho g eta_I(0) is -int psi n_x dS (n out of the body).
The script prints conj of both. On the lid eta = i omega phi/(g - i omega E) = psi K/nu, and dphi/dz = -i omega eta.

Energy. With phi = (g/(i omega)) psi A, the mean power into the fluid inside a cylinder of radius R is
1/2 rho omega int Im(phi conj(dphi/dr)) dS and the power the lid absorbs is 1/2 rho E int abs(dphi/dz)**2 dA
= 1/2 rho omega**2 E int abs(eta)**2 dA. Both are printed for A = 1 m, rho = 1000 kg/m**3, in W.

Self-checks; the script exits with 1 if one fails. What each covers, and what it does not:

1. b = a: X1 and the elevations equal MacCamy-Fuchs to 1e-8, for every projection. It covers the outer modes, the incident
   term, the force and elevation sums and the output convention, not the lid modes.
2. E = 0 with the lid b = 4: the lid region is open water, so X1 and the elevations must equal MacCamy-Fuchs, at every N and
   projection. It runs the lid roots, M, the lid radial functions, gamma and the lid force sum at the real lid, with a lid
   condition of the full strength nu. It does not test a complex K.
3. Roots: for E = 0.2, 1, 1e6 and the conjugate lid at E = 1, the residual of k tanh(kh) - K of every root, that the roots are
   distinct, and for E = 1e6 that root j is within 1e-3/h of the rigid root j pi/h (the branch is followed).
4. E = 1e6 against the rigid lid E = infinity (p_0 = 0, its own radial functions), N = 80: the difference is of the order of
   K/nu ~ g/(omega E) and must be below 1e-4. It covers the complex roots and radial functions at small K against a separate
   code path.
5. Pointwise matching at r = b, at the depths ZS, orders 0, 1, 2, 5, E = 1: the mismatch of phi and of dphi/dr between the
   outer and the lid expansions, each evaluated from its own modes, falls with N (N = 20, 80, 320).
6. Energy balance in the diffraction form for E = 0.2 and 1: the net power into a large cylinder (r = 1000 m, the incident
   plus scattered propagating wave) equals the power absorbed on the lid, 1/2 rho E int abs(dphi/dz)**2 dA, from the
   lid-region solution, to 1e-6 relative. Lossless lids (E = 0, E = infinity) give abs(1 + 2 a_0/(eps_m i**m H_m(kb))) = 1.
   With the conjugate lid K' = nu g/(g + i omega E), the same power must come out with the opposite sign (the lid then
   makes energy): that check shows the energy balance tells the two signs apart.
7. Orders: 0..29 against 0..59 at the two points, E = 1, below 1e-10.

The spread between the four projections and between extrapolations from different N (printed) bounds the truncation error of
the reference values; it is not a check of the formulation.

Run: python -I annular_lid_reference.py [--epsilon E] [--conjugate]. Without --epsilon it prints the references for
E = 0.2, 1 and 1e6; --conjugate uses K' = nu g/(g + i omega E) (the sign check of the HAMS-MREL tests). It takes about 7
minutes (measured 7 min 12 s on the laptop), most of it in the lid power integrals.
"""
import argparse
import math
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import h1vp, hankel1, ive, jv, jvp, kve, roots_legendre

A_RADIUS = 2.0   # cylinder radius a
B_RADIUS = 4.0   # outer radius of the lid b
DEPTH = 10.0     # h
KH = 1.5         # k h; HAMS is given the wave number k = 0.15 directly (Input_frequency_type 2)
G = 9.80665      # g of src/WavDynMods.f90
RHO = 1000.0     # rho of src/WavDynMods.f90, for the powers only
ORDERS = 30      # azimuthal orders m = 0..ORDERS-1 for the elevation and the powers (the force needs only m = 1)
POINTS = [(0.0, 8.0), (8.0, 0.0)]
ZS = [-1.0, -2.5, -5.0, -7.5, -9.0]   # depths of the pointwise matching check
PROJECTIONS = ("lid", "outer", "phi_lid", "phi_outer")
DEFAULT = PROJECTIONS[0]   # the projection of the reference values
EPSILONS = (0.2, 1.0, 1.0e6)
STEPS = 2000     # continuation steps from E = 0 to E
R_FAR = 1000.0   # radius of the control cylinder of the energy balance


def evanescent_roots(nu, h, count):
    """The roots k_n, n = 1..count, of k tan(k h) = -nu, k_n h in ((n - 1/2) pi, n pi). Solved as y = n pi - atan(nu h/y),
    y = k_n h, which has no pole in the bracket."""
    out = []
    for n in range(1, count + 1):
        y = brentq(lambda y: y - n * math.pi + math.atan(nu * h / y), (n - 0.5) * math.pi, n * math.pi, xtol=1e-15, rtol=1e-15,
                   maxiter=200)
        out.append(y / h)
    return np.array(out)


def lid_K(nu, omega, eps, sign=1):
    """K = nu g/(g - sign i omega E); sign = -1 is the conjugate lid. E = inf gives 0."""
    if math.isinf(eps):
        return 0j
    return nu * G / (G - sign * 1j * omega * eps)


def lid_roots(k, kn, h, nu, omega, eps, sign=1):
    """p_j, j = 0..N-1, of the lid modes cos(p (z + h)), p = -i k with k tanh(kh) = K, Re p >= 0, and the largest residual of
    k tanh(kh) - K over the roots. Newton on f(k) = k sinh(kh) - K cosh(kh), started from the open-water roots k_0 and i k_n
    and continued from E = 0 to E in STEPS steps of the angle atan(omega E'/g), along which K = nu cos(t) exp(sign i t)."""
    n = len(kn) + 1
    if eps == 0.0:
        return np.concatenate(([-1j * k], kn)).astype(complex), 0.0
    if math.isinf(eps):
        return (np.arange(n) * math.pi / h).astype(complex), 0.0
    kk = np.concatenate(([k], 1j * kn)).astype(complex)
    t_end = math.atan(omega * eps / G)
    for t in np.linspace(0.0, t_end, STEPS + 1)[1:]:
        K = nu * math.cos(t) * complex(math.cos(t), sign * math.sin(t))
        for _ in range(50):
            ch, sh = np.cosh(kk * h), np.sinh(kk * h)
            f = kk * sh - K * ch
            df = sh + kk * h * ch - K * h * sh
            step = f / df
            kk = kk - step
            if np.max(np.abs(step) / np.maximum(np.abs(kk), 1e-300)) < 1e-15:
                break
    K = lid_K(nu, omega, eps, sign)
    residual = np.max(np.abs(kk * np.tanh(kk * h) - K))
    p = -1j * kk
    p = np.where(p.real < 0, -p, p)
    return p, residual


def S(x, h):
    """sin(x h)/x, h at x = 0, for complex x."""
    x = np.asarray(x, dtype=complex)
    small = np.abs(x) < 1e-8
    safe = np.where(small, 1.0, x)
    return np.where(small, h * (1 - (x * h) ** 2 / 6), np.sin(safe * h) / safe)


class Problem:
    """One geometry, frequency and lid. b == a gives the bare cylinder. eps in m/s, math.inf for the rigid lid; sign = -1 for
    the conjugate lid condition."""

    def __init__(self, a, b, h, k, n, eps=0.0, sign=1):
        self.a, self.b, self.h, self.k, self.n, self.eps, self.sign = a, b, h, k, n, eps, sign
        self.nu = k * math.tanh(k * h)
        self.omega = math.sqrt(G * self.nu)
        self.kn = evanescent_roots(self.nu, h, n - 1)
        self.K = lid_K(self.nu, self.omega, eps, sign)
        self.p, self.root_residual = lid_roots(k, self.kn, h, self.nu, self.omega, eps, sign)
        self.rigid0 = math.isinf(eps)                                         # p_0 = 0 exactly
        q = np.concatenate(([-1j * k], self.kn)).astype(complex)              # outer modes cos(q (z + h))/cos(q h)
        p = self.p
        cq, cp = np.cos(q * h), np.cos(p * h)
        self.M = 0.5 * (S(p[None, :] - q[:, None], h) + S(p[None, :] + q[:, None], h)) / (cq[:, None] * cp[None, :])  # M[n, j]
        nrm = np.empty(n)
        nrm[0] = (h / 2 + math.sinh(2 * k * h) / (4 * k)) / math.cosh(k * h) ** 2
        nrm[1:] = (h / 2 + np.sin(2 * self.kn * h) / (4 * self.kn)) / np.cos(self.kn * h) ** 2
        self.N = nrm                                                           # int Z_n**2
        self.P = 0.5 * (h + S(2 * p, h)) / cp ** 2                             # int L_j**2 (bilinear)
        self.zint_lid = S(p, h) / cp                                           # int L_j dz
        self.zint_out = np.concatenate(([math.tanh(k * h) / k], np.tan(self.kn * h) / self.kn))   # int Z_n dz
        self.cache = {}

    # radial functions, each normalised to 1 at r = b
    def lid_radial(self, m, r):
        """R_j(r)/R_j(b) and R_j'(b)/R_j(b) for every j, from exponentially scaled Bessel functions. r may be an array (then
        the result is [len(r), N])."""
        p, a, b = self.p, self.a, self.b
        r = np.asarray(r, dtype=float)
        scalar = r.ndim == 0
        r = np.atleast_1d(r)[:, None]
        pr = p.real                                                            # >= 0
        ipe_a = 0.5 * (ive(m - 1, p * a) + ive(m + 1, p * a))                  # I_m'(p a) exp(-Re(p) a)
        kpe_a = -0.5 * (kve(m - 1, p * a) + kve(m + 1, p * a))                 # K_m'(p a) exp(p a)

        def s(x):                                                              # R(x) exp(p a - Re(p) x)
            return ive(m, p * x) * kpe_a - np.exp(-(p + pr) * (x - a)) * kve(m, p * x) * ipe_a

        def sp(x):                                                             # R'(x) exp(p a - Re(p) x)/p
            ip = 0.5 * (ive(m - 1, p * x) + ive(m + 1, p * x))
            kp = -0.5 * (kve(m - 1, p * x) + kve(m + 1, p * x))
            return ip * kpe_a - np.exp(-(p + pr) * (x - a)) * kp * ipe_a

        with np.errstate(invalid="ignore", divide="ignore"):
            gamma = p * sp(b) / s(b)
            sb = s(b)
            if r.shape[0] == 1:
                ratio = np.exp(pr * (r - b)) * s(r) / sb
            else:                                                              # skip the modes below exp(-60) at r
                ratio = np.zeros((r.shape[0], len(p)), dtype=complex)
                for i, x in enumerate(r[:, 0]):
                    use = pr * (b - x) < 60.0
                    pu, pru = p[use], pr[use]
                    su = (ive(m, pu * x) * kpe_a[use]
                          - np.exp(-(pu + pru) * (x - a)) * kve(m, pu * x) * ipe_a[use])
                    ratio[i, use] = np.exp(pru * (x - b)) * su / sb[use]
        if self.rigid0:                                                        # p_0 = 0: r**m + a**(2m) r**-m
            q = (a / b) ** (2 * m)
            ratio[:, 0] = ((r[:, 0] / b) ** m + (a * a / (r[:, 0] * b)) ** m) / (1 + q)
            gamma[0] = m / b * (1 - q) / (1 + q)
        return (ratio[0] if scalar else ratio), gamma

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
        lid modes (radial functions 1 at r = b). Cached per order and projection."""
        key = (m, projection)
        if key not in self.cache:
            self.cache[key] = self._solve(m, projection)
        return self.cache[key]

    def _solve(self, m, projection):
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
            W = (M * (gamma / P)[None, :]) @ M.T                               # W[p, n] = sum_q gamma_q M_pq M_nq / P_q
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

    def lid_elevation(self, r, th, projection=DEFAULT, orders=ORDERS):
        """Elevation over eta_I(0) at (r cos th, r sin th, 0), a <= r <= b, on the lid: psi K/nu, HAMS convention."""
        r, th = np.asarray(r, dtype=float), np.asarray(th, dtype=float)
        total = np.zeros(r.shape, dtype=complex)
        for m in range(orders):
            _, c, _ = self.solve(m, projection)
            total += (self.lid_radial(m, r.ravel())[0] @ c).reshape(r.shape) * np.cos(m * th)
        return np.conj(total * self.K / self.nu)

    def powers(self, projection=DEFAULT, orders=ORDERS):
        """(power into the cylinder r = R_FAR, power absorbed on the lid), W for an incident amplitude of 1 m. The first from
        the propagating outer mode only (the evanescent ones are below 1e-100 at R_FAR), the second from the lid modes."""
        k, b = self.k, self.b
        c2 = (G / self.omega) ** 2                                             # abs(phi/psi)**2
        x, w = radial_rule(self.a, b)
        far = 0.0
        lid = 0.0
        for m in range(orders):
            a, c, _ = self.solve(m, projection)
            s = (1.0 if m == 0 else 2.0) * 1j ** m
            wm = 2 * math.pi if m == 0 else math.pi
            f = s * jv(m, k * R_FAR) + a[0] * hankel1(m, k * R_FAR) / hankel1(m, k * b)
            df = k * (s * jvp(m, k * R_FAR) + a[0] * h1vp(m, k * R_FAR) / hankel1(m, k * b))
            far += wm * R_FAR * self.N[0] * (f * np.conj(df)).imag
            if self.b > self.a:
                g = self.lid_radial(m, x)[0] @ c                                # psi(r, z = 0) of order m
                lid += wm * np.sum(w * x * np.abs(g) ** 2)
        far *= 0.5 * RHO * self.omega * c2
        e = 0.0 if math.isinf(self.eps) else self.eps
        lid *= 0.5 * RHO * e * abs(self.K) ** 2 * c2
        return far, lid

    def matching_residual(self, m, zs, projection=DEFAULT):
        """Pointwise mismatch at r = b of phi and of dphi/dr between the outer and the lid expansions of order m, at the depths
        zs, each over the largest outer value there."""
        a, c, inc = self.solve(m, projection)
        beta = self.outer_beta(m)
        gamma = self.lid_radial(m, self.b)[1]
        z = np.asarray(zs)[:, None]
        zo = np.concatenate([np.cosh(self.k * (z + self.h)) / math.cosh(self.k * self.h),
                             np.cos(self.kn[None, :] * (z + self.h)) / np.cos(self.kn * self.h)[None, :]], axis=1)
        lo = np.cos(self.p[None, :] * (z + self.h)) / np.cos(self.p * self.h)[None, :]
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


def radial_rule(a, b, per=12, levels=24):
    """Gauss-Legendre nodes and weights on [a, b], graded geometrically towards both ends (the lid modes of large p vary
    over a length 1/p next to r = b)."""
    half = 0.5 * (b - a)
    d = half * np.concatenate(([0.0], np.geomspace(1e-7, 1.0, levels)))
    t, wt = roots_legendre(per)
    xs, ws = [], []
    for lo, hi in zip(d[:-1], d[1:]):
        mid, rad = 0.5 * (lo + hi), 0.5 * (hi - lo)
        for side in (a + mid + rad * t, b - mid - rad * t):
            xs.append(side)
            ws.append(rad * wt)
    return np.concatenate(xs), np.concatenate(ws)


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


LABELS = ["surge X1", "eta(0, 8, 0)", "eta(8, 0, 0)"]
NS = (20, 40, 80, 160, 320)


def extrapolate(t, n, count):        # linear in 1/N through n and 2 n, componentwise
    return [2 * t[2 * n][i] - t[n][i] for i in range(count)]


def reference(eps, sign, mf, require):
    """The convergence table, the extrapolations, the projections and the powers for one lid. Returns the reference values
    (X1, two elevations) and the extrapolated lid power."""
    a, b, h, k = A_RADIUS, B_RADIUS, DEPTH, KH / DEPTH
    name = "E = %g m/s%s" % (eps, ", conjugate K' = nu g/(g + i omega E)" if sign < 0 else "")
    print("\nLid a < r < b, %s, projection '%s', convergence in N (powers in W for A = 1 m):" % (name, DEFAULT))
    table, power = {}, {}
    for n in NS:
        prob = Problem(a, b, h, k, n, eps, sign)
        table[n] = values(prob)
        far, lid = prob.powers()
        power[n] = [far, lid]
        print("  N = %3d  root residual %.1e  power into r = %g: %.8e  lid power: %.8e  balance %.1e" %
              (n, prob.root_residual, R_FAR, far, lid, abs(far - lid) / abs(lid)))
        for lab, v in zip(LABELS, table[n]):
            print("      %-13s %s" % (lab, fmt(v)))
    ext = {n: extrapolate(table, n, 3) for n in (20, 40, 160)}
    pext = {n: extrapolate(power, n, 2) for n in (40, 160)}
    print("  Extrapolated in 1/N, linear through N and 2N:")
    for i, lab in enumerate(LABELS):
        print("    %-13s from 20, 40: %.6f   from 40, 80: %.6f   from 160, 320: %.6f" %
              (lab, abs(ext[20][i]), abs(ext[40][i]), abs(ext[160][i])))
    print("    lid power     from 40, 80: %.6e   from 160, 320: %.6e" % (pext[40][1], pext[160][1]))
    best = ext[160]
    spread = [abs(ext[40][i] - best[i]) / abs(best[i]) for i in range(3)]
    print("  Other projections, extrapolated from N = 160, 320, relative difference to '%s':" % DEFAULT)
    for proj in PROJECTIONS[1:]:
        t = {n: values(Problem(a, b, h, k, n, eps, sign), proj) for n in (160, 320)}
        e = extrapolate(t, 160, 3)
        d = [abs(e[i] - best[i]) / abs(best[i]) for i in range(3)]
        spread = [max(x, y) for x, y in zip(spread, d)]
        print("    %-9s " % proj + "  ".join("%s %.6f (%.1e)" % (lab, abs(e[i]), d[i]) for i, lab in enumerate(LABELS)))
    print("  Reference, %s (extrapolated from N = 160, 320). Spread: the largest relative difference to the extrapolation"
          " from N = 40, 80 and to the other projections." % name)
    for i, lab in enumerate(LABELS):
        print("    %-13s %s  re %.8f  im %.8f  spread %.1e" % (lab, fmt(best[i]), best[i].real, best[i].imag, spread[i]))
    pspread = abs(pext[40][1] - pext[160][1]) / abs(pext[160][1]) if pext[160][1] != 0 else 0.0
    print("    lid power     %.6e W  spread %.1e (from 40, 80)" % (pext[160][1], pspread))
    print("  Over no lid (MacCamy-Fuchs): " + ", ".join("%s %.4f" % (lab, abs(x) / abs(y)) for lab, x, y in zip(LABELS, best, mf)))
    return best, pext[160][1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--epsilon", type=float, default=None, help="damping coefficient E of the lid, m/s (default: 0.2, 1, 1e6)")
    parser.add_argument("--conjugate", action="store_true", help="use the conjugate lid K' = nu g/(g + i omega E)")
    args = parser.parse_args(argv)
    sign = -1 if args.conjugate else 1

    a, b, h = A_RADIUS, B_RADIUS, DEPTH
    k = KH / h
    nu = k * math.tanh(k * h)
    omega = math.sqrt(G * nu)
    print("a = %g, b = %g, h = %g, k = %.6f (kh = %g), omega = %.10f rad/s, nu = %.10f" % (a, b, h, k, KH, omega, nu))
    failed = []

    def require(cond, what):
        if not cond:
            failed.append(what)

    print("\nCheck 1, b = a against MacCamy-Fuchs (outer modes only):")
    mf = [maccamy_fuchs_surge(a, h, k)] + [maccamy_fuchs_elevation(a, k, x, y) for x, y in POINTS]
    for proj in PROJECTIONS:
        got = values(Problem(a, a, h, k, 20), proj)
        err = max(abs(g - r) / abs(r) for g, r in zip(got, mf))
        require(err < 1e-8, "b = a, " + proj)
        print("  %-9s worst relative difference %.2e" % (proj, err))
    for lab, r in zip(LABELS, mf):
        print("  MacCamy-Fuchs %-13s %s  re %.8f  im %.8f" % (lab, fmt(r), r.real, r.imag))

    print("\nCheck 2, E = 0 with the lid b = %g (lid region = open water) against MacCamy-Fuchs, through the lid modes:" % b)
    for proj in PROJECTIONS:
        errs = []
        for n in (20, 80, 320):
            got = values(Problem(a, b, h, k, n, 0.0), proj)
            errs.append(max(abs(g - r) / abs(r) for g, r in zip(got, mf)))
        require(max(errs) < 1e-8, "E = 0, " + proj)
        print("  %-9s worst relative difference at N = 20, 80, 320: %s" % (proj, ", ".join("%.1e" % e for e in errs)))

    print("\nCheck 3, lid roots (N = 320): residual of k tanh(kh) - K, smallest distance between two roots, K:")
    for eps, sg in [(0.2, 1), (1.0, 1), (1.0e6, 1), (1.0, -1)]:
        prob = Problem(a, b, h, k, 320, eps, sg)
        p = prob.p
        gap = np.min(np.abs(p[:, None] - p[None, :]) + np.eye(len(p)) * 1e9)
        line = "  E = %-7g %s residual %.1e  smallest gap %.3e  K = %+.6e %+.6ei  p_0 h = %+.6f %+.6fi" % (
            eps, "conj" if sg < 0 else "    ", prob.root_residual, gap, prob.K.real, prob.K.imag, (p[0] * h).real,
            (p[0] * h).imag)
        require(prob.root_residual < 1e-10 and gap > 1e-3, "roots, E = %g, sign %d" % (eps, sg))
        if eps == 1.0e6:
            branch = np.max(np.abs(p[1:] - np.arange(1, len(p)) * math.pi / h)) * h
            line += "  max abs(p_j h - j pi), j >= 1: %.1e" % branch
            require(branch < 1e-3, "branch, E = 1e6")
        print(line)

    print("\nCheck 4, E = 1e6 against the rigid lid E = infinity (p_0 = 0), N = 80:")
    rig = values(Problem(a, b, h, k, 80, math.inf))
    big = values(Problem(a, b, h, k, 80, 1.0e6))
    d = max(abs(x - y) / abs(y) for x, y in zip(big, rig))
    print("  worst relative difference %.1e (g/(omega E) = %.1e)" % (d, G / (omega * 1e6)))
    require(d < 1e-4, "E = 1e6 against rigid")

    print("\nCheck 5, pointwise matching at r = b, z = %s, E = 1: largest mismatch of phi and of dphi/dr over m = 0, 1, 2, 5,"
          "\n  each relative to the largest outer value:" % ZS)
    for proj in PROJECTIONS:
        rows = []
        for n in (20, 80, 320):
            res = [Problem(a, b, h, k, n, 1.0).matching_residual(m, ZS, proj) for m in (0, 1, 2, 5)]
            rows.append((max(r[0] for r in res), max(r[1] for r in res)))
        print("  %-9s " % proj + "  ".join("N = %d: phi %.1e, dphi/dr %.1e" % (n, r[0], r[1]) for n, r in zip((20, 80, 320), rows)))
        if proj == DEFAULT:
            require(rows[0][0] > rows[1][0] > rows[2][0] and rows[0][1] > rows[1][1] > rows[2][1], "matching, " + proj)

    print("\nCheck 6, energy balance (diffraction form): power into the cylinder r = %g equals the lid power, 1/2 rho E int"
          " abs(dphi/dz)**2 dA, W for A = 1 m:" % R_FAR)
    for eps, sg in [(0.2, 1), (1.0, 1), (1.0, -1)]:
        for n in (80, 320):
            far, lid = Problem(a, b, h, k, n, eps, sg).powers()
            rel = abs(far - lid) / abs(lid)
            print("  E = %-4g %s N = %3d  into r = %g: %+.10e  lid: %+.10e  relative residual %.1e" %
                  (eps, "conj" if sg < 0 else "    ", n, R_FAR, far, lid, rel))
            if sg > 0:
                require(rel < 1e-6, "energy balance, E = %g, N = %d" % (eps, n))
            else:
                require(far < 0 < lid, "conjugate lid makes energy, N = %d" % n)
    for eps in (0.0, math.inf):
        defect = max(abs(Problem(a, b, h, k, 80, eps).energy_defect(m)) for m in range(ORDERS))
        print("  E = %g, lossless: max abs(abs(1 + 2 a_0/(eps_m i**m H_m(kb))) - 1) %.1e (N = 80)" % (eps, defect))
        require(defect < 1e-8, "lossless, E = %g" % eps)

    tail = max(abs(x - y) for x, y in zip(values(Problem(a, b, h, k, 80, 1.0))[1:],
                                          values(Problem(a, b, h, k, 80, 1.0), orders=2 * ORDERS)[1:]))
    print("\nCheck 7, orders 0..%d against 0..%d, E = 1, N = 80, worst elevation: %.1e" % (ORDERS - 1, 2 * ORDERS - 1, tail))
    require(tail < 1e-10, "orders with the lid")

    for eps in ((args.epsilon,) if args.epsilon is not None else EPSILONS):
        reference(eps, sign, mf, require)

    print("\nSelf-checks " + ("passed" if not failed else "FAILED: " + ", ".join(failed)))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
