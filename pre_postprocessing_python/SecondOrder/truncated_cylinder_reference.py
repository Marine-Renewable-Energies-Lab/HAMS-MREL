"""Semi-analytical reference for a truncated vertical cylinder, by eigenfunction matching, to validate HAMS-MREL's mean drift.

Problem (design: docs/superpowers/specs/2026-11-06-truncated-cylinder-drift-reference-design.md, sections 2 and 3). A vertical
cylinder of radius a and draft d, sharp bilge, in water of depth h, with the gap g = h - d below it. Linear potential flow,
time factor exp(-i omega t), nu = omega**2/G. G and RHO are make_mesh.G and make_mesh.RHO, the values of src/WavDynMods.f90.
z is up, z = 0 is the still waterline, z = -h the sea bed. Each azimuthal order m is solved on its own.

Vertical modes. Exterior, r > a, -h < z < 0:

    Z_0(z) = cosh k (z + h)/cosh k h,       k tanh(k h) = nu, k > 0,
    Z_n(z) = cos k_n (z + h)/cos k_n h,     k_n tan(k_n h) = -nu, k_n h in ((n - 1/2) pi, n pi), n = 1 ... N,

so that Z_n(0) = 1. Interior, r < a, -h < z < -d: cos(lambda_n (z + h)), lambda_n = n pi/g, n = 0 ... N, not normalised.
k is drift_convergence.wave_number (Newton). k_n is solved as y = n pi - atan(nu h/y), y = k_n h, by brentq; this form has
no pole in the bracket (as evanescent_roots in annular_lid_reference.py on feature/damping-lid).

Corner basis. The fluid corner at (a, -d) is a 270 degree wedge, so the radial velocity on the gap behaves like
(-d - z)**(-1/3) (design section 3). Reflecting the gap about the sea bed, where dphi/dz = 0, gives t = (z + h)/g in [-1, 1],
with the corner at t = +-1 and the bed at t = 0, and the even basis

    u_p(t) = (1 - t**2)**(-1/3) C_2p^(1/6)(t),   p = 0 ... P,

with C^(1/6) the Gegenbauer polynomial of scipy.special.eval_gegenbauer. The weight (1 - t**2)**(-1/3) is the Gegenbauer
weight (1 - t**2)**(lambda - 1/2) for lambda = 1/6, in which the C_2p^(1/6) are orthogonal, and the Jacobi weight for
alpha = beta = -1/3. gegenbauer_basis gives Gauss-Jacobi nodes and weights for it (scipy.special.roots_jacobi).

Inner products. For f even about the bed, the integral over the gap of u_p f dz is (g/2) int_-1^1 u_p(t) f dt. The modes are
cos(a t) or cosh(a t) in t, and Gegenbauer's integral (DLMF 18.17.17; x -> i x for cosh) gives the closed form

    int_-1^1 (1 - t**2)**(-1/3) C_2p^(1/6)(t) cos(a t) dt  = c_p (-1)**p J_(2p+1/6)(a) a**(-1/6),
    int_-1^1 (1 - t**2)**(-1/3) C_2p^(1/6)(t) cosh(a t) dt = c_p I_(2p+1/6)(a) a**(-1/6),
    c_p = pi 2**(5/6) Gamma(2p + 1/3)/((2p)! Gamma(1/6)),

with the limit a -> 0 equal to c_0 2**(-1/6)/Gamma(7/6) for p = 0 and zero otherwise. inner_exterior and inner_interior use
it. The propagating column uses ive: I(kg)/cosh(kh) = ive(kg) 2 exp(kg - kh)/(1 + exp(-2kh)), which does not overflow at
large kh. The Gauss-Jacobi rule of the basis gives the same integrals on the reflected interval; the checks compare both.

HAMS output normalisation: later tasks fill this in from src/PrintOutput.f90, src/PotentWavForce.f90 and
src/SecondOrderOutput.f90, with the routine each rule comes from. Nothing here depends on it yet.

Self-checks: python -I truncated_cylinder_reference.py check [NAME ...] runs the named checks, or all of them without NAME.
The script exits with 1 if one fails. What each covers, and what it does not:

- roots: at each omega of drift_convergence.WIDE_FREQUENCIES, h = 100, N = CHECK_N: abs(k tanh kh - nu) and
  abs(k_n tan k_n h + nu) are below 1e-12 nu, and every k_n lies strictly inside ((n - 1/2) pi/h, n pi/h). It covers
  wave_numbers, not the use of the roots. The residual of k_n is limited by the rounding of k_n h, about (k_n h)**2 eps/h,
  so it holds to 1e-12 nu only for the first few tens of roots; the check uses the first CHECK_N.
- inner_products: default cylinder, p <= CHECK_P, n <= CHECK_N, the exterior at every WIDE_FREQUENCIES. inner_exterior and
  inner_interior equal scipy.integrate.quad of the same integrands on the unreflected gap, with the corner weight handled by
  weight="alg", to 1e-10 (max abs difference over max abs quad value, per matrix). The Gauss-Jacobi sums of the basis are
  checked against the same quad values, which tests the nodes, the weights and the reflection. It covers the closed form,
  the reflection and the mode normalisation. It does not test that u_p is the right basis for the corner; the matching
  check of a later task does.
"""
import argparse
import math
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer, gammaln, ive, jv, roots_jacobi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drift_convergence  # noqa: E402
import make_mesh  # noqa: E402

LAMBDA = 1.0 / 6.0     # Gegenbauer parameter of the corner basis
ALPHA = LAMBDA - 0.5   # Jacobi parameter alpha = beta = -1/3 of its weight
GAUSS_JACOBI_NODES = 96   # default size of the Gauss-Jacobi rule: exact for polynomials of degree < 192

# Resolution of the root and inner product checks
CHECK_P = 10
CHECK_N = 20


@dataclass
class Cylinder:
    """Truncated vertical cylinder, lengths in m. Defaults: the study's sharp hull (HULLS["sharp"] in drift_convergence.py)."""
    radius: float = 10.0
    draft: float = 10.0
    depth: float = 100.0

    @property
    def gap(self):
        return self.depth - self.draft


def wave_numbers(omega, depth, count):
    """(k, kn): k > 0 with k tanh(k h) = nu, and kn, the first count roots of kn tan(kn h) = -nu, increasing."""
    nu = omega ** 2 / make_mesh.G
    k = drift_convergence.wave_number(omega, depth)
    kn = np.empty(count)
    for n in range(1, count + 1):
        y = brentq(lambda y: y - n * math.pi + math.atan(nu * depth / y), (n - 0.5) * math.pi, n * math.pi,
                   xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=200)
        kn[n - 1] = y / depth
    return k, kn


@dataclass
class Basis:
    """The corner basis u_p(t) = (1 - t**2)**(-1/3) C_2p^(1/6)(t), p = 0 ... p_max, with a Gauss-Jacobi rule for its weight:
    int_-1^1 (1 - t**2)**(-1/3) f(t) dt ~ sum(weights * f(nodes)). poly[p, j] = C_2p^(1/6)(nodes[j])."""
    p_max: int
    nodes: np.ndarray
    weights: np.ndarray
    poly: np.ndarray

    def values(self, t):
        """u_p(t), shape (p_max + 1,) + shape(t), for abs(t) < 1."""
        t = np.asarray(t, dtype=float)
        return np.array([eval_gegenbauer(2 * p, LAMBDA, t) for p in range(self.p_max + 1)]) * (1.0 - t ** 2) ** ALPHA


def gegenbauer_basis(p_max, nodes=GAUSS_JACOBI_NODES):
    """The basis u_p, p = 0 ... p_max, with a Gauss-Jacobi rule of `nodes` points for alpha = beta = -1/3."""
    t, w = roots_jacobi(nodes, ALPHA, ALPHA)
    poly = np.array([eval_gegenbauer(2 * p, LAMBDA, t) for p in range(p_max + 1)])
    return Basis(p_max, t, w, poly)


def _gegenbauer_constant(p):
    """c_p = pi 2**(5/6) Gamma(2p + 1/3)/((2p)! Gamma(1/6)), elementwise in p."""
    p = np.asarray(p, dtype=float)
    return math.pi * 2.0 ** (1.0 - LAMBDA) * np.exp(gammaln(2 * p + 2 * LAMBDA) - gammaln(2 * p + 1) - gammaln(LAMBDA))


def _cos_integral(p_max, a):
    """M[p, j] = int_-1^1 (1 - t**2)**(-1/3) C_2p^(1/6)(t) cos(a_j t) dt, for a_j >= 0."""
    p = np.arange(p_max + 1)[:, None]
    a = np.asarray(a, dtype=float)[None, :]
    safe = np.where(a > 0, a, 1.0)
    out = _gegenbauer_constant(p) * (-1.0) ** p * jv(2 * p + LAMBDA, safe) * safe ** (-LAMBDA)
    at_zero = np.where(p == 0, _gegenbauer_constant(0) * 2.0 ** (-LAMBDA) / math.gamma(1 + LAMBDA), 0.0)
    return np.where(a > 0, out, at_zero)


def inner_exterior(cyl, k, kn, basis):
    """E[p, n] = int_-h^-d u_p(t(z)) Z_n(z) dz, column 0 the propagating Z_0, columns 1 ... len(kn) the evanescent Z_n;
    shape (p_max + 1, len(kn) + 1)."""
    g, h = cyl.gap, cyl.depth
    p = np.arange(basis.p_max + 1)
    a = k * g
    ratio = 2.0 * math.exp(a - k * h) / (1.0 + math.exp(-2.0 * k * h))   # exp(a)/cosh(kh), without overflow
    propagating = _gegenbauer_constant(p) * ive(2 * p + LAMBDA, a) * a ** (-LAMBDA) * ratio
    kn = np.asarray(kn, dtype=float)
    evanescent = _cos_integral(basis.p_max, kn * g) / np.cos(kn * h)[None, :]
    return 0.5 * g * np.column_stack([propagating, evanescent])


def inner_interior(cyl, n_max, basis):
    """F[p, n] = int_-h^-d u_p(t(z)) cos(lambda_n (z + h)) dz, lambda_n = n pi/g, n = 0 ... n_max;
    shape (p_max + 1, n_max + 1)."""
    return 0.5 * cyl.gap * _cos_integral(basis.p_max, np.arange(n_max + 1) * math.pi)


def _gauss_jacobi_exterior(cyl, k, kn, basis):
    """inner_exterior by the Gauss-Jacobi rule of the basis on the reflected interval, z + h = g abs(t). For the checks."""
    s = cyl.gap * np.abs(basis.nodes)
    h = cyl.depth
    modes = [np.cosh(k * s) / math.cosh(k * h)] + [np.cos(q * s) / math.cos(q * h) for q in kn]
    return 0.5 * cyl.gap * (basis.poly * basis.weights) @ np.array(modes).T


def _gauss_jacobi_interior(cyl, n_max, basis):
    """inner_interior by the Gauss-Jacobi rule of the basis. For the checks."""
    t = np.abs(basis.nodes)
    modes = np.cos(np.arange(n_max + 1)[:, None] * math.pi * t[None, :])
    return 0.5 * cyl.gap * (basis.poly * basis.weights) @ modes.T


def check_roots():
    """Residuals and brackets of wave_numbers; see the module docstring."""
    cyl = Cylinder()
    h = cyl.depth
    worst_k = worst_kn = 0.0
    bracket = True
    for omega in drift_convergence.WIDE_FREQUENCIES:
        nu = omega ** 2 / make_mesh.G
        k, kn = wave_numbers(omega, h, CHECK_N)
        worst_k = max(worst_k, abs(k * math.tanh(k * h) - nu) / nu)
        worst_kn = max(worst_kn, np.max(np.abs(kn * np.tan(kn * h) + nu)) / nu)
        n = np.arange(1, CHECK_N + 1)
        bracket &= bool(np.all(((n - 0.5) * math.pi / h < kn) & (kn < n * math.pi / h)))
    ok = worst_k < 1e-12 and worst_kn < 1e-12 and bracket
    return ok, "max |k tanh kh - nu|/nu %.1e, max |kn tan kn h + nu|/nu %.1e, brackets %s" % (worst_k, worst_kn, bracket)


def _quad_gap(cyl, p, mode):
    """int over -h < z < -d of u_p(t(z)) mode(z) dz by QUADPACK on the unreflected gap, the corner weight (-d - z)**(-1/3)
    handled by weight="alg"; (1 - t)**(-1/3) = g**(1/3) (-d - z)**(-1/3), (1 + t)**(-1/3) is smooth on the gap."""
    g, h, d = cyl.gap, cyl.depth, cyl.draft

    def f(z):
        t = (z + h) / g
        return eval_gegenbauer(2 * p, LAMBDA, t) * (1.0 + t) ** ALPHA * mode(z)

    # epsabs: entries with p large against a slowly varying mode are near zero by cancellation, and a smaller absolute
    # tolerance is below their rounding floor (QUADPACK warns). The matrices' largest entries are above 60.
    value, _ = quad(f, -h, -d, weight="alg", wvar=(0.0, ALPHA), epsabs=1e-12, epsrel=1e-12, limit=400)
    return g ** (-ALPHA) * value


def _relative(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


def check_inner_products():
    """inner_exterior, inner_interior and the Gauss-Jacobi sums against quad; see the module docstring."""
    cyl = Cylinder()
    h = cyl.depth
    basis = gegenbauer_basis(CHECK_P)
    worst = {"exterior": 0.0, "interior": 0.0, "exterior GJ": 0.0, "interior GJ": 0.0}
    for omega in drift_convergence.WIDE_FREQUENCIES:
        k, kn = wave_numbers(omega, h, CHECK_N)
        modes = [lambda z: math.cosh(k * (z + h)) / math.cosh(k * h)]
        modes += [lambda z, q=q: math.cos(q * (z + h)) / math.cos(q * h) for q in kn]
        ref = np.array([[_quad_gap(cyl, p, m) for m in modes] for p in range(CHECK_P + 1)])
        worst["exterior"] = max(worst["exterior"], _relative(inner_exterior(cyl, k, kn, basis), ref))
        worst["exterior GJ"] = max(worst["exterior GJ"], _relative(_gauss_jacobi_exterior(cyl, k, kn, basis), ref))
    lam = np.arange(CHECK_N + 1) * math.pi / cyl.gap
    modes = [lambda z, q=q: math.cos(q * (z + h)) for q in lam]
    ref = np.array([[_quad_gap(cyl, p, m) for m in modes] for p in range(CHECK_P + 1)])
    worst["interior"] = _relative(inner_interior(cyl, CHECK_N, basis), ref)
    worst["interior GJ"] = _relative(_gauss_jacobi_interior(cyl, CHECK_N, basis), ref)
    ok = all(v < 1e-10 for v in worst.values())
    return ok, ", ".join("%s %.1e" % item for item in worst.items()) + " (max |difference|/max |quad|)"


CHECKS = {"roots": check_roots, "inner_products": check_inner_products}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p_check = sub.add_parser("check")
    p_check.add_argument("names", nargs="*", metavar="NAME")
    args = parser.parse_args(argv)
    names = args.names or list(CHECKS)
    unknown = [n for n in names if n not in CHECKS]
    if unknown:
        parser.error("unknown check(s) %s; known: %s" % (", ".join(unknown), ", ".join(CHECKS)))
    failed = 0
    for name in names:
        try:
            ok, message = CHECKS[name]()
        except Exception as exc:  # a check that cannot run has failed
            ok, message = False, "%s: %s" % (type(exc).__name__, exc)
        failed += not ok
        print("%s %s: %s" % ("PASS" if ok else "FAIL", name, message))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
