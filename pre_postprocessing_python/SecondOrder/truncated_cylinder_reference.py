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

Matching (solve_diffraction, _solve_order). Each order m is matched on its own, phi = sum_m phi_m(r, z) cos(m theta). The
unknown is the radial velocity on the gap, dphi_m/dr(a, z) = sum_p c_p u_p(t), shared by both sides, so the velocity is
continuous by construction. The exterior coefficients come from projecting dphi_m/dr at r = a on Z_n over the whole depth
(the body velocity on the wall, zero for diffraction, and u on the gap), the interior ones from projecting it on
cos(lambda_n (z + h)) over the gap. Continuity of phi on the gap, projected on u_p, p = 0 ... P, is then a complex system of
size P + 1 in c_p (a Galerkin method):

    sum_q c_q [sum_n E_pn E_qn/(N_n kappa_n) - sum_n F_pn F_qn/(nu_n mu_n)] = -E_p0 inc_m 2i/(pi ka H_m'(ka)) + (radiation),

E = inner_exterior, F = inner_interior, N_n = int_-h^0 Z_n**2 dz, kappa_0 = k H_m'(ka)/H_m(ka), kappa_n = k_n K_m'/K_m(k_n a),
mu_0 = m/a, mu_n = lambda_n I_m'/I_m(lambda_n a), nu_0 = g, nu_n = g/2. The incident term is J_m - J_m' H_m/H_m' at ka, by the
Wronskian. Ratios use hankel1/h1vp, kve and ive only, so nothing overflows at small gap, large P, N or m. For m = 0 the
constant interior mode has no radial velocity: F_q0 = 0 for q >= 1 (orthogonality of C_2q^(1/6) to C_0), so the zero-flux
condition fixes c_0 (zero for diffraction) and the constant takes its place among the unknowns; the system stays P + 1.

Radiation (solve_radiation). Unit velocity in mode j, dphi_j/dn = n_j on the body with n out of the fluid, which is HAMS's
right-hand side (RBC_RIGHT in src/BodyIntgr.f90, the same normal as the -dphi_I/dn of DBC_RIGHT), n_5 = z n_x - x n_z.
Surge and pitch are order 1, heave order 0. On the wall, n_x = -cos(theta), so dphi_1/dr = cos(theta) and
dphi_5/dr = z cos(theta): they enter _solve_order as wall = int_-d^0 Z_n dz and int_-d^0 z Z_n dz (_wall_moments). On the
bottom, n_z = 1, so dphi_3/dz = 1 and dphi_5/dz = -r cos(theta), which no homogeneous interior mode carries (each has
dphi/dz = 0 at z = -d). The interior particular solutions, s = z + h,

    heave: phi_p = (s**2 - r**2/2)/(2 g),       pitch: phi_p = -r (s**2 - r**2/4)/(2 g),

are harmonic for their order, meet the bottom condition and dphi/dz = 0 on the bed. They enter the matching through their
potential on the gap projected on u_p (Gauss-Jacobi, exact for this even quadratic in t) and their radial velocity on the
gap projected on cos(lambda_n s), in closed form; the homogeneous interior series then carries u minus that velocity. For
heave the m = 0 flux condition gives c_0 = flux_gap[0]/F_00: the gap carries the outflow int dphi_p/dr dz = -a/2 that the
rising bottom displaces. The bottom integral of phi_p is closed form.

Added mass and damping, HAMS's normalisation. RFORCE in src/PotentWavForce.f90 sums, over the panels, with the mesh normal
and n_5 as above, AMAS(i, j) = Im(i rho int phi_j n_i dS) = rho Re int phi_j n_i dS and BDMP(i, j) =
-Re(i rho omega int phi_j n_i dS) = rho omega Im int phi_j n_i dS, i the force, j the motion. This is the force on the body,
i omega rho int phi_j n_i dS per unit velocity, written as i omega A - B in exp(-i omega t). It writes to unit 61,
AmssDamp.1 (opened in src/InputFiles.f90), the rows OUFR, i, j, AMAS(i, j)/RHO, BDMP(i, j)/(RHO W1): the WAMIT form
A/(rho L**k) and B/(rho omega L**k) with L = 1 (it never divides by REFL; the study's ControlFile has REFL = 1.0 anyway).
added_mass_damping returns Re and Im of int phi_j n_i dS, moments about the origin, the XR of the study.

Convergence. The exterior and interior series of dphi/dr at r = a are expansions of a velocity that is singular at the
corner and jumps there, so their pointwise value converges slowly in N, about N**(-2/3) a few metres from the corner; phi
converges about as N**(-1.6), and so do the excitation forces, the added mass and the damping. The truncation in N of the
sums in the matrix grows with P, so N must grow with P (N >= 100 P or so); at fixed N, larger P gives slightly worse
forces. N does not grow the linear system (size P + 1), so N in the thousands is cheap: a solve at N = 16000, P = 20 takes
about 0.5 s.

Incident wave and HAMS's excitation normalisation. HAMS's incident potential is VINP in src/PatclVelct.f90: for finite
depth, VINP = DUM cosh(k (z + h)) e^(i k x), DUM = -A omega/(k sinh kh), which is phi_I = -(g A/omega) Z_0(z) e^(i k x) by the
dispersion relation. Its sign with e^(i k x) for heading 0 along +x makes the time factor exp(-i omega t), and the incident
elevation, eta = (i omega/g) phi(z = 0), is eta_I = -i A e^(i k x). solve_diffraction uses this phi_I with A = 1. EFORCE in
src/PotentWavForce.f90 computes X = i omega rho int phi n_j dS with the normals of the mesh, which point out of the fluid
(make_mesh.check_normals), and n_5 = (z - z_r) n_x - (x - x_r) n_z about the rotation centre (CalRotNormals in
src/NormalProcess.f90), here the origin. It writes to unit 62, ExcForce.3 (opened in src/InputFiles.f90), the columns
MOD = abs(X)/(rho g A), NREL = -Im X/(rho g A), NIMG = -Re X/(rho g A). NREL + i NIMG = conj(X/(rho g eta_I(0))): the WAMIT
convention, exp(+i omega t) with the phase relative to the incident elevation at the origin. excitation returns
NREL + i NIMG, moments in m (HAMS divides moments by rho g A only). EFORCE's PHS column applies 180/pi to ATAN2D, which is
already in degrees, so it is not used. In a one-off comparison (not a check, and not reproducible from the repository:
the HAMS output is outside it) with the ExcForce.3 of the study's finest sharp-hull HAMS run, source method, the three
components agree in sign, phase and modulus to about 1e-3, the mesh error.

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
  check does.
- bottom_limit_x1: gap 1e-3 h (BOTTOM_GAP), CHECK_RESOLUTION, every WIDE_FREQUENCIES: X_1 against MacCamy-Fuchs,
  4 tanh(kh)/(k**2 conj(H_1'(ka))), to 1e-3 in modulus and as a complex number (which tests the conjugation). The gap
  itself changes X_1 by O(g/h): 1.2e-3 at 0.40 rad/s, linear in g from g = 1 m down to 1 mm, converged in N and P. So X_1 is
  extrapolated to g = 0 from g and 2g, 2 X(g) - X(2g), and that is compared; the message also prints the difference at
  g = 1e-3 h itself. The extrapolation assumes X - MF is linear in g, so the check also asserts
  abs(X(2g) - MF)/abs(X(g) - MF) in [1.8, 2.2] (measured 1.955) at each omega where abs(X(g) - MF) is above 1e-13 abs(MF);
  at 1.80 rad/s (kh = 33) the gap changes X_1 by about 1e-15, round-off, and the ratio is not tested there (the message
  says so). The extrapolated value is blind to errors of the interior and gap terms: they enter X_1 at O(g), and the
  extrapolation removes O(g); the linearity assertion and the raw difference at g = 1e-3 h are what see them. It covers
  the solve and the excitation at the gap where the interior modes are steepest, and the phase convention of excitation.
  It does not cover X_3 or X_5, nor the corner of a deep gap.
- matching: omega = 0.95, orders 0, 1, 2, 5, at 10, 50 and 90 % of the gap below the corner: the largest pointwise
  difference between the exterior and the interior series at r = a, of phi and of dphi/dr, each over the largest exterior
  value, falls monotonically for P = 10, 20, 40. N is tied to P, N = 200 P (MATCHING_MODES_PER_P): the velocity is matched by
  construction, so its pointwise mismatch is the N truncation of the two series, and at fixed N the mismatch of phi stops
  falling with P (see Convergence). It covers the matching and the corner basis. It does not measure the error at the corner.
- scattering_energy: fixed body, default cylinder, CHECK_RESOLUTION, every WIDE_FREQUENCIES: the scattered energy flux
  2 pi sum_m abs(a_m)**2 equals the optical-theorem term -sqrt(8 pi/k) Re(e^(3 i pi/4) a(0)), a(0) = sum_m a_m, from
  far_field_amplitude, to 1e-6. Energy is conserved order by order, so it tests each order's propagating coefficient and the
  far-field normalisation, not the convergence in M, N or P.
- reciprocity: default cylinder, every WIDE_FREQUENCIES: A15 = A51 and B15 = B51 to 1e-8, and X_j from diffraction equal
  to haskind_excitation of the mode j radiation solution to 1e-6, j = 1, 3, 5. Values at N = CHECK_RESOLUTION's 8000 and
  2 N = 16000, P = 20, extrapolated in N with the exponent 5/3 (Richardson); the message also prints the residuals at
  N = 16000 itself (about 5e-7). Surge's identity holds to round-off at any truncation (symmetric Galerkin system); heave
  and pitch converge as N**(-5/3) through the particular solution. check_reciprocity's docstring says what each part tests.
- radiation_energy: same resolution and extrapolation: B11, B33 and B55 from the body pressure equal
  radiated_energy_damping, the energy flux at infinity, to 1e-6. It tests the flux formula and the sign and normalisation of
  B, and for heave and pitch the particular solutions. Neither radiation check bounds the error of A, B or X themselves.
"""
import argparse
import functools
import math
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer, gammaln, h1vp, hankel1, ive, jv, jvp, kve, roots_jacobi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drift_convergence  # noqa: E402
import make_mesh  # noqa: E402

LAMBDA = 1.0 / 6.0     # Gegenbauer parameter of the corner basis
ALPHA = LAMBDA - 0.5   # Jacobi parameter alpha = beta = -1/3 of its weight
GAUSS_JACOBI_NODES = 96   # default size of the Gauss-Jacobi rule: exact for polynomials of degree < 192

# Resolution of the root and inner product checks
CHECK_P = 10
CHECK_N = 20
# Resolution of the solution checks: N exterior and interior modes, corner basis size P, orders 0 ... M. N = 400 P, above the
# N >= 100 P rule (Convergence); the radiation checks extrapolate from N and 2 N (RICHARDSON_EXPONENT).
CHECK_RESOLUTION = dict(n_modes=8000, p_max=20, orders=15)
RICHARDSON_EXPONENT = 5.0 / 3.0   # the radiation identities' residual falls as N**(-5/3); see check_reciprocity


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


def _cosh_ratio(k, s, h):
    """cosh(k s)/cosh(k h) and sinh(k s)/cosh(k h) for 0 <= s <= h, without overflow."""
    s = np.asarray(s, dtype=float)
    den = 1.0 + math.exp(-2.0 * k * h)
    plus, minus = np.exp(k * (s - h)), np.exp(-k * (s + h))
    return (plus + minus) / den, (plus - minus) / den


def _mode_norms(k, kn, h):
    """N_n = int_-h^0 Z_n**2 dz, n = 0 ... len(kn)."""
    e = math.exp(-2.0 * k * h)
    n0 = h * 2.0 * e / (1.0 + e) ** 2 + math.tanh(k * h) / (2.0 * k)    # (h/2 + sinh(2kh)/(4k))/cosh(kh)**2
    nn = (0.5 * h + np.sin(2.0 * kn * h) / (4.0 * kn)) / np.cos(kn * h) ** 2
    return np.concatenate([[n0], nn])


def _wall_moments(cyl, k, kn):
    """(I0, I1): I0[n] = int_-d^0 Z_n dz and I1[n] = int_-d^0 z Z_n dz on the wall, n = 0 ... len(kn), in closed form
    (s = z + h from g to h)."""
    g, h, d = cyl.gap, cyl.depth, cyl.draft
    ch_g, sh_g = _cosh_ratio(k, g, h)
    i0 = [(math.tanh(k * h) - sh_g) / k]
    i1 = [d * sh_g / k - (1.0 - ch_g) / k ** 2]
    c = np.cos(kn * h)
    i0 = np.concatenate([i0, (np.sin(kn * h) - np.sin(kn * g)) / kn / c])
    i1 = np.concatenate([i1, (d * np.sin(kn * g) / kn + (np.cos(kn * h) - np.cos(kn * g)) / kn ** 2) / c])
    return i0, i1


@dataclass
class Order:
    """The solution of one azimuthal order m; phi_m(r, z) multiplies cos(m theta).

    Exterior, r >= a: phi_m = Z_0(z) (incident J_m(k r) + scattered H_m(k r)/H_m(k a))
                              + sum_n evanescent[n-1] Z_n(z) K_m(k_n r)/K_m(k_n a).
    Interior, r <= a, z <= -d: phi_m = interior[0] (r/a)**m
                              + sum_n interior[n] cos(lambda_n (z + h)) I_m(lambda_n r)/I_m(lambda_n a) + particular.
    Gap, r = a: dphi_m/dr = sum_p gap[p] u_p(t). The coefficients are amplitudes at r = a; the ratios are evaluated with
    hankel1, kve and ive and never overflow."""
    m: int
    incident: complex
    scattered: complex
    evanescent: np.ndarray
    interior: np.ndarray
    gap: np.ndarray
    particular: object = None


@dataclass
class Solution:
    """A first-order potential phi(r, theta, z) = sum_m phi_m(r, z) cos(m theta), in m**2/s per unit wave amplitude, time
    factor exp(-i omega t). For diffraction it is the total potential, incident plus scattered, with HAMS's incident
    potential (see the module docstring). A radiation solution (solve_radiation) uses the same type with incident = 0, one
    order, per unit body velocity, and for heave and pitch a particular solution in the interior.

    particular, per order, is the hook for the interior particular solution of heave and pitch (_Particular). It is an
    object with
        fields(r, z) -> (phi, dphi/dr, dphi/dz) of the particular solution, r <= a, -h <= z <= -d, arrays,
        potential_gap[p] = int over the gap of phi_p(a, z) u_p(t(z)) dz, p = 0 ... P,
        flux_gap[n] = int over the gap of dphi_p/dr(a, z) cos(lambda_n (z + h)) dz, n = 0 ... N,
        bottom = int_0^a phi_p(r, -d) r**(m+1) dr, the weight of the order's bottom integral in _body_integrals,
    and the radial body velocity on the wall enters _solve_order as wall[n] = int_-d^0 v_r(z) Z_n(z) dz."""
    cyl: Cylinder
    omega: float
    k: float
    kn: np.ndarray
    p_max: int
    orders: dict

    @property
    def lam(self):
        return np.arange(len(self.kn) + 1) * math.pi / self.cyl.gap

    def _exterior(self, m, r, z):
        """(phi_m, dphi_m/dr, dphi_m/dz) of the exterior series at r >= a, -h <= z <= 0, arrays of one shape."""
        o, cyl, k, kn = self.orders[m], self.cyl, self.k, self.kn
        a, h = cyl.radius, cyl.depth
        r, z = np.asarray(r, dtype=float), np.asarray(z, dtype=float)
        s = z + h
        z0, dz0 = _cosh_ratio(k, s, h)
        dz0 = k * dz0
        ha = hankel1(m, k * a)
        radial = o.incident * jv(m, k * r) + o.scattered * hankel1(m, k * r) / ha
        dradial = k * (o.incident * jvp(m, k * r) + o.scattered * h1vp(m, k * r) / ha)
        phi, dr, dz = z0 * radial, z0 * dradial, dz0 * radial
        if len(kn):
            q = kn.reshape((-1,) + (1,) * r.ndim)
            x = q * r
            decay = np.exp(-q * (r - a)) / kve(m, q * a)
            kr = kve(m, x) * decay
            dkr = -0.5 * q * (kve(m - 1, x) + kve(m + 1, x)) * decay
            zn = np.cos(q * s) / np.cos(q * h)
            dzn = -q * np.sin(q * s) / np.cos(q * h)
            b = o.evanescent.reshape(q.shape)
            phi = phi + np.sum(b * zn * kr, axis=0)
            dr = dr + np.sum(b * zn * dkr, axis=0)
            dz = dz + np.sum(b * dzn * kr, axis=0)
        return phi, dr, dz

    def _interior(self, m, r, z):
        """(phi_m, dphi_m/dr, dphi_m/dz) of the interior series at r <= a, -h <= z <= -d, arrays of one shape."""
        o, cyl = self.orders[m], self.cyl
        a, h = cyl.radius, cyl.depth
        r, z = np.asarray(r, dtype=float), np.asarray(z, dtype=float)
        c = o.interior
        phi = c[0] * (r / a) ** m + 0j * r
        dr = c[0] * m * r ** max(m - 1, 0) / a ** m if m else 0j * r
        dz = 0j * r
        lam = self.lam[1:].reshape((-1,) + (1,) * r.ndim)
        if lam.size:
            x = lam * r
            grow = np.exp(lam * (r - a)) / ive(m, lam * a)
            ir = ive(m, x) * grow
            dir_ = 0.5 * lam * (ive(m - 1, x) + ive(m + 1, x)) * grow
            s = z + h
            cn, sn = np.cos(lam * s), -lam * np.sin(lam * s)
            b = c[1:].reshape(lam.shape)
            phi = phi + np.sum(b * cn * ir, axis=0)
            dr = dr + np.sum(b * cn * dir_, axis=0)
            dz = dz + np.sum(b * sn * ir, axis=0)
        if o.particular is not None:
            pp, pr, pz = o.particular.fields(r, z)
            phi, dr, dz = phi + pp, dr + pr, dz + pz
        return phi, dr, dz

    def _fields(self, r, theta, z):
        """phi and its Cartesian gradient at broadcast points, nan inside the body."""
        r, theta, z = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (r, theta, z)))
        a, d, h = self.cyl.radius, self.cyl.draft, self.cyl.depth
        outside = (r >= a) & (z <= 0.0) & (z >= -h)
        inside = (r < a) & (z <= -d) & (z >= -h)
        rs = np.maximum(r, 1e-9 * a)   # the gradient's (1/r) d/dtheta term at the axis, as its limit to O(1e-9 a)
        phi, gr, gt, gz = (np.full(r.shape, np.nan + 0j) for _ in range(4))
        for mask, series in ((outside, self._exterior), (inside, self._interior)):
            if not mask.any():
                continue
            rr, tt, zz = rs[mask], theta[mask], z[mask]
            p = np.zeros(rr.shape, complex)
            fr, ft, fz = p.copy(), p.copy(), p.copy()
            for m in self.orders:
                pm, prm, pzm = series(m, rr, zz)
                cm, sm = np.cos(m * tt), np.sin(m * tt)
                p += pm * cm
                fr += prm * cm
                ft += -m * pm * sm / rr
                fz += pzm * cm
            phi[mask], gr[mask], gt[mask], gz[mask] = p, fr, ft, fz
        ct, st = np.cos(theta), np.sin(theta)
        return phi, np.array([ct * gr - st * gt, st * gr + ct * gt, gz])

    def potential(self, r, theta, z):
        """phi at (r, theta, z) anywhere in the fluid (r >= a, or r < a below the bottom), nan inside the body; broadcasts."""
        return self._fields(r, theta, z)[0]

    def gradient(self, r, theta, z):
        """Cartesian gradient (dphi/dx, dphi/dy, dphi/dz), x = r cos(theta), y = r sin(theta), z up; shape (3,) + the
        broadcast shape of the arguments, nan inside the body."""
        return self._fields(r, theta, z)[1]


def _solve_order(cyl, m, k, kn, basis, incident=0.0, wall=None, particular=None):
    """One order of the matching; see the module docstring. A complex system of size P + 1."""
    a, g = cyl.radius, cyl.gap
    n_modes = len(kn)
    big_e = inner_exterior(cyl, k, kn, basis)                 # (P+1, N+1)
    big_f = inner_interior(cyl, n_modes, basis)              # (P+1, N+1)
    norms = _mode_norms(k, kn, cyl.depth)
    ka = k * a
    kappa = np.empty(n_modes + 1, complex)                   # Z_n' radial log derivatives at r = a
    kappa[0] = k * h1vp(m, ka) / hankel1(m, ka)
    x = kn * a
    kappa[1:] = -0.5 * kn * (kve(m - 1, x) + kve(m + 1, x)) / kve(m, x)
    lam = np.arange(n_modes + 1) * math.pi / g
    mu = np.empty(n_modes + 1)
    mu[0] = m / a
    y = lam[1:] * a
    mu[1:] = 0.5 * lam[1:] * (ive(m - 1, y) + ive(m + 1, y)) / ive(m, y)
    nu_int = np.full(n_modes + 1, 0.5 * g)
    nu_int[0] = g
    wall = np.zeros(n_modes + 1) if wall is None else np.asarray(wall)
    flux_p = np.zeros(n_modes + 1) if particular is None else np.asarray(particular.flux_gap)
    pot_p = np.zeros(basis.p_max + 1) if particular is None else np.asarray(particular.potential_gap)

    ext = big_e / (norms * kappa)                           # E[p, n]/(N_n kappa_n)
    live = mu != 0.0                                        # m = 0, n = 0: no radial velocity, constant unknown
    intr = np.where(live, big_f / np.where(live, nu_int * mu, 1.0), 0.0)
    matrix = ext @ big_e.T - intr @ big_f.T
    # incident: E_p0 inc (J_m - J_m' H_m/H_m') = E_p0 inc 2i/(pi ka H_m'(ka)), by the Wronskian
    incident_at_wall = incident * 2j / (math.pi * ka * h1vp(m, ka))
    rhs = pot_p - big_e[:, 0] * incident_at_wall - ext @ wall - intr @ flux_p
    if m == 0:
        # F[q, 0] = 0 for q >= 1, so the flux condition gives gap[0]; the constant interior[0] takes its place
        c0 = flux_p[0] / big_f[0, 0]
        rhs = rhs - matrix[:, 0] * c0
        matrix = matrix.copy()
        matrix[:, 0] = -big_f[:, 0]
        sol = np.linalg.solve(matrix, rhs)
        gap, const = np.concatenate([[c0], sol[1:]]), sol[0]
    else:
        gap = np.linalg.solve(matrix, rhs)
        const = None
    exterior = (gap @ big_e + wall) / (norms * kappa)       # amplitudes at r = a of the scattered parts
    scattered = exterior[0] - incident * k * jvp(m, ka) / kappa[0]
    interior = np.where(live, (gap @ big_f - flux_p) / np.where(live, nu_int * mu, 1.0), 0.0).astype(complex)
    if const is not None:
        interior[0] = const
    return Order(m, complex(incident), complex(scattered), exterior[1:], interior, gap, particular)


def _order_list(orders):
    return list(range(orders + 1)) if isinstance(orders, (int, np.integer)) else [int(m) for m in orders]


def solve_diffraction(cyl, omega, orders, n_modes, p_max):
    """Fixed body in HAMS's incident wave of unit amplitude along +x. orders: an int M (orders 0 ... M) or the orders to
    solve. n_modes: evanescent exterior modes and interior modes N. p_max: corner basis size P."""
    k, kn = wave_numbers(omega, cyl.depth, n_modes)
    basis = gegenbauer_basis(p_max)
    amplitude = -make_mesh.G / omega                       # phi_I = -(g A/omega) Z_0 e^(i k x), A = 1
    solved = {}
    for m in _order_list(orders):
        inc = amplitude * (1.0 if m == 0 else 2.0) * 1j ** m   # e^(ikx) = sum eps_m i^m J_m(kr) cos(m theta)
        solved[m] = _solve_order(cyl, m, k, kn, basis, incident=inc)
    return Solution(cyl, omega, k, kn, p_max, solved)


class _Particular:
    """Interior particular solution of heave (mode 3, order 0) or pitch (mode 5, order 1), unit velocity, s = z + h:

        heave: phi_p = (s**2 - r**2/2)/(2 g),        dphi_p/dz = s/g,       1 on z = -d, 0 on the bed,
        pitch: phi_p = -r (s**2 - r**2/4)/(2 g),     dphi_p/dz = -r s/g,   -r on z = -d, 0 on the bed,

    both harmonic for their order. The hook of Solution: fields, potential_gap (Gauss-Jacobi, exact, phi_p(a, z) is an even
    quadratic in t), flux_gap (closed form), and bottom, the bottom integral of phi_p(r, -d) r**(m+1) dr over 0 < r < a, in
    the weight _body_integrals uses for the order."""

    def __init__(self, cyl, mode, basis, n_modes):
        if mode not in (3, 5):
            raise ValueError("particular solution: mode 3 or 5, got %s" % mode)
        self.cyl, self.mode = cyl, mode
        a, g, h = cyl.radius, cyl.gap, cyl.depth
        phi_a = self.fields(np.full(basis.nodes.shape, a), g * np.abs(basis.nodes) - h)[0]
        self.potential_gap = 0.5 * g * (basis.poly * basis.weights) @ phi_a
        lam = np.arange(1, n_modes + 1) * math.pi / g
        one = np.concatenate([[g], np.zeros(n_modes)])                                # int_0^g cos(lambda_n s) ds
        s2 = np.concatenate([[g ** 3 / 3.0], 2.0 * g * (-1.0) ** np.arange(1, n_modes + 1) / lam ** 2])   # s**2 cos
        if mode == 3:
            self.flux_gap = -a / (2.0 * g) * one
            self.bottom = (g ** 2 * a ** 2 / 2.0 - a ** 4 / 8.0) / (2.0 * g)
        else:
            self.flux_gap = -(s2 - 0.75 * a ** 2 * one) / (2.0 * g)
            self.bottom = -(g ** 2 * a ** 4 / 4.0 - a ** 6 / 24.0) / (2.0 * g)

    def fields(self, r, z):
        r, z = np.asarray(r, dtype=float), np.asarray(z, dtype=float)
        g, s = self.cyl.gap, z + self.cyl.depth
        if self.mode == 3:
            return (s ** 2 - r ** 2 / 2.0) / (2.0 * g), -r / (2.0 * g), s / g
        return -r * (s ** 2 - r ** 2 / 4.0) / (2.0 * g), -(s ** 2 - 0.75 * r ** 2) / (2.0 * g), -r * s / g


RADIATION_ORDER = {1: 1, 3: 0, 5: 1}   # azimuthal order of each radiation mode


def solve_radiation(cyl, omega, mode, n_modes, p_max):
    """Radiation potential of mode 1 (surge), 3 (heave) or 5 (pitch about the origin), unit velocity: dphi/dn = n_j on the
    body, n out of the fluid, as HAMS's right-hand side (RBC_RIGHT in src/BodyIntgr.f90). Wall: dphi_1/dr = cos(theta),
    dphi_5/dr = z cos(theta); bottom: dphi_3/dz = 1, dphi_5/dz = -r cos(theta), by the particular solutions (_Particular).
    One order, solved by _solve_order as diffraction is."""
    if mode not in RADIATION_ORDER:
        raise ValueError("solve_radiation: mode 1, 3 or 5, got %s" % mode)
    k, kn = wave_numbers(omega, cyl.depth, n_modes)
    basis = gegenbauer_basis(p_max)
    m = RADIATION_ORDER[mode]
    i0, i1 = _wall_moments(cyl, k, kn)
    wall = {1: i0, 3: None, 5: i1}[mode]
    particular = None if mode == 1 else _Particular(cyl, mode, basis, n_modes)
    return Solution(cyl, omega, k, kn, p_max, {m: _solve_order(cyl, m, k, kn, basis, wall=wall, particular=particular)})


def _body_integrals(cyl, sol):
    """Integrals of phi n_j over the wetted surface, n out of the fluid (into the body, as in HAMS), j = 1, 3, 5, moments
    about the origin."""
    a, d = cyl.radius, cyl.draft
    i0, i1 = _wall_moments(cyl, sol.k, sol.kn)
    lam = sol.lam[1:]
    sign = (-1.0) ** np.arange(1, len(lam) + 1)            # cos(lambda_n g)
    out = {}
    if 1 in sol.orders:
        o = sol.orders[1]
        at_wall = np.concatenate([[o.incident * jv(1, sol.k * a) + o.scattered], o.evanescent])
        wall_1 = -math.pi * a * (at_wall @ i0)                 # n_x = -cos(theta) on the wall
        wall_5 = -math.pi * a * (at_wall @ i1)                 # z n_x
        y = lam * a
        bottom = o.interior[0] * a ** 3 / 4.0 + np.sum(o.interior[1:] * sign * a ** 2 * ive(2, y) / (lam * ive(1, y)))
        if o.particular is not None:
            bottom = bottom + o.particular.bottom
        out[1], out[5] = wall_1, wall_5 - math.pi * bottom   # -x n_z, n_z = 1 on the bottom
    if 0 in sol.orders:
        o = sol.orders[0]
        y = lam * a
        bottom = o.interior[0] * a ** 2 / 2.0 + np.sum(o.interior[1:] * sign * a * ive(1, y) / (lam * ive(0, y)))
        if o.particular is not None:
            bottom = bottom + o.particular.bottom
        out[3] = 2.0 * math.pi * bottom
    return out


def excitation(cyl, sol):
    """X_j, j = 1, 3, 5, of a diffraction solution, in HAMS's ExcForce.3 normalisation (module docstring):
    conj(X/(rho g eta_I(0))), X = i omega rho int phi n_j dS, eta_I(0) = -i A. Over rho g A, moments in m."""
    if 0 not in sol.orders or 1 not in sol.orders:
        raise ValueError("excitation needs orders 0 and 1")
    if cyl != sol.cyl:
        raise ValueError("excitation: the cylinder is not the one the solution was solved for")
    return {j: _hams_excitation(1j * sol.omega * make_mesh.RHO * value) for j, value in _body_integrals(cyl, sol).items()}


def _hams_excitation(force):
    """A complex force or moment per unit wave amplitude, time factor exp(-i omega t), in ExcForce.3's normalisation."""
    return complex(np.conj(force / (make_mesh.RHO * make_mesh.G * -1j)))


def added_mass_damping(cyl, omega, rads):
    """(A, B), 6 x 6 real, in HAMS's AmssDamp.1 normalisation (module docstring): A[i-1, j-1] = A_ij/rho,
    B[i-1, j-1] = B_ij/(rho omega), force i, motion j, rotations with L = 1, moments about the origin. rads: {j: Solution}
    of solve_radiation for modes j among 1, 3, 5, all at this cylinder and omega. Rows and columns 2, 4, 6 stay zero, and so
    do the couplings of heave with surge and pitch (other azimuthal order)."""
    a, b = np.zeros((6, 6)), np.zeros((6, 6))
    for j, sol in rads.items():
        if j not in RADIATION_ORDER or set(sol.orders) != {RADIATION_ORDER[j]}:
            raise ValueError("added_mass_damping: rads[%s] is not the radiation solution of mode %s" % (j, j))
        if sol.cyl != cyl or sol.omega != omega:
            raise ValueError("added_mass_damping: rads[%s] was solved for another cylinder or frequency" % j)
        for i, value in _body_integrals(cyl, sol).items():
            a[i - 1, j - 1], b[i - 1, j - 1] = value.real, value.imag
    return a, b


def _propagating_norm(sol):
    """N_0 = int_-h^0 Z_0**2 dz and the azimuthal weight int_0^2pi cos(m theta)**2 dtheta of the single order of sol."""
    (m,) = sol.orders
    return m, _mode_norms(sol.k, np.empty(0), sol.cyl.depth)[0], (2.0 * math.pi if m == 0 else math.pi)


def haskind_excitation(rad):
    """X_j from the radiation solution of mode j by the Haskind relation, in ExcForce.3's normalisation. With Green's second
    identity on a vertical cylinder r = R in the far field, where only the propagating terms survive,
    X_j = i omega rho int_body (phi_I n_j - phi_j dphi_I/dn) dS = -(2 rho g/pi) eps_m i**m c_m N_0 S_m/H_m(ka), with
    phi_I = -(g/omega) Z_0 sum eps_m i**m J_m(kr) cos(m theta), c_m = int cos(m theta)**2 dtheta and S_m the propagating
    amplitude of phi_j at r = a."""
    m, n0, weight = _propagating_norm(rad)
    eps = 1.0 if m == 0 else 2.0
    s = rad.orders[m].scattered
    force = -(2.0 * make_mesh.RHO * make_mesh.G / math.pi) * eps * 1j ** m * weight * n0 * s
    force /= hankel1(m, rad.k * rad.cyl.radius)
    return _hams_excitation(force)


def radiated_energy_damping(rad):
    """B_jj/(rho omega) from the energy flux at infinity of the radiation solution of mode j: the mean power through r = R is
    -(rho omega/2) int Im(phi conj(dphi/dr)) dS = (rho omega/pi) c_m N_0 abs(S_m/H_m(ka))**2 by the Wronskian, and B_jj is
    twice it at unit velocity."""
    m, n0, weight = _propagating_norm(rad)
    return 2.0 * weight * n0 * abs(rad.orders[m].scattered / hankel1(m, rad.k * rad.cyl.radius)) ** 2 / math.pi


def far_field_amplitude(sol):
    """a_m, m = 0 ... M, so that eta_scattered ~ sum_(m=-M)^M a_|m| e^(i m theta) e^(i k r)/sqrt(r) as k r -> infinity, with
    eta = (i omega/g) phi(z = 0) (time factor exp(-i omega t)). a_-m = a_m by symmetry about the x axis. For a diffraction
    solution it is the scattered wave only. Needs the orders 0 ... M."""
    orders = sorted(sol.orders)
    if orders != list(range(len(orders))):
        raise ValueError("far_field_amplitude needs the orders 0 ... M, got %s" % orders)
    k, a = sol.k, sol.cyl.radius
    out = np.empty(len(orders), complex)
    for m in orders:
        c = (1j * sol.omega / make_mesh.G) * sol.orders[m].scattered / hankel1(m, k * a) * math.sqrt(2.0 / (math.pi * k))
        c *= np.exp(-1j * (m * math.pi / 2.0 + math.pi / 4.0))
        out[m] = c if m == 0 else 0.5 * c
    return out


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


def maccamy_fuchs_x1(omega, radius, depth):
    """Surge excitation of the bottom-mounted cylinder (MacCamy and Fuchs) in HAMS's normalisation, over rho g A:
    4 tanh(kh)/(k**2 conj(H_1'(ka))). Its modulus is 4 tanh(kh)/(k**2 abs(H_1'(ka)))."""
    k = drift_convergence.wave_number(omega, depth)
    return 4.0 * math.tanh(k * depth) / (k ** 2 * np.conj(h1vp(1, k * radius)))


BOTTOM_GAP = 1e-3        # gap of the bottom-mounted limit, over h
LINEARITY_RATIO = (1.8, 2.2)   # bounds of |X(2g) - MF|/|X(g) - MF|
LINEARITY_FLOOR = 1e-13  # |X(g) - MF|/|MF| below which the gap's effect is round-off and the ratio is not tested


def check_bottom_limit_x1():
    """Gap 1e-3 h: X_1 against MacCamy-Fuchs, after extrapolation in the gap, and the linearity of X_1 - MF in the gap;
    see the module docstring."""
    h = Cylinder().depth
    res = CHECK_RESOLUTION
    worst_raw = worst_mod = worst_complex = 0.0
    ratios, skipped = [], []
    for omega in drift_convergence.WIDE_FREQUENCIES:
        x = []
        for gap in (BOTTOM_GAP * h, 2.0 * BOTTOM_GAP * h):
            cyl = Cylinder(draft=h - gap)
            x.append(excitation(cyl, solve_diffraction(cyl, omega, (0, 1), res["n_modes"], res["p_max"]))[1])
        limit = 2.0 * x[0] - x[1]                            # removes the O(g) term of the gap
        ref = maccamy_fuchs_x1(omega, Cylinder().radius, h)
        worst_raw = max(worst_raw, abs(abs(x[0]) - abs(ref)) / abs(ref))
        worst_mod = max(worst_mod, abs(abs(limit) - abs(ref)) / abs(ref))
        worst_complex = max(worst_complex, abs(limit - ref) / abs(ref))
        if abs(x[0] - ref) > LINEARITY_FLOOR * abs(ref):
            ratios.append(abs(x[1] - ref) / abs(x[0] - ref))
        else:
            skipped.append(omega)
    linear = all(LINEARITY_RATIO[0] <= q <= LINEARITY_RATIO[1] for q in ratios)
    ok = worst_mod < 1e-3 and worst_complex < 1e-3 and linear
    return ok, ("extrapolated to g = 0: max ||X1| - |MF||/|MF| %.1e, max |X1 - MF|/|MF| %.1e; at g = 1e-3 h itself "
                "max ||X1| - |MF||/|MF| %.1e; |X(2g) - MF|/|X(g) - MF| %.3f ... %.3f%s"
                % (worst_mod, worst_complex, worst_raw, min(ratios), max(ratios),
                   ", not tested at %s rad/s (below the round-off floor)" % ", ".join("%.2f" % w for w in skipped)
                   if skipped else ""))


MATCHING_ORDERS = (0, 1, 2, 5)
MATCHING_P = (10, 20, 40)
MATCHING_FRACTIONS = (0.1, 0.5, 0.9)   # depth below the corner, as a fraction of the gap
MATCHING_MODES_PER_P = 200             # N = 200 P in the matching check


def check_matching():
    """Pointwise mismatch of the exterior and interior series at r = a on the gap; see the module docstring."""
    cyl = Cylinder()
    z = -cyl.draft - np.array(MATCHING_FRACTIONS) * cyl.gap
    r = np.full_like(z, cyl.radius)
    errors = []
    for p_max in MATCHING_P:
        sol = solve_diffraction(cyl, 0.95, MATCHING_ORDERS, MATCHING_MODES_PER_P * p_max, p_max)
        d_phi = d_dr = scale_phi = scale_dr = 0.0
        for m in MATCHING_ORDERS:
            pe, re_, _ = sol._exterior(m, r, z)
            pi, ri, _ = sol._interior(m, r, z)
            d_phi, d_dr = max(d_phi, np.max(np.abs(pe - pi))), max(d_dr, np.max(np.abs(re_ - ri)))
            scale_phi, scale_dr = max(scale_phi, np.max(np.abs(pe))), max(scale_dr, np.max(np.abs(re_)))
        errors.append((d_phi / scale_phi, d_dr / scale_dr))
    falls = all(errors[i + 1][j] < errors[i][j] for i in range(len(errors) - 1) for j in range(2))
    text = ", ".join("P = %d, N = %d: phi %.1e, dphi/dr %.1e" % (p, MATCHING_MODES_PER_P * p, e[0], e[1])
                     for p, e in zip(MATCHING_P, errors))
    return falls, text + " (max |difference| over max |exterior|)"


def check_scattering_energy():
    """Fixed body: scattered energy flux against the optical theorem; see the module docstring."""
    cyl = Cylinder()
    res = CHECK_RESOLUTION
    worst = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        sol = solve_diffraction(cyl, omega, res["orders"], res["n_modes"], res["p_max"])
        a = far_field_amplitude(sol)
        full = np.concatenate([a[:0:-1], a])                     # a_m, m = -M ... M
        scattered = 2.0 * math.pi * np.sum(np.abs(full) ** 2)
        optical = -math.sqrt(8.0 * math.pi / sol.k) * np.real(np.exp(0.75j * math.pi) * np.sum(full))
        worst = max(worst, abs(scattered - optical) / abs(optical))
    return worst < 1e-6, "max |scattered - optical|/|optical| %.1e" % worst


@functools.lru_cache(maxsize=None)
def _first_order(omega, n_modes, p_max):
    """Default cylinder: A, B, X from diffraction, X by Haskind and B_jj from the energy flux, j = 1, 3, 5 (arrays indexed
    j = 1, 3, 5 in that order). Cached, so the two radiation checks solve each case once."""
    cyl = Cylinder()
    rads = {j: solve_radiation(cyl, omega, j, n_modes, p_max) for j in (1, 3, 5)}
    a, b = added_mass_damping(cyl, omega, rads)
    x = excitation(cyl, solve_diffraction(cyl, omega, (0, 1), n_modes, p_max))
    return dict(a=a, b=b, x=np.array([x[j] for j in (1, 3, 5)]),
                haskind=np.array([haskind_excitation(rads[j]) for j in (1, 3, 5)]),
                flux=np.array([radiated_energy_damping(rads[j]) for j in (1, 3, 5)]))


def _first_order_extrapolated(omega):
    """_first_order at N and 2 N of CHECK_RESOLUTION, extrapolated in N: (2**p Q(2N) - Q(N))/(2**p - 1), p =
    RICHARDSON_EXPONENT. Also returns the values at 2 N."""
    n, p_max = CHECK_RESOLUTION["n_modes"], CHECK_RESOLUTION["p_max"]
    lo, hi = _first_order(omega, n, p_max), _first_order(omega, 2 * n, p_max)
    w = 2.0 ** RICHARDSON_EXPONENT
    return {key: (w * hi[key] - lo[key]) / (w - 1.0) for key in lo}, hi


def _reciprocity_errors(q):
    a, b = q["a"], q["b"]
    sym = max(abs(a[0, 4] - a[4, 0]) / abs(a[0, 4]), abs(b[0, 4] - b[4, 0]) / abs(b[0, 4]))
    return sym, np.max(np.abs(q["haskind"] - q["x"]) / np.abs(q["x"]))


def check_reciprocity():
    """A15 = A51 and B15 = B51 to 1e-8, and X_j from diffraction equal to the Haskind value from the radiation solution of
    mode j to 1e-6, j = 1, 3, 5, at every WIDE_FREQUENCIES, after extrapolation in N (CHECK_RESOLUTION, N and 2 N).

    What holds by construction. The matching is a symmetric Galerkin system (the matrix of _solve_order is complex
    symmetric, the exterior and interior traces are orthogonal projections, the gap velocity is shared), so where no
    particular solution enters, Green's identities hold for the truncated series exactly: surge's Haskind relation holds
    to round-off at any N and P (1e-16 at N = 50, P = 5). There, the check tests the Haskind formula, the far-field
    normalisation and the sign and conjugation of excitation and added_mass_damping, not convergence. Pitch and heave have
    the interior particular solution, which is not in the span of the N truncated interior modes; the residual of the
    identities is the tail of its cosine projection (coefficients n**-2) against that of the corner-singular gap velocity
    (n**-2/3), so it falls as N**(-5/3): measured exponent 1.68 to 1.73 for N = 1000 ... 8000, 1.3e-7 for A15 = A51 still at
    N = 32000. One Richardson step with p = 5/3 removes that term; what remains falls about as N**-2.9 (2e-8 at N = 4000,
    8000). So A15 = A51, B15 = B51 and Haskind for modes 3 and 5 test the particular solutions, their gap projections and
    bottom integrals, and the pitch wall term. They do not bound the error of A, B or X themselves, which converge as about
    N**-1.6 (Convergence)."""
    worst_sym = worst_haskind = raw_sym = raw_haskind = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        q, hi = _first_order_extrapolated(omega)
        sym, haskind = _reciprocity_errors(q)
        worst_sym, worst_haskind = max(worst_sym, sym), max(worst_haskind, haskind)
        sym, haskind = _reciprocity_errors(hi)
        raw_sym, raw_haskind = max(raw_sym, sym), max(raw_haskind, haskind)
    ok = worst_sym < 1e-8 and worst_haskind < 1e-6
    return ok, ("extrapolated in N: max |A15 - A51|/|A15|, |B15 - B51|/|B15| %.1e, max |X Haskind - X diffraction|/|X| "
                "%.1e, modes 1, 3, 5; at N = %d itself %.1e and %.1e"
                % (worst_sym, worst_haskind, 2 * CHECK_RESOLUTION["n_modes"], raw_sym, raw_haskind))


def check_radiation_energy():
    """B11, B33 and B55 from body pressure (added_mass_damping) equal the energy flux at infinity
    (radiated_energy_damping) to 1e-6, at every WIDE_FREQUENCIES, after extrapolation in N as in check_reciprocity. Surge
    holds to round-off at any truncation (Green's identity for the symmetric Galerkin solution, as check_reciprocity says),
    so for it the check tests the flux formula and the normalisation and sign of B. Heave and pitch have the residual
    N**(-5/3) of the particular solution and test it. It does not test the convergence of B."""
    def error(v):
        return np.max(np.abs(np.diag(v["b"])[[0, 2, 4]] - v["flux"]) / np.abs(v["flux"]))

    worst = raw = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        q, hi = _first_order_extrapolated(omega)
        worst, raw = max(worst, error(q)), max(raw, error(hi))
    return worst < 1e-6, ("extrapolated in N: max |B_jj - flux|/|flux| %.1e, j = 1, 3, 5; at N = %d itself %.1e"
                          % (worst, 2 * CHECK_RESOLUTION["n_modes"], raw))


CHECKS = {"roots": check_roots, "inner_products": check_inner_products, "bottom_limit_x1": check_bottom_limit_x1,
          "matching": check_matching, "scattering_energy": check_scattering_energy, "reciprocity": check_reciprocity,
          "radiation_energy": check_radiation_energy}


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
