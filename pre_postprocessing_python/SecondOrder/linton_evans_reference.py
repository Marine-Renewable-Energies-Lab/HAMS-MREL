"""Exact reference for wave scattering by an array of bottom-mounted vertical circular cylinders, to validate HAMS-MREL's
multi-body solve: first-order exciting loads per cylinder, the wave field, and the mean drift loads.

Usage (from this folder):
    python -I linton_evans_reference.py check [NAME ...]
    python -I linton_evans_reference.py reference --out le_reference.csv [--layout square]
    python -I linton_evans_reference.py compare RESULTS [--out le_compare.csv] [--reference le_reference.csv]
    python -I linton_evans_reference.py plot RESULTS --out FOLDER [--reference le_reference.csv]

Method (Linton & Evans 1990, after Spring & Monkmeyer 1974). N cylinders of radii a_j, centres c_j = (x_j, y_j), each from
the sea bed z = -h to above the free surface, in water of depth h. Linear potential flow, time factor exp(-i omega t) in this
module, z up, z = 0 the still waterline. Every boundary is vertical, so the potential separates exactly:

    phi(x, y, z) = -(i g A/omega) Z(z) psi(x, y),     Z(z) = cosh k(z + h)/cosh kh,     omega**2 = g k tanh kh,

with only the propagating vertical mode (no evanescent modes, since no boundary has a horizontal part), the elevation
eta = A psi, and psi solving the 2D Helmholtz equation with d psi/dn = 0 on every circle. The incident wave of heading beta is
psi_I = exp(i k (x cos beta + y sin beta)): unit elevation at the origin. Near cylinder j, with polar coordinates (r_j, theta_j)
about c_j,

    psi_I = I_j sum_m i**m exp(-i m beta) J_m(k r_j) exp(i m theta_j),   I_j = exp(i k (x_j cos beta + y_j sin beta)),

and the scattered wave is psi_S = sum_j sum_n A^j_n Z^j_n H_n(k r_j) exp(i n theta_j), Z^j_n = J_n'(k a_j)/H_n'(k a_j),
H = H^(1) (outgoing in exp(-i omega t)). Graf's addition theorem (DLMF 10.23.7), for r_k < R_jk,

    H_n(k r_j) exp(i n theta_j) = sum_m H_(n-m)(k R_jk) exp(i (n-m) alpha_jk) J_m(k r_k) exp(i m theta_k),

with (R_jk, alpha_jk) the polar form of c_k - c_j, expands the other cylinders' waves about c_k. The no-flux condition on
cylinder k, order m, then gives the linear system of Spring & Monkmeyer (1974) and Linton & Evans (1990)

    A^k_m + sum_(j != k) sum_n Z^j_n H_(n-m)(k R_jk) exp(i (n-m) alpha_jk) A^j_n = -I_k i**m exp(-i m beta),

truncated at |m|, |n| <= M (solve). On cylinder k the total field is, by the Wronskian J_m H_m' - J_m' H_m = 2i/(pi x),

    psi(a_k, theta) = sum_m -A^k_m 2i/(pi k a_k H_m'(k a_k)) exp(i m theta)          (surface_coefficients).

The Graf form is not trusted on its own: check bc_direct sums every cylinder's Hankel series directly (no addition
theorem) at points on each circle and requires d psi/dn = 0 there, and surface_direct requires the surface values above.

First-order loads (exciting_loads). Pressure p = i omega rho phi = rho g A Z psi, force on the body int p n dS with n out of the
fluid (into the body, n = -(cos theta, sin theta, 0) on the wall, as HAMS's mesh normal), moments about (x_k, y_k, 0):

    F_x + i F_y terms: F_x = -rho g A (tanh kh/k) a int psi cos theta d theta, F_y the same with sin theta,
    M_x = +rho g A Zz a int psi sin theta d theta,  M_y = -rho g A Zz a int psi cos theta d theta,
    Zz = int_-h^0 z Z dz = -(cosh kh - 1)/(k**2 cosh kh).

F_z = M_z = 0 (no horizontal surface; n radial). For one cylinder F_x = 4 rho g A tanh(kh)/(k**2 H_1'(ka)), MacCamy & Fuchs
(1954), the closed form of check maccamy_fuchs_force. HAMS's ExcForce.3 holds conj(X/(rho g A)) (exp(+i omega t), phase
relative to the elevation at the origin; truncated_cylinder_reference.py, "Exciting force"), so to_hams conjugates.

Mean drift, near-field (near_field_drift). The time average over a period, n_b out of the body (into the fluid):

    F = -(rho g/4) oint |eta|**2 n_b dl + (rho/4) int int |grad phi|**2 n_b dS,

the waterline term and the quadratic-velocity term (WAMIT v7.5 manual, eqs. 15.54 to 15.57, terms 1 and 2; for a fixed body
the other terms vanish). With phi above, |grad phi|**2 = (g A/omega)**2 [|grad_h psi|**2 Z**2 + k**2 |psi|**2 S**2],
S = sinh k(z + h)/cosh kh, and d psi/dr = 0 on the wall, so the depth integrals are closed form:

    C2 = int Z**2 dz = (h/2 + sinh(2kh)/(4k))/cosh(kh)**2,  S2 = int S**2 dz = (sinh(2kh)/(4k) - h/2)/cosh(kh)**2,

and Gauss-Legendre in z for the moment integrals int z Z**2 dz, int z S**2 dz (smooth, 64 points, exact to rounding). Per
rho g A**2, g/omega**2 = 1/(k tanh kh). The waterline term has no lever arm about (x_k, y_k, 0) for M_x, M_y, and n is radial,
so the mean M_z about the centre is 0 and F_z = 0.

Mean drift, far-field (far_field_drift). The scattered wave behaves as psi_S ~ sqrt(2/(pi k r)) exp(i (k r - pi/4)) K(theta),

    K(theta) = sum_j exp(-i k (x_j cos theta + y_j sin theta)) sum_n A^j_n Z^j_n (-i)**n exp(i n theta),

and the mean horizontal force on the whole array, per rho g A**2, is (Maruo 1960; Newman 1967, in the form of
drift_convergence.maccamy_fuchs_fx, there for beta = 0)

    F = (C_g/C)/(pi k) int |K|**2 (e_beta - e_theta) d theta,   C_g/C = (1 + 2kh/sinh 2kh)/2.

The optical theorem, int |K|**2 d theta = -2 pi Re K(beta), checks energy conservation (check optical_theorem).

Mean drift, control surface (control_surface_drift). The momentum balance of the fluid between the bodies inside a vertical
circle C (radius R, centre P, from the sea bed to the free surface) and C, n out of that fluid at C, gives

    F_i = -oint_C dl [(rho g/4) |eta|**2 n_i + int_-h^0 (-(rho/4) |grad phi|**2 n_i + (rho/2) Re(d_i phi conj(d_n phi))) dz],
    M_z(P) = -oint_C dl int_-h^0 (rho/2) Re(d_theta phi conj(d_r phi)) dz,

(Chen 2007, eqs. (39) to (41), without the hull terms, which vanish for fixed bodies; the sea bed adds nothing horizontal). The
field is summed directly from the Hankel series of every cylinder, without Graf. A circle around one cylinder only must give
its near-field load, and a circle around the whole array the far-field load and sum_k (c_k x F_k + M_z,k) for M_z.

Truncation (default_order). M = ceil(k a_max) + max(M_EXTRA, the order at which (a_j/(R_jk - a_k))**M < 1e-16): the loads
converge fast, but the field of each cylinder's series on its neighbours needs the second. check truncation requires the
loads at M and M + 10 to agree to 1e-12. The unknowns are scaled (Solution), which keeps the condition number of the
system below 10 for the square layout.

Layouts (LAYOUTS). `square`: the HAMS validation case, four cylinders of radius 5 m at the corners of a square, centres
(+-10, +-10) m (centre spacing 20 m = 4a), depth 10 m, headings 0 and 45 deg, the frequencies of FREQUENCIES (k a from about
0.19 to 1.9, below the first irregular frequency of a bottom-mounted circle, k a = j_0,1 = 2.405, the first zero of J_0).

Published values (docs/second_order_validation.md, "Linton & Evans array"). Linton & Evans (1990), Table 1: two cylinders in
the long-wave limit, eight values (check linton_evans_table1, LE_TABLE1). Linton & Evans (1990), Fig. 3: the total drift of the
square of their Fig. 2, the HAMS case scaled by 1/5, digitised (check linton_evans_fig3, LE_FIG3_TOTAL), compared outside the
near-trapping band 1.38 < k a < 2.07, where the published curves are smoother and lower than this reference (its boundary
condition holds to 2e-15 there, the paper's truncation M = 6 gives the same values, and HAMS converges to this reference).
Their Fig. 2 plots, for the cylinder off the diagonal, (|X_x| cos beta + |X_y| sin beta)/|F|, not |X . e_beta|/|F| (compared
by overlay only). Kim (1993), Fig. 3: |F_x| and |F_y| on each of 19 bottom-mounted cylinders in a row from the same exact
method, agreeing there with the hybrid finite elements of Kagemoto & Yue (1986) (check kim_fig3, KIM_FIG3_X, KIM_FIG3_Y).

Checks (CHECKS; `check` runs all, 3 min 46 s on a laptop (measured once), and exits with 1 if one fails): maccamy_fuchs_force,
maccamy_fuchs_drift, depth_integrals, bc_direct, surface_direct, truncation, optical_theorem, near_far, control_surface,
invariance, kim_fig3, linton_evans_table1, linton_evans_fig3. Each docstring states what it compares and its limit.

Units SI. G and RHO are make_mesh.G and make_mesh.RHO, the values of src/WavDynMods.f90.

Sources:
- Linton, C.M. and Evans, D.V. (1990). The interaction of waves with arrays of vertical circular cylinders. J. Fluid Mech. 215,
  549-569. doi:10.1017/S0022112090002750
- Spring, B.H. and Monkmeyer, P.L. (1974). Interaction of plane waves with vertical cylinders. Proc. 14th Int. Conf. Coastal
  Eng., 1828-1847.
- MacCamy, R.C. and Fuchs, R.A. (1954). Wave forces on piles: a diffraction theory. Technical Memorandum 69, Beach Erosion Board.
- Maruo, H. (1960). The drift of a body floating on waves. J. Ship Res. 4, 1-10.
- Newman, J.N. (1967). The drift force and moment on ships in waves. J. Ship Res. 11(1), 51-60.
- Chen, X.B. (2007). Middle-field formulation for the computation of wave-drift loads. J. Eng. Math. 59, 61-82.
- Kim, M.H. (1993). Interaction of waves with N vertical circular cylinders. J. Waterway, Port, Coastal, Ocean Eng. 119(6),
  671-689. doi:10.1061/(ASCE)0733-950X(1993)119:6(671). Fig. 3 (check kim_fig3).
- Kagemoto, H. and Yue, D.K.P. (1986). Wave forces on a platform supported on a large number of floating legs. Proc. 5th
  Int. Offshore Mech. and Arctic Eng. Symp., ASME (cited by Kim 1993 as 1986b, for the hybrid finite elements of Fig. 3).
- NIST DLMF 10.23.7 (Graf's addition theorem).
- WAMIT User Manual, Version 7.5, section 15.8.
"""

import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.special import h1vp, hankel1, jv, jvp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # run with python -I: the sibling modules

import make_mesh

G = make_mesh.G
RHO = make_mesh.RHO
M_EXTRA = 12


def wave_number(omega, depth):
    """k of omega**2 = g k tanh(k h), by Newton to rounding."""
    nu = omega ** 2 / G
    k = max(nu, math.sqrt(nu / depth))
    for _ in range(200):
        t = math.tanh(k * depth)
        step = (k * t - nu) / (t + k * depth * (1.0 - t * t))
        k -= step
        if abs(step) < 1e-15 * k:
            break
    return k


class Array:
    """Bottom-mounted cylinders: centres (N, 2) in m, radii (N,) in m, depth h in m."""

    def __init__(self, centres, radii, depth):
        self.centres = np.asarray(centres, dtype=float).reshape(-1, 2)
        self.radii = np.broadcast_to(np.asarray(radii, dtype=float), (len(self.centres),)).copy()
        self.depth = float(depth)
        for i in range(len(self.centres)):
            for j in range(i):
                gap = np.hypot(*(self.centres[i] - self.centres[j])) - self.radii[i] - self.radii[j]
                if gap <= 0.0:
                    raise ValueError("cylinders %d and %d overlap" % (j + 1, i + 1))

    def __len__(self):
        return len(self.centres)

    def moved(self, shift=(0.0, 0.0), angle=0.0, order=None):
        """The same array rotated by angle (rad) about the origin, then shifted, cylinders in the given order."""
        c, s = math.cos(angle), math.sin(angle)
        rotated = self.centres @ np.array([[c, s], [-s, c]]) + np.asarray(shift, dtype=float)
        idx = np.arange(len(self)) if order is None else np.asarray(order)
        return Array(rotated[idx], self.radii[idx], self.depth)


def _layouts():
    """The layouts: `square` from drift_convergence's LE_* constants (the HAMS case), `single` for the MacCamy-Fuchs limit."""
    import drift_convergence as dc
    return {
        "square": dict(centres=dc.LE_CENTRES, radii=dc.BOTTOM_RADIUS, depth=dc.BOTTOM_DEPTH),
        "single": dict(centres=[(0.0, 0.0)], radii=dc.BOTTOM_RADIUS, depth=dc.BOTTOM_DEPTH),
    }


def _case_constants():
    import drift_convergence as dc
    return dc.LE_FREQUENCIES, dc.LE_HEADINGS, dc.LE_FIELD_POINTS


LAYOUTS = _layouts()
# rad/s; at h = 10 m, a = 5 m these are k a = 0.19 ... 1.84. HEADINGS in degrees. FIELD_POINTS (x, y, 0) in m
FREQUENCIES, HEADINGS, FIELD_POINTS = _case_constants()


# Truncation order of the load curves of `plot`: the loads converge far faster than the field between the cylinders. For the
# square, X and the near-field drift at order 16 are within 7e-15 of the default order (0.40, 1.10, 1.78, 1.90 rad/s, both
# headings; measured), and the default order makes the curves take an hour
LOAD_ORDER = 16


def layout(name):
    spec = LAYOUTS[name]
    return Array(spec["centres"], spec["radii"], spec["depth"])


def default_order(array, k):
    """Truncation order: ceil(k a_max) + M_EXTRA, and at least what makes (a_j/(R_jk - a_k))**M < 1e-16 for the closest pair,
    so that each cylinder's Hankel series converges to rounding on its neighbours (checks bc_direct, control_surface)."""
    q = 0.0
    for j in range(len(array)):
        for i in range(len(array)):
            if i != j:
                R = math.hypot(*(array.centres[i] - array.centres[j]))
                q = max(q, array.radii[j] / (R - array.radii[i]))
    geometric = int(math.ceil(math.log(1e-16) / math.log(q))) if q > 0.0 else 0
    return int(math.ceil(k * array.radii.max())) + max(M_EXTRA, geometric)


class Solution:
    """Scattering coefficients A[j, m + M] of one frequency and heading."""

    def __init__(self, array, omega, beta, order=None):
        self.array = array
        self.omega = float(omega)
        self.beta = float(beta)
        self.k = wave_number(omega, array.depth)
        k = self.k
        self.M = default_order(array, k) if order is None else int(order)
        M = self.M
        m = np.arange(-M, M + 1)
        self.m = m
        n_cyl = len(array)
        size = 2 * M + 1
        ka = k * array.radii
        hp = np.array([h1vp(m, x) for x in ka])   # H_m'(k a_j), (N, 2M+1)
        jp = np.array([jvp(m, x) for x in ka])
        self.Z = jp / hp
        # Scaling. A^j_n grows like |H_n'(k a_j)| (a/R)**|n| at small k a, and the entries Z H_(n-m)(k R) span hundreds of
        # orders of magnitude, so the system is solved for x^j_n = A^j_n/|H_n'(k a_j)|, each row divided by |H_m'(k a_k)|:
        # entry J_n' (|H_n'|/H_n') H_(n-m)(k R) e^(i (n-m) alpha)/|H_m'(k a_k)|, of the size (a/R)**(|m|+|n|).
        scale = np.abs(hp)
        self.scale = scale
        cb, sb = math.cos(beta), math.sin(beta)
        incident = np.exp(1j * k * (array.centres[:, 0] * cb + array.centres[:, 1] * sb))
        self.incident_factor = incident
        rhs = np.zeros((n_cyl, size), dtype=complex)
        for kk in range(n_cyl):
            rhs[kk] = -incident[kk] * (1j ** m) * np.exp(-1j * m * beta) / scale[kk]
        mat = np.eye(n_cyl * size, dtype=complex)
        diff = m[None, :] - m[:, None]   # [row m, column n] = n - m
        for kk in range(n_cyl):
            for jj in range(n_cyl):
                if jj == kk:
                    continue
                d = array.centres[kk] - array.centres[jj]
                R = math.hypot(d[0], d[1])
                alpha = math.atan2(d[1], d[0])
                block = (hankel1(diff, k * R) * np.exp(1j * diff * alpha) * (jp[jj] * scale[jj] / hp[jj])[None, :]
                         / scale[kk][:, None])
                mat[kk * size:(kk + 1) * size, jj * size:(jj + 1) * size] = block
        x = np.linalg.solve(mat, rhs.reshape(-1))
        self.residual = np.linalg.norm(mat @ x - rhs.reshape(-1)) / np.linalg.norm(rhs)
        self._matrix = mat
        self.x = x.reshape(n_cyl, size)
        self.A = self.x * scale
        # A Z = x J' |H'|/H', the scattered-wave coefficients, computed without the large A and small Z
        self.AZ = self.x * jp * scale / hp

    @property
    def condition(self):
        """The 2-norm condition number of the scaled system (an SVD, so computed on demand)."""
        return np.linalg.cond(self._matrix)

    # Field ----------------------------------------------------------------------------------------------------------------
    def field(self, x, y, gradient=False):
        """psi (= eta/A) at points, summed directly from every cylinder's Hankel series (no Graf). With gradient, also
        d psi/dx, d psi/dy. Points inside a cylinder give the series' analytic continuation, not a field."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        k, m = self.k, self.m
        cb, sb = math.cos(self.beta), math.sin(self.beta)
        inc = np.exp(1j * k * (x * cb + y * sb))
        psi = inc.astype(complex)
        gx = 1j * k * cb * inc
        gy = 1j * k * sb * inc
        for j in range(len(self.array)):
            dx = x - self.array.centres[j, 0]
            dy = y - self.array.centres[j, 1]
            r = np.hypot(dx, dy)
            th = np.arctan2(dy, dx)
            coef = self.AZ[j]
            e = np.exp(1j * np.multiply.outer(th, m))
            H = hankel1(m, np.multiply.outer(k * r, np.ones_like(m)))
            psi = psi + np.sum(coef * H * e, axis=-1)
            if gradient:
                Hp = h1vp(m, np.multiply.outer(k * r, np.ones_like(m)))
                dr = np.sum(coef * k * Hp * e, axis=-1)
                dth = np.sum(coef * 1j * m * H * e, axis=-1) / r
                c, s = np.cos(th), np.sin(th)
                gx = gx + dr * c - dth * s
                gy = gy + dr * s + dth * c
        if gradient:
            return psi, gx, gy
        return psi

    # Surface ----------------------------------------------------------------------------------------------------------------
    def surface_coefficients(self, j):
        """Fourier coefficients s_m of psi(a_j, theta) = sum_m s_m exp(i m theta), by the Wronskian."""
        x = self.k * self.array.radii[j]
        return -self.x[j] * self.scale[j] * 2j / (math.pi * x * h1vp(self.m, x))

    def surface(self, j, theta):
        s = self.surface_coefficients(j)
        e = np.exp(1j * np.multiply.outer(theta, self.m))
        return e @ s, e @ (1j * self.m * s)   # psi, d psi/d theta


def depth_integrals(k, h, points=64):
    """int_-h^0 of Z, z Z, Z**2, S**2, z Z**2, z S**2 (Z = cosh k(z+h)/cosh kh, S = sinh k(z+h)/cosh kh)."""
    ch = math.cosh(k * h)
    z1 = math.tanh(k * h) / k
    zz = -(ch - 1.0) / (k * k * ch)
    c2 = (h / 2.0 + math.sinh(2 * k * h) / (4 * k)) / ch ** 2
    s2 = (math.sinh(2 * k * h) / (4 * k) - h / 2.0) / ch ** 2
    t, w = np.polynomial.legendre.leggauss(points)
    z = -h / 2.0 * (1.0 - t)
    w = w * h / 2.0
    Z = np.cosh(k * (z + h)) / ch
    S = np.sinh(k * (z + h)) / ch
    return dict(Z=z1, zZ=zz, Z2=c2, S2=s2, zZ2=float(np.sum(w * z * Z * Z)), zS2=float(np.sum(w * z * S * S)),
                Z2_quad=float(np.sum(w * Z * Z)), S2_quad=float(np.sum(w * S * S)), zZ_quad=float(np.sum(w * z * Z)))


LOADS = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


def exciting_loads(sol, j):
    """X_1 ... X_6 of cylinder j over rho g A (moments in m, about (x_j, y_j, 0)), exp(-i omega t), per unit elevation at the
    origin."""
    s = sol.surface_coefficients(j)
    M = sol.M
    a = sol.array.radii[j]
    d = depth_integrals(sol.k, sol.array.depth)
    sp, sm = s[M + 1], s[M - 1]
    int_cos = math.pi * (sp + sm)          # int psi cos theta d theta
    int_sin = 1j * math.pi * (sp - sm)     # int psi sin theta d theta: sin = (e^(i t) - e^(-i t))/(2i)
    return np.array([-d["Z"] * a * int_cos, -d["Z"] * a * int_sin, 0.0, d["zZ"] * a * int_sin, -d["zZ"] * a * int_cos, 0.0],
                    dtype=complex)


def to_hams(value):
    """ExcForce.3's NREL + i NIMG of a load of exciting_loads (module docstring)."""
    return np.conj(value)


def _theta_grid(M, extra=8):
    n = 4 * M + 4 * extra + 8
    return 2.0 * math.pi * np.arange(n) / n, 2.0 * math.pi / n


def near_field_drift(sol, j):
    """Mean F_x ... M_z of cylinder j over rho g A**2 (moments in m about (x_j, y_j, 0)), fixed body."""
    k, h = sol.k, sol.array.depth
    a = sol.array.radii[j]
    d = depth_integrals(k, h)
    t = math.tanh(k * h)
    theta, dth = _theta_grid(sol.M)
    psi, dpsi = sol.surface(j, theta)
    p2 = np.abs(psi) ** 2
    h2 = np.abs(dpsi / a) ** 2
    c, s = np.cos(theta), np.sin(theta)
    # force density on the wall per unit length of waterline: waterline + velocity term
    f = -0.25 * p2 + 0.25 * (h2 * d["Z2"] + k * k * p2 * d["S2"]) / (k * t)
    fm = 0.25 * (h2 * d["zZ2"] + k * k * p2 * d["zS2"]) / (k * t)   # z-weighted velocity term
    out = np.zeros(6)
    out[0] = np.sum(f * c) * a * dth
    out[1] = np.sum(f * s) * a * dth
    out[3] = -np.sum(fm * s) * a * dth
    out[4] = np.sum(fm * c) * a * dth
    return out


def near_field_terms(sol, j):
    """The waterline and velocity terms of F_x, F_y of cylinder j, as near_field_drift sums them."""
    k, h = sol.k, sol.array.depth
    a = sol.array.radii[j]
    d = depth_integrals(k, h)
    t = math.tanh(k * h)
    theta, dth = _theta_grid(sol.M)
    psi, dpsi = sol.surface(j, theta)
    p2 = np.abs(psi) ** 2
    h2 = np.abs(dpsi / a) ** 2
    wl = -0.25 * p2
    vel = 0.25 * (h2 * d["Z2"] + k * k * p2 * d["S2"]) / (k * t)
    c, s = np.cos(theta), np.sin(theta)
    return dict(Fx_wl=np.sum(wl * c) * a * dth, Fx_vel=np.sum(vel * c) * a * dth,
                Fy_wl=np.sum(wl * s) * a * dth, Fy_vel=np.sum(vel * s) * a * dth)


def kochin(sol, theta):
    theta = np.asarray(theta, dtype=float)
    k, m = sol.k, sol.m
    out = np.zeros(theta.shape, dtype=complex)
    e = np.exp(1j * np.multiply.outer(theta, m)) * ((-1j) ** m)
    for j in range(len(sol.array)):
        cx, cy = sol.array.centres[j]
        phase = np.exp(-1j * k * (cx * np.cos(theta) + cy * np.sin(theta)))
        out += phase * (e @ sol.AZ[j])
    return out


def far_field_drift(sol, count=None):
    """Mean F_x, F_y of the whole array over rho g A**2, from the Kochin function."""
    k, h = sol.k, sol.array.depth
    extent = sol.k * (np.hypot(sol.array.centres[:, 0], sol.array.centres[:, 1]) + sol.array.radii).max()
    n = count or int(8 * (extent + sol.M) + 64)
    theta = 2.0 * math.pi * np.arange(n) / n
    K2 = np.abs(kochin(sol, theta)) ** 2
    ratio = 0.5 * (1.0 + 2 * k * h / math.sinh(2 * k * h))
    w = 2.0 * math.pi / n
    fx = ratio / (math.pi * k) * np.sum(K2 * (math.cos(sol.beta) - np.cos(theta))) * w
    fy = ratio / (math.pi * k) * np.sum(K2 * (math.sin(sol.beta) - np.sin(theta))) * w
    return np.array([fx, fy])


def control_surface_drift(sol, centre, radius, count=None):
    """Mean F_x, F_y and M_z about centre, over rho g A**2, from the momentum balance on a vertical circle (module docstring),
    the field summed directly."""
    k, h = sol.k, sol.array.depth
    d = depth_integrals(k, h)
    t = math.tanh(k * h)
    n = count or int(8 * (k * radius + sol.M) + 64)
    theta = 2.0 * math.pi * np.arange(n) / n
    c, s = np.cos(theta), np.sin(theta)
    x = centre[0] + radius * c
    y = centre[1] + radius * s
    psi, gx, gy = sol.field(x, y, gradient=True)
    gn = gx * c + gy * s
    gt = (-gx * s + gy * c) * radius   # d psi/d theta
    p2 = np.abs(psi) ** 2
    h2 = np.abs(gx) ** 2 + np.abs(gy) ** 2
    vel2 = (h2 * d["Z2"] + k * k * p2 * d["S2"]) / (k * t)
    dl = radius * 2.0 * math.pi / n
    fx = -np.sum(0.25 * p2 * c - 0.25 * vel2 * c + 0.5 * np.real(gx * np.conj(gn)) * d["Z2"] / (k * t)) * dl
    fy = -np.sum(0.25 * p2 * s - 0.25 * vel2 * s + 0.5 * np.real(gy * np.conj(gn)) * d["Z2"] / (k * t)) * dl
    mz = -np.sum(0.5 * np.real(gt * np.conj(gn)) * d["Z2"] / (k * t)) * dl
    return np.array([fx, fy, mz])


def array_mz(sol, near=None):
    """Mean M_z of the whole array about the origin, from the near-field loads of the cylinders."""
    total = 0.0
    for j in range(len(sol.array)):
        f = near_field_drift(sol, j) if near is None else near[j]
        cx, cy = sol.array.centres[j]
        total += cx * f[1] - cy * f[0] + f[5]
    return total


# Checks ---------------------------------------------------------------------------------------------------------------------
def _rel(a, b, scale=None):
    a = np.asarray(a)
    b = np.asarray(b)
    s = np.max(np.abs(b)) if scale is None else scale
    return float(np.max(np.abs(a - b)) / s)


def check_maccamy_fuchs_force():
    """One cylinder: X_1 against 4 tanh(kh)/(k**2 H_1'(ka)) (MacCamy & Fuchs 1954), complex, and X_5 against its depth
    moment, 1e-13, X_2 and X_4 zero; and at heading 90 deg X_2 against the same value and X_4 against minus that moment
    (the sin part of exciting_loads)."""
    arr = layout("single")
    worst = 0.0
    for omega in FREQUENCIES:
        sol = Solution(arr, omega, 0.0)
        k, h, a = sol.k, arr.depth, arr.radii[0]
        exact = 4.0 * math.tanh(k * h) / (k * k * h1vp(1, k * a))
        X = exciting_loads(sol, 0)
        my = exact * depth_integrals(k, h)["zZ"] / (math.tanh(k * h) / k)
        worst = max(worst, abs(X[0] - exact) / abs(exact), abs(X[4] - my) / abs(my), abs(X[1]) / abs(exact),
                    abs(X[3]) / abs(my))
        Y = exciting_loads(Solution(arr, omega, math.pi / 2), 0)
        worst = max(worst, abs(Y[1] - exact) / abs(exact), abs(Y[3] + my) / abs(my), abs(Y[0]) / abs(exact))
    return worst, 1e-13


def check_maccamy_fuchs_drift():
    """One cylinder: near-field F_x against far-field F_x of drift_convergence.maccamy_fuchs_fx (an independent code, the one
    the bottom-mounted HAMS case was checked with), and against this module's far-field, 1e-10; F_y zero."""
    import drift_convergence
    arr = layout("single")
    worst = 0.0
    for omega in FREQUENCIES:
        sol = Solution(arr, omega, 0.0)
        near = near_field_drift(sol, 0)
        other = drift_convergence.maccamy_fuchs_fx(omega, arr.radii[0], arr.depth, count=2048)
        far = far_field_drift(sol)
        worst = max(worst, abs(near[0] - other) / abs(other), abs(near[0] - far[0]) / abs(far[0]), abs(near[1]) / abs(other))
    return worst, 1e-10


def _cases():
    arr = layout("square")
    for omega in (0.40, 0.95, 1.60, 1.90):
        for beta in HEADINGS:
            yield arr, omega, math.radians(beta)
    # an irregular array of unequal cylinders, off-centre, so no symmetry hides an error
    odd = Array([(3.0, -2.0), (17.0, 4.0), (6.0, 21.0)], [5.0, 3.0, 4.0], 12.0)
    for omega in (0.7, 1.3):
        for beta in (20.0, 200.0):
            yield odd, omega, math.radians(beta)


def check_bc_direct():
    """d psi/dn on every cylinder, summing every cylinder's Hankel series directly (no Graf), over k: 1e-11."""
    worst = 0.0
    for arr, omega, beta in _cases():
        sol = Solution(arr, omega, beta)
        theta = np.linspace(0, 2 * math.pi, 97)[:-1]
        for j in range(len(arr)):
            c, s = np.cos(theta), np.sin(theta)
            x = arr.centres[j, 0] + arr.radii[j] * c
            y = arr.centres[j, 1] + arr.radii[j] * s
            _, gx, gy = sol.field(x, y, gradient=True)
            worst = max(worst, float(np.max(np.abs(gx * c + gy * s))) / sol.k)
    return worst, 1e-11


def check_surface_direct():
    """psi on every cylinder from the Wronskian form against the direct sum: 1e-11."""
    worst = 0.0
    for arr, omega, beta in _cases():
        sol = Solution(arr, omega, beta)
        theta = np.linspace(0, 2 * math.pi, 61)[:-1]
        for j in range(len(arr)):
            x = arr.centres[j, 0] + arr.radii[j] * np.cos(theta)
            y = arr.centres[j, 1] + arr.radii[j] * np.sin(theta)
            worst = max(worst, float(np.max(np.abs(sol.field(x, y) - sol.surface(j, theta)[0]))))
    return worst, 1e-11


def check_truncation():
    """Loads and drift at M and M + 10: 1e-12 relative to the largest of each kind."""
    worst = 0.0
    for arr, omega, beta in _cases():
        s1 = Solution(arr, omega, beta)
        s2 = Solution(arr, omega, beta, order=s1.M + 10)
        x1 = np.array([exciting_loads(s1, j) for j in range(len(arr))])
        x2 = np.array([exciting_loads(s2, j) for j in range(len(arr))])
        d1 = np.array([near_field_drift(s1, j) for j in range(len(arr))])
        d2 = np.array([near_field_drift(s2, j) for j in range(len(arr))])
        worst = max(worst, _rel(x1, x2), _rel(d1, d2))
    return worst, 1e-12


def check_optical_theorem():
    """int |K|**2 = -2 pi Re K(beta) (energy conservation), relative: 1e-12."""
    worst = 0.0
    for arr, omega, beta in _cases():
        sol = Solution(arr, omega, beta)
        n = int(8 * (sol.k * 40 + sol.M) + 256)
        theta = 2 * math.pi * np.arange(n) / n
        lhs = np.sum(np.abs(kochin(sol, theta)) ** 2) * 2 * math.pi / n
        rhs = -2 * math.pi * np.real(kochin(sol, np.array([beta]))[0])
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
    return worst, 1e-12


def check_near_far():
    """Near-field F_x, F_y summed over the cylinders against far-field of the array: 1e-10 relative to |F|."""
    worst = 0.0
    for arr, omega, beta in _cases():
        sol = Solution(arr, omega, beta)
        near = sum(near_field_drift(sol, j)[:2] for j in range(len(arr)))
        far = far_field_drift(sol)
        worst = max(worst, _rel(near, far, np.hypot(*far)))
    return worst, 1e-10


def check_control_surface():
    """Control surface around each cylinder alone (radius a + 40% of the gap to the nearest neighbour) against its near-field
    F_x, F_y, M_z = 0; and around the whole array against far-field F_x, F_y and the near-field M_z of the array: 1e-9 relative
    to the largest |F| (forces) or |F| times the array radius (M_z)."""
    worst = 0.0
    for arr, omega, beta in _cases():
        sol = Solution(arr, omega, beta)
        near = [near_field_drift(sol, j) for j in range(len(arr))]
        scale = max(np.hypot(f[0], f[1]) for f in near)
        for j in range(len(arr)):
            gaps = [np.hypot(*(arr.centres[i] - arr.centres[j])) - arr.radii[i] - arr.radii[j]
                    for i in range(len(arr)) if i != j]
            r = arr.radii[j] + 0.4 * (min(gaps) if gaps else arr.radii[j])
            cs = control_surface_drift(sol, arr.centres[j], r)
            worst = max(worst, abs(cs[0] - near[j][0]) / scale, abs(cs[1] - near[j][1]) / scale,
                        abs(cs[2]) / (scale * arr.radii[j]))
        R = (np.hypot(arr.centres[:, 0], arr.centres[:, 1]) + arr.radii).max() * 1.3
        cs = control_surface_drift(sol, (0.0, 0.0), R)
        far = far_field_drift(sol)
        total = np.hypot(*far)
        worst = max(worst, abs(cs[0] - far[0]) / total, abs(cs[1] - far[1]) / total,
                    abs(cs[2] - array_mz(sol, near)) / (total * R))
    return worst, 1e-9


def check_invariance():
    """The square rotated by 90 deg with the heading, and its cylinders renumbered: the same loads per cylinder (forces
    rotated with it), 1e-11."""
    arr = layout("square")
    worst = 0.0
    for omega in (0.6, 1.6):
        beta = math.radians(15.0)
        s0 = Solution(arr, omega, beta)
        # rotate by 90 deg: cylinder j goes to the position of j + 1 (the layout is listed counter-clockwise)
        order = [3, 0, 1, 2]
        rot = arr.moved(angle=0.0, order=order)   # same positions, renumbered
        s1 = Solution(rot, omega, beta)
        for new, old in enumerate(order):
            worst = max(worst, _rel(exciting_loads(s1, new), exciting_loads(s0, old)),
                        _rel(near_field_drift(s1, new), near_field_drift(s0, old)))
        s2 = Solution(arr, omega, beta + math.pi / 2)   # heading rotated by 90 deg
        for j in range(4):
            jj = (j + 1) % 4   # cylinder j's position rotated by 90 deg is cylinder j + 1's
            x0 = exciting_loads(s0, j)
            x2 = exciting_loads(s2, jj)
            # F rotated by 90 deg: (Fx, Fy) -> (-Fy, Fx); (Mx, My) likewise. Phase: the incident wave at the rotated
            # position has the same phase, since the origin is the centre of rotation
            expect = np.array([-x0[1], x0[0], 0, -x0[4], x0[3], 0])
            worst = max(worst, _rel(x2, expect, np.max(np.abs(x0))))
            d0 = near_field_drift(s0, j)
            d2 = near_field_drift(s2, jj)
            worst = max(worst, _rel(d2, np.array([-d0[1], d0[0], 0, -d0[4], d0[3], 0]), np.max(np.abs(d0))))
    return worst, 1e-11


def check_depth_integrals():
    """The closed forms of int Z, int z Z, int Z**2, int S**2 against Gauss-Legendre: 1e-13."""
    worst = 0.0
    for k, h in ((0.05, 10.0), (0.4, 10.0), (1.2, 12.0)):
        d = depth_integrals(k, h)
        worst = max(worst, abs(d["Z2"] - d["Z2_quad"]) / d["Z2"], abs(d["S2"] - d["S2_quad"]) / d["S2"],
                    abs(d["zZ"] - d["zZ_quad"]) / abs(d["zZ"]))
    return worst, 1e-13


# Kim (1993), Fig. 3: |F_x| and |F_y|/(pi rho g a h A) on each of 19 bottom-mounted cylinders in a row along x, centre spacing
# 4a, h = 12a, k a = 0.5, heading 45 deg (the case of Kagemoto & Yue 1986, whose hybrid finite elements agree with it there).
# Digitised in November 2026 from the scanned figure: the centre of each marker's vertical extent, against the axis ticks
# (81 px per 0.01 on (a), 56 px per 0.01 on (b), at 260 dpi). A marker is about 0.0015 high in these units; the scan is
# slightly tilted, and the differences from the reference grow from cylinder 1 to 19 on both panels.
KIM_FIG3_X = (0.1210, 0.1018, 0.1119, 0.0942, 0.1060, 0.0941, 0.1041, 0.1016, 0.1041, 0.1105, 0.1010, 0.1141, 0.0954, 0.1109,
              0.0925, 0.1076, 0.0974, 0.1056, 0.1099)
KIM_FIG3_Y = (0.1098, 0.1222, 0.1264, 0.1284, 0.1280, 0.1276, 0.1262, 0.1255, 0.1247, 0.1247, 0.1248, 0.1250, 0.1258, 0.1257,
              0.1272, 0.1252, 0.1283, 0.1225, 0.1344)


# Linton & Evans (1990), Table 1 ("based on numerical computation of (2.15)"): two cylinders of equal radius a, centre
# spacing R = a/lambda, the line of centres along y (alpha_12 = -pi/2), heading beta = pi/6, in the long-wave limit k a -> 0:
# |X_x|/|F| and |X_y|/|F|, F the force on one cylinder alone (their eq. 3.3). Printed to three decimals.
LE_TABLE1 = {0.1: (0.875, 0.495), 0.2: (0.902, 0.481), 0.3: (0.953, 0.459), 0.4: (1.046, 0.435)}
# k a of the long-wave limit: the values change by under 3e-5 between k a = 3e-3 and 1e-3; smaller k a overflows the
# Hankel functions of the high orders
LE_TABLE1_KA = 1e-3
LE_TABLE1_ORDER = 16   # the loads agree to 1e-5 with order 24


def linton_evans_table1(lam):
    """|X_x|/|F|, |X_y|/|F| on either cylinder of the pair of Table 1 of Linton & Evans (1990)."""
    a, h = 1.0, 10.0
    k = LE_TABLE1_KA / a
    omega = math.sqrt(G * k * math.tanh(k * h))
    arr = Array([(0.0, 0.5 * a / lam), (0.0, -0.5 * a / lam)], a, h)
    sol = Solution(arr, omega, math.pi / 6.0, order=LE_TABLE1_ORDER)
    F = 4.0 * math.tanh(k * h) / (k * k * abs(h1vp(1, k * a)))
    X = exciting_loads(sol, 0)
    return abs(X[0]) / F, abs(X[1]) / F


def check_linton_evans_table1():
    """Published values: Linton & Evans (1990), Table 1, all eight values against the three printed decimals: within 1e-3,
    one unit of the last printed digit. Measured 4.98e-4: lambda = 0.3, |X_x|/|F| = 0.95350 against 0.953, on the rounding
    boundary; the other seven round to the printed values."""
    worst = 0.0
    for lam, values in LE_TABLE1.items():
        got = linton_evans_table1(lam)
        worst = max(worst, max(abs(g - v) for g, v in zip(got, values)))
    return worst, 1e-3


# Linton & Evans (1990), Fig. 3: the total drift force on the four cylinders of their Fig. 2 (a/h = 1/2 with h half the side
# of the square, R/h = 2, beta = pi/4) in the direction of wave advance, over the drift force f_1 on one cylinder alone, against
# k a. Digitised in November 2026 from the published figure: the centre of the uppermost curve (the total is above the others
# for k a >= 0.8) in each pixel column, against the axis ticks (123.5 px per unit at 300 dpi), at k a = 0.80, 0.85, ..., 3.10.
# The layout is the HAMS case scaled by 1/5 (their cylinders 1 to 4 are this module's 2, 1, 4, 3).
LE_FIG3_KA = tuple(round(0.80 + 0.05 * i, 2) for i in range(47))
LE_FIG3_TOTAL = (2.685, 3.025, 3.187, 3.183, 3.057, 2.883, 2.741, 2.656, 2.668, 2.741, 2.870, 3.028, 3.336, 3.534, 3.870,
                 4.121, 4.166, 4.053, 3.729, 3.522, 3.328, 3.389, 3.551, 3.822, 4.081, 4.291, 4.381, 4.336, 4.211, 4.045,
                 3.899, 3.785, 3.729, 3.741, 3.834, 4.020, 4.146, 4.344, 4.522, 4.684, 4.810, 4.887, 4.891, 4.810, 4.623,
                 4.474, 4.121)
# The near-trapping band of this layout, where the published curves are smoother and lower than this reference (module
# docstring, "Published values"); outside it the digitised total is compared
LE_FIG3_BAND = (1.38, 2.07)


def linton_evans_fig3_total(ka_values=LE_FIG3_KA):
    """f/f_1 of the whole array of Linton & Evans (1990), Fig. 3, in the direction of wave advance, at each k a."""
    a, h, depth, beta = 1.0, 2.0, 10.0, math.pi / 4.0
    arr = Array([(-h, h), (h, h), (h, -h), (-h, -h)], a, depth)
    one = Array([(0.0, 0.0)], a, depth)
    e = np.array([math.cos(beta), math.sin(beta)])
    out = []
    for ka in ka_values:
        k = ka / a
        omega = math.sqrt(G * k * math.tanh(k * depth))
        f1 = float(np.dot(near_field_drift(Solution(one, omega, beta), 0)[:2], e))
        # order 16: the loads, which need far fewer orders than the field between the cylinders, agree with the default
        # order to 1e-12 (check truncation), and the default takes minutes for this list
        sol = Solution(arr, omega, beta, order=16)
        out.append(sum(float(np.dot(near_field_drift(sol, j)[:2], e)) for j in range(4)) / f1)
    return np.array(out)


def check_linton_evans_fig3():
    """Published values: the total drift of Linton & Evans (1990), Fig. 3, digitised, outside the near-trapping band
    LE_FIG3_BAND: within 0.15 (measured 0.112 at k a = 2.55, where a cylinder's curve touches the total and the digitiser's
    column centre moves; 0.05, 1.5 %, elsewhere)."""
    ka = np.array(LE_FIG3_KA)
    keep = (ka < LE_FIG3_BAND[0]) | (ka > LE_FIG3_BAND[1])
    ref = linton_evans_fig3_total(ka[keep])
    return float(np.max(np.abs(ref - np.array(LE_FIG3_TOTAL)[keep]))), 0.15


def kim_fig3():
    """The reference's |F_x|, |F_y|/(pi rho g a h A) for the case of KIM_FIG3_X (a = 1 m, h = 12 m)."""
    a, h, ka = 1.0, 12.0, 0.5
    k = ka / a
    omega = math.sqrt(G * k * math.tanh(k * h))
    arr = Array([(4.0 * a * (i - 9), 0.0) for i in range(19)], a, h)
    sol = Solution(arr, omega, math.radians(45.0))
    X = np.array([exciting_loads(sol, j) for j in range(19)])
    return np.abs(X[:, 0]) / (math.pi * a * h), np.abs(X[:, 1]) / (math.pi * a * h)


def check_kim_fig3():
    """Published values: the 19 cylinders of Kim (1993), Fig. 3, against the digitised figure, largest absolute difference:
    within 0.0015, the height of a marker (measured 0.0011, at cylinder 19 on (b))."""
    fx, fy = kim_fig3()
    return float(max(np.max(np.abs(fx - KIM_FIG3_X)), np.max(np.abs(fy - KIM_FIG3_Y)))), 1.5e-3


CHECKS = dict(maccamy_fuchs_force=check_maccamy_fuchs_force, maccamy_fuchs_drift=check_maccamy_fuchs_drift,
              depth_integrals=check_depth_integrals, bc_direct=check_bc_direct, surface_direct=check_surface_direct,
              truncation=check_truncation, optical_theorem=check_optical_theorem, near_far=check_near_far,
              control_surface=check_control_surface, invariance=check_invariance, kim_fig3=check_kim_fig3,
              linton_evans_table1=check_linton_evans_table1, linton_evans_fig3=check_linton_evans_fig3)


def run_checks(names=None):
    failed = 0
    for name in names or CHECKS:
        value, limit = CHECKS[name]()
        ok = value <= limit
        failed += not ok
        print("%-22s %10.2e  (limit %.1e)  %s" % (name, value, limit, "ok" if ok else "FAIL"))
    return failed


# Reference rows, HAMS outputs, comparison ----------------------------------------------------------------------------------
COLUMNS = ["kind", "omega", "heading", "body", "ka", "X1_re", "X1_im", "X2_re", "X2_im", "X4_re", "X4_im", "X5_re", "X5_im",
           "Fx", "Fy", "Fz", "Mx", "My", "Mz", "point", "x", "y", "eta_re", "eta_im"]


def reference_rows(array=None, frequencies=None, headings=None, points=None):
    """Rows of the reference in HAMS's conventions: kind `body` (X_j of ExcForce.3 as NREL + i NIMG, the near-field drift of
    MeanDriftPressure.9, per rho g A**2, moments about the body centre), kind `array` (far-field F_x, F_y and M_z about the
    origin of MeanDriftMomentum.8, M_z from the near-field loads), kind `field` (the total elevation of the diffraction rows
    of PressureElevation.6p at z = 0, conj(psi))."""
    array = array or layout("square")
    rows = []
    for omega in frequencies or FREQUENCIES:
        for heading in headings or HEADINGS:
            sol = Solution(array, omega, math.radians(heading))
            ka = sol.k * array.radii.max()
            near = []
            for j in range(len(array)):
                X = to_hams(exciting_loads(sol, j))
                d = near_field_drift(sol, j)
                near.append(d)
                row = dict(kind="body", omega=omega, heading=heading, body=j + 1, ka=ka)
                for idx in (1, 2, 4, 5):
                    row["X%d_re" % idx] = X[idx - 1].real
                    row["X%d_im" % idx] = X[idx - 1].imag
                for name, value in zip(LOADS, d):
                    row[name] = value
                rows.append(row)
            far = far_field_drift(sol)
            rows.append(dict(kind="array", omega=omega, heading=heading, body=0, ka=ka, Fx=far[0], Fy=far[1],
                             Mz=array_mz(sol, near)))
            pts = np.array(points or FIELD_POINTS, dtype=float)
            eta = np.conj(sol.field(pts[:, 0], pts[:, 1]))
            for n, (pt, e) in enumerate(zip(pts, eta), start=1):
                rows.append(dict(kind="field", omega=omega, heading=heading, ka=ka, point=n, x=pt[0], y=pt[1],
                                 eta_re=e.real, eta_im=e.imag))
    return rows


def write_rows(path, rows, columns=None):
    columns = columns or COLUMNS
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({key: ("%.12g" % value if isinstance(value, float) else value) for key, value in row.items()})


def read_rows(path):
    rows = []
    for row in csv.DictReader(open(path)):
        rows.append({key: value if key in ("kind", "case") else float(value) for key, value in row.items() if value != ""})
    return rows


def _key(omega, heading):
    return (round(omega, 4), round(heading, 3))


def read_hams(folder):
    """From a case's output/Wamit_format folder: X (NREL + i NIMG of ExcForce.3), the near-field drift (Re of
    MeanDriftPressure.9), the far-field drift (Re of MeanDriftMomentum.8, I = 1, 2, 6) and the total elevation at the field
    points (the diffraction rows of PressureElevation.6p), each by (omega, heading) and then index or point."""
    X, near, far, field = {}, {}, {}, {}
    for line in open(os.path.join(folder, "ExcForce.3")):
        c = line.split()
        if len(c) >= 7:
            X.setdefault(_key(float(c[0]), float(c[1])), {})[int(c[2])] = complex(float(c[5]), float(c[6]))
    for name, table in (("MeanDriftPressure.9", near), ("MeanDriftMomentum.8", far)):
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            continue
        for line in open(path):
            c = line.split()
            if c and not line.startswith("#") and abs(float(c[1]) - float(c[2])) < 1e-9:
                table.setdefault(_key(float(c[0]), float(c[1])), {})[int(c[3])] = float(c[6])
    path = os.path.join(folder, "PressureElevation.6p")
    if os.path.exists(path):
        for line in open(path):
            c = line.split()
            # the diffraction rows: OUFR, heading, point, Re, Im (OutputPressureElevation_DiffractionMulti in
            # src/PressureElevationMulti.f90); the radiation rows have 14 fields and their continuation lines 12
            if len(c) == 5:
                field.setdefault(_key(float(c[0]), float(c[1])), {})[int(c[2])] = complex(float(c[3]), float(c[4]))
    return X, near, far, field


def compare_case(folder, ref_rows, case=""):
    """Errors of one HAMS case against the reference rows, per frequency, heading and body (or the array):
    X_force_err = max over X_1, X_2 of |X_hams - X_ref| / max(|X_1|, |X_2|), X_moment_err the same for X_4, X_5 (complex, so
    phase errors count); near_F_err = max over F_x, F_y of |F_hams - F_ref| over the largest |F_x|, |F_y| of all the cylinders
    at that frequency and heading (one cylinder's drift can be near zero), near_M_err likewise for M_x, M_y; near_FzMz_abs =
    max(|F_z|, |M_z|) of HAMS over the force scale (exact zero); far_err over F_x, F_y
    of the array over |F|, and Mz_err = |M_z,hams - M_z,ref| / (|F| 10 m); eta_err = max |eta_hams - eta_ref| over the field
    points, per unit amplitude. Also the signed percentages of |X_1|, |X_2| and near F_x, F_y where the reference is not small."""
    X, near, far, field = read_hams(folder)
    out = []
    by = {}
    for row in ref_rows:
        by.setdefault((row["kind"], _key(row["omega"], row["heading"])), []).append(row)
    for (kind, key), rows in sorted(by.items()):
        if kind == "body" and key in X:
            fsd = max(max(abs(r["Fx"]), abs(r["Fy"])) for r in rows)
            msd = max(max(abs(r["Mx"]), abs(r["My"])) for r in rows)
            for row in rows:
                b = int(row["body"])
                ref = {j: complex(row["X%d_re" % j], row["X%d_im" % j]) for j in (1, 2, 4, 5)}
                got = {j: X[key].get(6 * (b - 1) + j, complex("nan")) for j in (1, 2, 4, 5)}
                fs = max(abs(ref[1]), abs(ref[2]))
                ms = max(abs(ref[4]), abs(ref[5]))
                e = dict(case=case, kind="body", omega=key[0], heading=key[1], body=b, ka=row["ka"],
                         X_force_err=max(abs(got[j] - ref[j]) for j in (1, 2)) / fs,
                         X_moment_err=max(abs(got[j] - ref[j]) for j in (4, 5)) / ms)
                for j in (1, 2):
                    e["X%d_mod_pct" % j] = (100.0 * (abs(got[j]) - abs(ref[j])) / abs(ref[j])
                                            if abs(ref[j]) > 1e-3 * fs else "")
                if key in near:
                    n = {i: near[key].get(6 * (b - 1) + i, float("nan")) for i in range(1, 7)}
                    e["near_F_err"] = max(abs(n[1] - row["Fx"]), abs(n[2] - row["Fy"])) / fsd
                    e["near_M_err"] = max(abs(n[4] - row["Mx"]), abs(n[5] - row["My"])) / msd
                    e["near_FzMz_abs"] = max(abs(n[3]), abs(n[6])) / fsd
                    for i, name in ((1, "Fx"), (2, "Fy")):
                        e["near_%s_pct" % name] = (100.0 * (n[i] - row[name]) / abs(row[name])
                                                   if abs(row[name]) > 1e-2 * fsd else "")
                        e["near_%s" % name] = n[i]
                        e["ref_%s" % name] = row[name]
                out.append(e)
        elif kind == "array" and key in far:
            row = rows[0]
            f = far[key]
            scale = math.hypot(row["Fx"], row["Fy"])
            nan = float("nan")
            out.append(dict(case=case, kind="array", omega=key[0], heading=key[1], body=0, ka=row["ka"],
                            far_err=max(abs(f.get(1, nan) - row["Fx"]), abs(f.get(2, nan) - row["Fy"])) / scale,
                            far_Fx=f.get(1), far_Fy=f.get(2), ref_Fx=row["Fx"], ref_Fy=row["Fy"], far_Mz=f.get(6),
                            ref_Mz=row["Mz"], Mz_err=abs(f.get(6, nan) - row["Mz"]) / (scale * 10.0)))
        elif kind == "field" and key in field:
            err = max(abs(field[key].get(int(r["point"]), complex("nan")) - complex(r["eta_re"], r["eta_im"])) for r in rows)
            out.append(dict(case=case, kind="field", omega=key[0], heading=key[1], body=0, ka=rows[0]["ka"], eta_err=err))
    return out


COMPARE_COLUMNS = ["case", "kind", "omega", "heading", "body", "ka", "X_force_err", "X_moment_err", "X1_mod_pct", "X2_mod_pct",
                   "near_F_err", "near_M_err", "near_FzMz_abs", "near_Fx_pct", "near_Fy_pct", "near_Fx", "near_Fy",
                   "far_err", "far_Fx", "far_Fy", "ref_Fx", "ref_Fy", "far_Mz", "ref_Mz", "Mz_err", "eta_err"]
ERRORS = (("X forces", "X_force_err"), ("X moments", "X_moment_err"), ("near-field F", "near_F_err"),
          ("near-field M", "near_M_err"), ("far-field F", "far_err"), ("far-field M_z", "Mz_err"), ("elevation", "eta_err"))


def case_folders(results):
    """(name, Wamit_format folder) of the l_square_* cases under RESULTS that have output."""
    out = []
    for name in sorted(os.listdir(results)):
        folder = os.path.join(results, name, "output", "Wamit_format")
        if name.startswith("l_square_") and os.path.exists(os.path.join(folder, "ExcForce.3")):
            out.append((name, folder))
    return out


def worst(rows, case, key):
    values = [r[key] for r in rows if r["case"] == case and r.get(key) not in (None, "")]
    return max(values) if values else None


def summary(rows):
    """Per case, the largest of each error over frequencies, headings and cylinders, as a Markdown table."""
    cases = sorted(set(r["case"] for r in rows))
    print("| Case | " + " | ".join(label for label, _ in ERRORS) + " |")
    print("| :--- |" + " ---: |" * len(ERRORS))
    for case in cases:
        cells = []
        for _, key in ERRORS:
            value = worst(rows, case, key)
            cells.append("-" if value is None else "%.2e" % value)
        print("| %s | %s |" % (case, " | ".join(cells)))


def _hams_points(table, index, heading, depth, radius):
    pts = [(key, table[key].get(index)) for key in sorted(table) if abs(key[1] - heading) < 1e-6]
    pts = [(wave_number(key[0], depth) * radius, value) for key, value in pts if value is not None]
    return [p[0] for p in pts], [p[1] for p in pts]


def plot(results, out, ref_rows):
    """PNGs of HAMS against the reference: |X_1|, |X_2| per cylinder against k a (the reference on a fine grid, HAMS on its
    meshes), the near-field drift per cylinder, the far-field drift of the array, and the largest errors against panels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import drift_convergence as dc
    os.makedirs(out, exist_ok=True)
    arr = layout("square")
    a, h = arr.radii[0], arr.depth
    folders = case_folders(results)
    hams = [(name.replace("l_square_diff_", "HAMS "), read_hams(folder)) for name, folder in folders]
    markers = ["o", "s", "^", "D"]
    omegas = np.linspace(0.36, 1.95, 160)
    written = []
    for heading in HEADINGS:
        sols = [Solution(arr, w, math.radians(heading), order=LOAD_ORDER) for w in omegas]
        ka = np.array([s.k * a for s in sols])
        X = np.array([[to_hams(exciting_loads(s, j)) for j in range(len(arr))] for s in sols])
        D = np.array([[near_field_drift(s, j) for j in range(len(arr))] for s in sols])
        F = np.array([far_field_drift(s) for s in sols])
        for kind in ("X", "near"):
            fig, axes = plt.subplots(2, len(arr), figsize=(15, 6.5), sharex=True)
            for j in range(len(arr)):
                for r, idx in enumerate((1, 2)):
                    ax = axes[r, j]
                    if kind == "X":
                        ax.plot(ka, np.abs(X[:, j, idx - 1]), "k-", lw=1.2, label="exact")
                        label = "|X_%d| / (rho g A) [m^2]" % idx
                    else:
                        ax.plot(ka, D[:, j, idx - 1], "k-", lw=1.2, label="exact")
                        ax.axhline(0.0, color="0.6", lw=0.6)
                        label = "mean F_%s / (rho g A^2) [m]" % "xy"[idx - 1]
                    for mi, (name, data) in enumerate(hams):
                        table = data[0] if kind == "X" else data[1]
                        xs, ys = _hams_points(table, 6 * j + idx, heading, h, a)
                        ax.plot(xs, [abs(v) for v in ys] if kind == "X" else ys, markers[mi % 4], ms=4, mfc="none",
                                label=name)
                    ax.set_title("cylinder %d at (%g, %g)" % (j + 1, arr.centres[j, 0], arr.centres[j, 1]), fontsize=9)
                    if r == 1:
                        ax.set_xlabel("k a")
                    if j == 0:
                        ax.set_ylabel(label)
                    ax.grid(alpha=0.3)
            axes[0, 0].legend(fontsize=7)
            what = "exciting force (ExcForce.3)" if kind == "X" else "near-field mean drift (MeanDriftPressure.9)"
            fig.suptitle("Four bottom-mounted cylinders, a = 5 m, h = 10 m, centres (+-10, +-10) m, heading %g deg: %s"
                         % (heading, what))
            fig.tight_layout()
            path = os.path.join(out, "le_%s_beta%d.png" % ("exciting_force" if kind == "X" else "drift_per_cylinder",
                                                            heading))
            fig.savefig(path, dpi=130)
            plt.close(fig)
            written.append(path)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        for r in range(2):
            ax = axes[r]
            ax.plot(ka, F[:, r], "k-", lw=1.2, label="exact far-field (Kochin)")
            ax.plot(ka, D[:, :, r].sum(axis=1), "k:", lw=1.0, label="exact near-field, summed")
            for mi, (name, data) in enumerate(hams):
                xs, ys = _hams_points(data[2], r + 1, heading, h, a)
                ax.plot(xs, ys, markers[mi % 4], ms=4, mfc="none", label=name + " far-field")
            ax.set_xlabel("k a")
            ax.set_ylabel("mean F_%s / (rho g A^2) [m]" % "xy"[r])
            ax.grid(alpha=0.3)
        axes[0].legend(fontsize=7)
        fig.suptitle("Four bottom-mounted cylinders, heading %g deg: mean drift of the array (MeanDriftMomentum.8)" % heading)
        fig.tight_layout()
        path = os.path.join(out, "le_drift_array_beta%d.png" % heading)
        fig.savefig(path, dpi=130)
        plt.close(fig)
        written.append(path)
    return written + plot_convergence(results, out, ref_rows)


def plot_convergence(results, out, ref_rows):
    """PNG of the largest error of each quantity of `compare` against the panel count. Far-field M_z is left out: it is
    zero by symmetry at both headings, so its value is round-off (under 3e-6 of |F| times 10 m), not a mesh error."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import drift_convergence as dc
    folders = case_folders(results)
    rows = []
    for name, folder in folders:
        rows += compare_case(folder, ref_rows, name)
    cases = {c.name: c for c in dc.linton_evans_cases()}
    named = [name for name, _ in folders if name in cases]
    if not (rows and named):
        return []
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for label, key in ERRORS:
        if key == "Mz_err":
            continue
        xs = [dc.panels_of(cases[name]) for name in named if worst(rows, name, key) is not None]
        ys = [worst(rows, name, key) for name in named if worst(rows, name, key) is not None]
        if xs:
            ax.loglog(xs, ys, "o-", label=label)
    x = np.array([dc.panels_of(cases[name]) for name in named], dtype=float)
    ax.loglog(x, 0.1 * (x / x[0]) ** -1.0, "k:", lw=1, label="slope of second order (error ~ 1/panels)")
    ax.set_xlabel("panels in all (four cylinders)")
    ax.set_ylabel("largest error over frequencies, headings, cylinders")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    ax.set_title("HAMS against the exact array solution")
    fig.tight_layout()
    path = os.path.join(out, "le_convergence.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return [path]


def plot_reference(out):
    """PNGs of the reference alone: |eta| around the square at heading 45 deg, at 1.10 and 1.70 rad/s, with the field points;
    and |X_1| and the near-field mean F_x per cylinder against k a at both headings."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(out, exist_ok=True)
    arr = layout("square")
    written = []
    xs = np.linspace(-40, 40, 201)
    gx, gy = np.meshgrid(xs, xs)
    inside = np.zeros(gx.shape, dtype=bool)
    for (cx, cy), a in zip(arr.centres, arr.radii):
        inside |= np.hypot(gx - cx, gy - cy) < a
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, omega in zip(axes, (1.10, 1.70)):
        sol = Solution(arr, omega, math.radians(45.0))
        eta = np.abs(sol.field(gx.ravel(), gy.ravel())).reshape(gx.shape)
        eta[inside] = np.nan
        im = ax.pcolormesh(gx, gy, eta, shading="auto", cmap="viridis")
        fig.colorbar(im, ax=ax, label="|eta| / A")
        pts = np.array(FIELD_POINTS)
        ax.plot(pts[:, 0], pts[:, 1], "w+", ms=7)
        ax.set_aspect("equal")
        ax.set_title("omega = %.2f rad/s, k a = %.2f, heading 45 deg" % (omega, sol.k * arr.radii[0]))
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    fig.suptitle("Exact wave field around the four cylinders (+: the HAMS field points)")
    fig.tight_layout()
    path = os.path.join(out, "le_reference_field.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path)
    omegas = np.linspace(0.36, 1.95, 160)
    single = layout("single")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for c, heading in enumerate(HEADINGS):
        sols = [Solution(arr, w, math.radians(heading), order=LOAD_ORDER) for w in omegas]
        ones = [Solution(single, w, math.radians(heading), order=LOAD_ORDER) for w in omegas]
        ka = np.array([s.k * arr.radii[0] for s in sols])
        for j in range(len(arr)):
            axes[0, c].plot(ka, [abs(exciting_loads(s, j)[0]) for s in sols], label="cylinder %d" % (j + 1))
            axes[1, c].plot(ka, [near_field_drift(s, j)[0] for s in sols], label="cylinder %d" % (j + 1))
        axes[0, c].plot(ka, [abs(exciting_loads(s, 0)[0]) for s in ones], "k--", lw=1, label="one cylinder alone")
        axes[1, c].plot(ka, [near_field_drift(s, 0)[0] for s in ones], "k--", lw=1, label="one cylinder alone")
        axes[0, c].set_title("heading %g deg" % heading)
        axes[1, c].set_xlabel("k a")
        for r in range(2):
            axes[r, c].grid(alpha=0.3)
    axes[0, 0].set_ylabel("|X_1| / (rho g A) [m^2]")
    axes[1, 0].set_ylabel("mean F_x / (rho g A^2) [m]")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Exact loads per cylinder in the square array, against one cylinder alone (MacCamy-Fuchs)")
    fig.tight_layout()
    path = os.path.join(out, "le_reference_loads.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path)
    fx, fy = kim_fig3()
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    n = np.arange(1, 20)
    for ax, ref, pub, label in ((axes[0], fx, KIM_FIG3_X, "|F_x|"), (axes[1], fy, KIM_FIG3_Y, "|F_y|")):
        ax.plot(n, pub, "o", mfc="none", ms=9, label="Kim (1993), Fig. 3, digitised")
        ax.plot(n, ref, "k+", ms=9, label="this reference")
        ax.set_ylabel(label + " / (pi rho g a h A)")
        ax.grid(alpha=0.3)
    axes[1].set_xlabel("cylinder number")
    axes[0].legend(fontsize=8)
    fig.suptitle("19 bottom-mounted cylinders in a row, spacing 4a, h = 12a, k a = 0.5, heading 45 deg\n"
                 "largest difference %.4f (a marker is about 0.0015 high)"
                 % max(np.max(np.abs(fx - KIM_FIG3_X)), np.max(np.abs(fy - KIM_FIG3_Y))), fontsize=10)
    fig.tight_layout()
    path = os.path.join(out, "le_reference_kim1993_fig3.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path)
    return written



def plot_linton_evans(results, out):
    """PNG of Linton & Evans (1990), Fig. 3 against the reference and HAMS: the total mean drift along the wave over that on
    one cylinder alone, f/f_1, for heading 45 deg (the HAMS case is their layout scaled by 5). The published values are the
    digitised LE_FIG3_TOTAL; HAMS is MeanDriftMomentum.8 (far-field) of the l_square_* cases under RESULTS, if given, with
    f_1 from the reference. A second panel gives |X . e_beta|/|F| of the two cylinders on the diagonal, reference and HAMS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(out, exist_ok=True)
    arr = layout("square")
    a, depth = arr.radii[0], arr.depth
    beta = math.pi / 4.0
    e = np.array([math.cos(beta), math.sin(beta)])
    one = Array([(0.0, 0.0)], a, depth)
    kas = np.linspace(0.15, 1.95, 241)
    total, diag = [], []
    for ka in kas:
        k = ka / a
        omega = math.sqrt(G * k * math.tanh(k * depth))
        s1 = Solution(one, omega, beta, order=16)
        f1 = float(np.dot(near_field_drift(s1, 0)[:2], e))
        F = abs(np.dot(exciting_loads(s1, 0)[:2], e))
        sol = Solution(arr, omega, beta, order=16)
        total.append(sum(float(np.dot(near_field_drift(sol, j)[:2], e)) for j in range(len(arr))) / f1)
        # this module's cylinders 1 and 3 lie on the diagonal along the wave: Linton & Evans's cylinders 2 and 4
        diag.append([abs(np.dot(exciting_loads(sol, j)[:2], e)) / F for j in (0, 2)])
    diag = np.array(diag)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    ax.plot(kas, total, "k-", lw=1.2, label="exact (this reference)")
    ax.plot(LE_FIG3_KA, LE_FIG3_TOTAL, "o", ms=4, mfc="none", color="tab:red",
            label="data from Linton & Evans (1990), Fig. 3, digitised")
    ax.axvspan(*LE_FIG3_BAND, color="0.92", zorder=0)
    ax2 = axes[1]
    ax2.plot(kas, diag[:, 0], "k-", lw=1.2, label="exact, lee cylinder (L&E 2)")
    ax2.plot(kas, diag[:, 1], "k--", lw=1.2, label="exact, weather cylinder (L&E 4)")
    markers = ["s", "^", "D", "v"]
    for mi, (name, folder) in enumerate(case_folders(results) if results else []):
        X, _, far, _ = read_hams(folder)
        xs, ys, xd, yd1, yd3 = [], [], [], [], []
        for key in sorted(far):
            if abs(key[1] - 45.0) > 1e-6:
                continue
            k = wave_number(key[0], depth)
            s1 = Solution(one, key[0], beta, order=16)
            f1 = float(np.dot(near_field_drift(s1, 0)[:2], e))
            F = abs(np.dot(exciting_loads(s1, 0)[:2], e))
            xs.append(k * a)
            ys.append((far[key][1] * e[0] + far[key][2] * e[1]) / f1)
            if key in X:
                # ExcForce.3 holds conj(X): the modulus of the component along e is unchanged
                xd.append(k * a)
                yd1.append(abs(X[key][1] * e[0] + X[key][2] * e[1]) / F)
                yd3.append(abs(X[key][13] * e[0] + X[key][14] * e[1]) / F)
        label = name.replace("l_square_diff_", "HAMS ")
        ax.plot(xs, ys, markers[mi % 4], ms=5, color="C%d" % mi, label=label + " (far-field)")
        ax2.plot(xd, yd1, markers[mi % 4], ms=5, color="C%d" % mi, label=label)
        ax2.plot(xd, yd3, markers[mi % 4], ms=5, mfc="none", color="C%d" % mi)
    ax.set_xlabel("k a")
    ax.set_ylabel("f / f_1, total mean drift along the wave")
    ax.set_title("total drift (grey: near-trapping band)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    ax2.set_xlabel("k a")
    ax2.set_ylabel("|X . e_beta| / |F|")
    ax2.set_title("first-order force on the diagonal cylinders (HAMS: filled lee, open weather)", fontsize=10)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=7)
    fig.suptitle("Four cylinders, a/h = 1/2 (h half the side), heading 45 deg: data from Linton, C.M. & Evans, D.V. (1990), "
                 "J. Fluid Mech. 215, 549-569, Fig. 3", fontsize=10)
    fig.tight_layout()
    path = os.path.join(out, "le_linton_evans_fig3_hams.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return [path]

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command")
    p = sub.add_parser("check", help="the self-checks (CHECKS); exit 1 if one fails")
    p.add_argument("names", nargs="*")
    p = sub.add_parser("reference", help="write the reference rows of the HAMS case to a CSV")
    p.add_argument("--out", required=True)
    p = sub.add_parser("compare", help="errors of the l_square_* cases under RESULTS")
    p.add_argument("results")
    p.add_argument("--out")
    p.add_argument("--reference", help="a CSV of `reference`, instead of computing the rows")
    p = sub.add_parser("plot", help="PNGs of the reference, and of HAMS against it when RESULTS is given")
    p.add_argument("results", nargs="?")
    p.add_argument("--out", required=True)
    p.add_argument("--reference")
    args = parser.parse_args(argv)
    if args.command == "check":
        return 1 if run_checks(args.names) else 0
    if args.command == "reference":
        rows = reference_rows()
        write_rows(args.out, rows)
        print("%d rows in %s" % (len(rows), args.out))
        return 0
    if args.command == "compare":
        ref = read_rows(args.reference) if args.reference else reference_rows()
        rows = []
        for name, folder in case_folders(args.results):
            rows += compare_case(folder, ref, name)
        summary(rows)
        if args.out:
            write_rows(args.out, rows, COMPARE_COLUMNS)
            print("%d rows in %s" % (len(rows), args.out))
        return 0
    if args.command == "plot":
        written = plot_reference(args.out)
        written += plot_linton_evans(args.results, args.out)
        if args.results:
            ref = read_rows(args.reference) if args.reference else reference_rows()
            written += plot(args.results, args.out, ref)
        print("\n".join(written))
        return 0
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
