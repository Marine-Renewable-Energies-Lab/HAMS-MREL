"""Semi-analytical reference for a truncated vertical cylinder, by eigenfunction matching, to validate HAMS-MREL's mean drift.

Problem (design: docs/superpowers/specs/2026-11-06-truncated-cylinder-drift-reference-design.md, sections 2 and 3). A vertical
cylinder of radius a and draft d, sharp bilge, in water of depth h, with the gap g = h - d below it. Linear potential flow,
time factor exp(-i omega t), nu = omega**2/G. G and RHO are make_mesh.G and make_mesh.RHO, the values of src/WavDynMods.f90.
z is up, z = 0 is the still waterline, z = -h the sea bed. Each azimuthal order m is solved on its own.

Vertical modes. Exterior, r > a, -h < z < 0:

    Z_0(z) = cosh k (z + h)/cosh k h,       k tanh(k h) = nu, k > 0,
    Z_n(z) = cos k_n (z + h)/cos k_n h,     k_n tan(k_n h) = -nu, k_n h in ((n - 1/2) pi, n pi), n = 1 ... N,

so that Z_n(0) = 1. Interior, r < a, -h < z < -d: cos(lambda_n (z + h)), lambda_n = n pi/g, n = 0 ... N, not normalised.
k is drift_convergence.wave_number (Newton). k_n is solved as y = n pi - atan(nu h/y), y = k_n h, by brentq for
n <= BRENTQ_ROOTS; this form has no pole in the bracket (as evanescent_roots in annular_lid_reference.py on
feature/damping-lid). Beyond, the same equation by fixed-point iteration, vectorised (contraction below 1.3e-5 there).

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

Free body and motions. body_matrices applies the formulas of make_mesh.hydrostatics, the rule that wrote the study's
Hydrostatic.in (a uniform solid prism at its draft, z_G = z_B = -d/2, mass and restoring about the waterplane centre, here
the origin), to the exact circle. HAMS reads MATX and CRS from Hydrostatic.in (ReadHydroStatic in src/HydroStatic.f90;
BLNR, BQDR and KSTF are zero in the study) and SolveMotion (src/SolveMotion.f90, called from src/HAMS_Prog.f90 after EFORCE)
solves [-omega**2 (MATX + AMAS) + i omega (BDMP + BLNR) + CRS + KSTF] DSPL = -Im EXFC - i Re EXFC = -i conj(EXFC), with AMAS
and BDMP dimensional and EXFC = X = i omega rho int phi n dS, AMP = 1 (Initialisation in src/HydroStatic.f90). That is the
conjugate of [-omega**2 (M + A) - i omega B + C] xi = X, the equation of motion in exp(-i omega t), times i: DSPL =
conj(i xi) = conj(xi/eta_I(0)), the convention of ExcForce.3. It writes to unit 63, Motion.4 (opened in src/InputFiles.f90),
MOD = abs(DSPL)/NFAC, NREL + i NIMG = DSPL/NFAC, NFAC = AMP/REFL**MEXP = 1 (REFL = 1). (PrintOutput.f90 writes only the
headers of the .rao files; the rows come from SolveMotion.) motions returns xi itself, per unit wave amplitude, rotations in
rad/m, in this module's convention, as total_solution and the near-field terms need it; hams_motion converts to Motion.4's
complex value, and abs(xi) is Motion.4's MOD. total_solution adds sum_j (-i omega xi_j) phi_j, phi_j per unit velocity.

Far-field drift. far_field_fx is Maruo's formula, F_x/(rho g A**2) = (1/2)(C_g/C_p) int abs(a)**2 (1 - cos theta) dtheta,
heading 0, with a(theta) the outgoing elevation amplitude of far_field_amplitude (spec section 5.3 of
docs/superpowers/specs/2026-11-04-second-order-waves-design.md, "Far-field."; there a(theta) is the same amplitude, in
exp(-i omega t)). That spec's general form, (1/2)(C_g/C_p)[X - int abs(a)**2 cos theta dtheta] with the optical term
X = -sqrt(8 pi/k) Re(e^(3 i pi/4) a(0)) for this module's a and eta_I = -i e^(i k x), is the form HAMS evaluates
(pair_core in src/FarFieldDrift.f90, with its own Kochin amplitude); the two are equal when no energy is lost,
X = int abs(a)**2 dtheta, which scattering_energy checks for the fixed body and free_energy for the free one. HAMS writes
it to MeanDriftMomentum.8 (write_far_rows in src/SecondOrderDriver.f90, write_drift_row and output_value in
src/SecondOrderOutput.f90) as conj(F)/(rho g REFL), REFL = 1, so the Re column is F_x/(rho g A**2) itself (F_x is real),
positive along +x, the direction the wave travels at heading 0: F e_x with e_x the unit vector of +x.

Gap basis (three families). The single family above is rho**(-1/3) times a power series in rho at the corner. The local
expansion of the gap velocity in the 270 degree wedge has every exponent k/3 - 1/3: rho**(-1/3), rho**0 (the rigid-body
radial velocity of a moving wall), rho**(1/3) (the second wedge eigenfunction), rho**(2/3), ... The near field squares the
velocity next to the corner and sees the missing ones: with the single family the free body's near-field F_x settles
(in N) 17 % (P = 20) and 4 % (P = 60) from the far field at 0.70 rad/s. So p_max may be a tuple (basis_spec): the
families (1 - t**2)**(lam - 1/2) C_2p^lam(t), lam = 1/6, 1/2, 5/6 (FAMILIES; exponents -1/3 + n, n, 1/3 + n), each even
in t, with their own p_max; an int P is the single family of Tasks 1 to 4, unchanged. Each family's inner products are
Gegenbauer's integral with its lam (c_p = pi 2**(1 - lam) Gamma(2p + 2 lam)/((2p)! Gamma(lam))), and check_inner_products
tests all three against quad. The families are not orthogonal and together nearly dependent (the Gram matrix
int u_p u_q dt of (20, 4, 4) has 8 eigenvalues below 1e-10 of the largest; _gram, exact in sigma, 1 - abs(t) = sigma**3),
so _solve_order solves the Galerkin system in the orthonormal combinations of the eigenvectors above GRAM_CUTOFF = 1e-10
(_gram_transform), gap = T c; for m = 0 the zero-flux condition borders the system, since each family's p = 0 carries
flux. Without that reduction the coefficients reach 3e4 and the near field is erratic in N and P.

Matching tail (_matching_tail, MATRIX_TAIL). The matrix sums sum_n E E/(N kappa) - F F/(nu mu) of the corner family
converge only like N**(-4/3) (E_pn ~ n**(-2/3)), and the wall projection like N**(-5/3); the Hankel expansion gives
their non-oscillating leading tails in closed form (Hurwitz zeta), checked against brute-force sums from N = 2000 to 2e5
to 1e-3 of the tail. Adding them makes the near-field identity converge like N**(-2) instead of N**(-1.3) and leaves the
single-family results N-converged at N = 4000 (fixed body, 0.70 rad/s, 6.21e-4 at N = 4000 and 6.25e-4 at N = 32000).
It applies to every solve, single family included.

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

Near field. near_field_terms writes the terms of spec 5.3 of docs/superpowers/specs/2026-11-04-second-order-waves-design.md
(WAMIT eq. 15.54-15.57) under the kernel T- of spec 5.4, P(f, g) = (f_i conj(g_j) + conj(f_j) g_i)/4, over rho g, n the
unit normal out of the fluid (into the body), x relative to xr, xi' = xi + alpha x xr the translation of xr and alpha =
(0, xi_5, 0) the rotation (the module's xi are about the origin), X = xi' + alpha x (x - xr), phi_t = -i omega phi:

    I    (1/2) int_WL r' [P(zeta, zeta) - 2 P(zeta, X_3)] a dtheta, zeta = (i omega/g) phi(a, theta, 0), r' = n' or
         (x - xr) x n', n' = (-cos theta, -sin theta, 0) the unit horizontal normal into the body,
    Ia   zero: the wall is vertical at the waterline (n_z = 0),
    II   -(1/g) int_S r (1/2) sum_k P(d_k phi, d_k phi) dS,  r = n or (x - xr) x n,
    III  -(1/g) int_S r sum_k P(X_k, d_k phi_t) dS,
    IV   force P(alpha x G), moment P(alpha x H) + P(xi' x G), G = -(1/g) int n phi_t dS, H = -(1/g) int (x - xr) x n phi_t dS,
         the first-order dynamic pressure only,
    V    F_S,z and M_S,y of Task 7 of docs/superpowers/plans/2026-11-04-mean-wave-drift.md, with A = pi a**2, (x_f, y_f) =
         -(x_r, y_r), L_ij about xr, V = pi a**2 d, (x_b, y_b) = -(x_r, y_r), Z_o = z_r (_hydrostatic_term).

The order is that of MeanDriftTerms.txt (T_WATERLINE ... T_HYDROSTATIC in src/NearFieldDrift.f90); near_field's total is in
MeanDriftPressure.9's normalisation, conj(F)/(rho g REFL**k), REFL = 1 (write_drift_row in src/SecondOrderOutput.f90), real
for i = j. F_S,z depends on Z_o (the second-order position of a body rotated about xr), so F_z is for xr = 0, as HAMS's XR.

Velocity on the body: the normal component is the body condition, d phi/dn = -i omega X.n; the tangential ones and phi come
from the series. The raw exterior series of d phi/dz on the wall converges like n**(-2/3) e^(i n pi rho/h) near the corner
(rho the distance from it), 1 to 5 % pointwise at N = 16000 within a metre of it, so _body_series (i) extends the series
beyond N: the coefficients of every mode follow in closed form from the gap velocity sum_p gap[p] u_p and the wall
velocity (Order.wall_velocity), exactly as _solve_order computes the first N; (ii) sums the wall series with a smooth
window (_smooth_window), which removes the partial-sum oscillation; (iii) adds, per node, the exact sum less the windowed
sum of the asymptotic tails, known in closed form from the Hankel expansion of E_pn: the corner term
(g/h) sqrt(2/pi) S (n pi g/h)**(-2/3) cos(k_n g - pi/3) sin(k_n s), S = sum_p gap[p] c_p, the jump of the wall velocity
-V_c [cos(k_n rho) - cos(k_n (s + g))]/(n pi), and the n**(-5/3) term of phi, by polylogarithms (_polylog, DLMF
25.12.12). The bottom (interior) series converges like exp(-lambda_n (a - r)); its corner tail
(1/2) sqrt(2/pi) S (n pi)**(-2/3) sqrt(a/r) exp(-lambda_n (a - r)) is added the same way. With the three-family basis
(Gap basis) each family lam has such tails with beta = lam + 1/2 and phase chi = (lam/2 + 1/4) pi (_corner_models); the
lam = 1/2 family's tail is that of a jump, so its corner value sum_p gap[p] is subtracted from the wall velocity V_c in
the jump term, and the lam = 5/6 tails (n**(-4/3)) are added as the lam = 1/6 ones. Measured at 0.70 rad/s, order 1,
N = 2000: rho**(1/3) d phi/dz on the wall at rho = 1e-3 ... 0.1 m changes by below 3e-5 between 16384 and 65536 modes, and
rho**(1/3) times the tangential velocity on the wall and on the bottom tend to the same limit at the corner (0.765 and
0.766, fixed body), as the 270 degree wedge requires.

Quadrature (body_quadrature). Near the corner abs(grad phi)**2 is rho**(-2/3) times a series in rho**(1/3) (the wedge
exponents 2k/3 and the regular terms), not in rho. Gauss-Jacobi for the weight rho**(-2/3) integrates only the first term
exactly and converges like n**(-4/3) for the rest: F_x of the fixed body at 0.70 rad/s, N = 4000, P = 20, changes by 6e-4
of F_x between 8 and 32 nodes and is not settled at 64. With rho = L s**3 every term is a polynomial in s, so the default
CORNER_RULE "cube" uses Gauss-Legendre in s = (rho/L)**(1/3) on the corner zone (min(g, d) on the wall, min(g, a/2) on the
bottom), which carries the weight rho**(-2/3) exactly in its Jacobian 3 L s**2: the same case agrees to 1e-7 between 8 and
32 nodes. Gauss-Legendre in r beyond. The trapezoidal rule in theta is exact.

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
- bottom_limit_x1: gap 1e-3 h (BOTTOM_GAP), LIGHT_RESOLUTION, every WIDE_FREQUENCIES: X_1 against MacCamy-Fuchs,
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
- scattering_energy: fixed body, default cylinder, LIGHT_RESOLUTION, every WIDE_FREQUENCIES: the scattered energy flux
  2 pi sum_m abs(a_m)**2 equals the optical-theorem term -sqrt(8 pi/k) Re(e^(3 i pi/4) a(0)), a(0) = sum_m a_m, from
  far_field_amplitude, to 1e-6. Energy is conserved order by order, so it tests each order's propagating coefficient and the
  far-field normalisation, not the convergence in M, N or P.
- reciprocity: default cylinder, every WIDE_FREQUENCIES: A15 = A51 and B15 = B51 to 1e-8, and X_j from diffraction equal
  to haskind_excitation of the mode j radiation solution to 1e-6, j = 1, 3, 5. Values at N = CHECK_RESOLUTION's 8000 and
  2 N = 16000, P = 20, extrapolated in N with the exponent 5/3 (Richardson); the message also prints the residuals at
  N = 16000 itself (about 1e-6). With MATRIX_TAIL no identity holds to round-off any more (the solve sees the tail of the
  sums, the forces only the N modes); all converge in N, and the extrapolated residuals are 2e-9. check_reciprocity's docstring says what each part tests.
- radiation_energy: same resolution and extrapolation: B11, B33 and B55 from the body pressure equal
  radiated_energy_damping, the energy flux at infinity, to 1e-6. It tests the flux formula and the sign and normalisation of
  B, and for heave and pitch the particular solutions. Neither radiation check bounds the error of A, B or X themselves.
- particular_bottom: _Particular.bottom, closed form, equals quad of fields(r, -d) r**(m+1) over 0 < r < a to 1e-12, heave
  and pitch, gaps 90 m and 0.1 m. The bottom term is real and enters only Re A33 and Re A55 (and so the heave and pitch
  RAOs), which neither radiation check sees. It does not test fields; the radiation checks do, through the gap.
- bottom_limit_drift: gap 1e-3 h, fixed, LIGHT_RESOLUTION, every WIDE_FREQUENCIES: far_field_fx against
  drift_convergence.maccamy_fuchs_fx to 1e-3, after the same extrapolation in the gap as bottom_limit_x1, 2 F(g) - F(2g)
  (the gap itself changes F_x by 8.7e-4 at 0.40 rad/s, linear in g, so the raw value is mostly the gap's physical
  effect, not an error); the message prints the raw difference at g = 1e-3 h too. Like bottom_limit_x1 it asserts
  abs(F(2g) - MF)/abs(F(g) - MF) in [1.8, 2.2] (measured 1.95) where abs(F(g) - MF) is above 1e-13 abs(MF); from about
  1.2 rad/s the gap's effect is round-off and the ratio is not tested there. The extrapolated value is blind to errors of
  the interior and gap terms, which enter F_x at O(g), and the extrapolation removes O(g); the ratio and the raw value see
  them. It covers far_field_amplitude summed over orders 0 ... 15, Maruo's formula, C_g/C_p and the normalisation over
  rho g A**2. It does not cover the radiated waves (free body), nor a deep gap; free_energy covers the free body's
  outgoing wave, not its drift.
- body_matrices: body_matrices(cyl) against make_mesh.hydrostatics on a 4096-sided polygon of the same radius, to 1e-5
  relative in every nonzero entry, with the same zero pattern (the polygon's area error is about 4e-7). It covers the
  transcription of the rule, not the rule itself.
- free_energy: free body, default cylinder, LIGHT_RESOLUTION, every WIDE_FREQUENCIES: for the total outgoing wave of
  total_solution with the motions, the optical-theorem term equals the energy flux 2 pi sum_m abs(a_m)**2 to 1e-10 of
  the fixed body's flux (the free body's own flux is 1.5e-5 at 0.40 rad/s, where it nearly follows the wave, and the
  difference of two O(1) sums has round-off 1e-16 against it). Measured to round-off at N = 1000 and 4000: a free body
  with real M, A and C does no mean work, so the outgoing wave carries the energy the body takes out of the incident
  wave. It tests the sign and phase of X against the radiation
  waves, the normalisation of B, the velocity factor -i omega and the combination in total_solution (a sign flip of xi
  gives 0.6, conj(xi) 0.16, 1.01 xi 7e-3). It cannot see an error in A, M or C, which change xi but keep the body
  lossless; body_matrices and particular_bottom cover M, C and the bottom term of A, and the rest of A only the
  near-field against far-field F_x of the free body (fx_identity) tests.
- long_wave_motions: free body at 0.05 rad/s (kh = 0.16), LIGHT_RESOLUTION: xi_3/eta_I(0) -> 1 and xi_5/(-d eta_I/dx(0))
  -> 1 (heave follows the elevation, pitch the slope; a positive xi_5 lowers the +x side), each within
  omega**2 (M_jj + A_jj)/C_jj (long_wave_tolerance: 4.6e-3 heave, 7.5e-3 pitch), the size of the O(omega**2) inertia
  correction; measured 2.1e-4 and 2.4e-3, falling like omega**2 in heave (0.1, 0.05, 0.025 rad/s: 8.8e-4, 2.1e-4,
  5.2e-5). It covers the real part of the RAO solve: C, the zero-frequency X and its phase, and the M15 + A15 coupling,
  which enters pitch at leading order through surge, xi_1 ~ X_1/omega**2 (a sign flip of M15 or A15 gives 4.0 or 1.7,
  dropping A15 0.86); a sign flip of the inertia gives 8.9e-3 in heave, C55 1 % high 7.6e-3 in pitch. It does not see
  B (1e-7 at this omega) nor an error of A33 or A55 below the tolerance.
- fx_identity: fixed and free, every WIDE_FREQUENCIES, NEAR_CHECK_RESOLUTION (gap basis (20, 4, 4), N = 64000):
  near_field's F_x against far_field_fx, to 5e-4 of abs(far F_x), with the floor of drift_convergence.percent (1e-4 of
  the motion's largest abs(far F_x)), which acts only at free 0.40 rad/s (far F_x = 3e-6). Measured: fixed below 2.4e-6,
  free 4.0e-4 (0.70), 4.1e-4 (0.80, near the pitch resonance), 1.9e-4 (0.40) and below 1.7e-4 elsewhere. The gate of
  the near field: it tests every term of F_x, the body velocities (the body condition and the corner-corrected series),
  the quadrature, and through the free body the motions, A and the radiation potentials near the body. Prescribed single
  motions at 0.70 rad/s, (20, 4, 4), N = 64000, against the momentum far field (X - int abs(a)**2 cos(theta) dtheta,
  valid without the lossless assumption): surge 1.8e-6, heave 5.6e-8, pitch -1.6e-4 absolute (F 1.92, 0.15, 0.80); the
  pitch residual falls slowly with the basis ((10, 1, 1): -2.0e-4) and is most of the free-body error at 0.70 rad/s.
  It does not test F_z or M_y; reference_point tests M_y's moment arms.
- reference_point: xr = (0, 0, -3) against the origin, fixed and free, 0.70 and 1.40 rad/s, LIGHT_RESOLUTION:
  M_y(xr) = M_y(0) - (xr x F)_y and F_x unchanged, to 1e-10. The motions are re-expressed about xr (xi' = xi + alpha x xr),
  so it tests that every moment arm, the IV moment (alpha x H + xi' x G) and M_S transform as a rigid shift. F_z is not
  compared: F_S,z depends on Z_o (Near field).
- pair_symmetry: free body, omega_i = 0.70, omega_j = 0.95: T-_ij = conj(T-_ji) for every term and component, to 1e-12. It
  tests that each state's frequency, motion and conjugation enter the kernel symmetrically, not the values.
"""
import argparse
import functools
import math
import os
import sys
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer, gamma, gammaln, h1vp, hankel1, ive, jv, jvp, kve, roots_jacobi, zeta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import drift_convergence  # noqa: E402
import make_mesh  # noqa: E402

LAMBDA = 1.0 / 6.0     # Gegenbauer parameter of the corner basis
ALPHA = LAMBDA - 0.5   # Jacobi parameter alpha = beta = -1/3 of its weight
GAUSS_JACOBI_NODES = 96   # default size of the Gauss-Jacobi rule: exact for polynomials of degree < 192

# Resolution of the root and inner product checks
CHECK_P = 10
CHECK_N = 20
# Resolution of the radiation checks: N exterior and interior modes, corner basis size P, orders 0 ... M. N = 400 P, above the
# N >= 100 P rule (Convergence); the radiation checks extrapolate from N and 2 N (RICHARDSON_EXPONENT).
CHECK_RESOLUTION = dict(n_modes=8000, p_max=20, orders=15)
RICHARDSON_EXPONENT = 5.0 / 3.0   # the radiation identities' residual falls as N**(-5/3); see check_reciprocity
# Resolution of the checks that need no convergence in N: scattering_energy (an identity that holds order by order at any
# truncation) and the bottom-limit checks (1e-3 against MacCamy-Fuchs; their extrapolated error is the O(g**2) remainder,
# the same at N = 500 and 8000). N = 100 P, the rule of Convergence.
LIGHT_RESOLUTION = dict(n_modes=1000, p_max=10, orders=15)


@dataclass
class Cylinder:
    """Truncated vertical cylinder, lengths in m. Defaults: the study's sharp hull (HULLS["sharp"] in drift_convergence.py)."""
    radius: float = 10.0
    draft: float = 10.0
    depth: float = 100.0

    @property
    def gap(self):
        return self.depth - self.draft


BRENTQ_ROOTS = 500   # wave_numbers: roots by brentq; beyond, by fixed-point iteration


def wave_numbers(omega, depth, count):
    """(k, kn): k > 0 with k tanh(k h) = nu, and kn, the first count roots of kn tan(kn h) = -nu, increasing."""
    nu = omega ** 2 / make_mesh.G
    k = drift_convergence.wave_number(omega, depth)
    kn = np.empty(count)
    for n in range(1, min(count, BRENTQ_ROOTS) + 1):
        y = brentq(lambda y: y - n * math.pi + math.atan(nu * depth / y), (n - 0.5) * math.pi, n * math.pi,
                   xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=200)
        kn[n - 1] = y / depth
    if count > BRENTQ_ROOTS:
        # y = n pi - atan(nu h/y) by fixed-point iteration: its derivative nu h/(y**2 + (nu h)**2) is below 1e-2 here
        n = np.arange(BRENTQ_ROOTS + 1, count + 1)
        y = (n - 0.5) * math.pi
        for _ in range(16):
            y = n * math.pi - np.arctan(nu * depth / y)
        kn[BRENTQ_ROOTS:] = y / depth
    return k, kn


# Gegenbauer parameters of the gap basis families: (1 - t**2)**(lambda - 1/2) C_2p^lambda(t) behaves like
# rho**(lambda - 1/2) times a power series in rho at the corner, so together they carry every corner exponent k/3 - 1/3
FAMILIES = (1.0 / 6.0, 0.5, 5.0 / 6.0)


@dataclass
class Family:
    """One family of the gap basis, u_p(t) = (1 - t**2)**(lam - 1/2) C_2p^lam(t), p = 0 ... p_max, with a Gauss-Jacobi rule
    for its weight: int_-1^1 (1 - t**2)**(lam - 1/2) f(t) dt ~ sum(weights * f(nodes)); poly[p, j] = C_2p^lam(nodes[j])."""
    lam: float
    p_max: int
    nodes: np.ndarray
    weights: np.ndarray
    poly: np.ndarray


@dataclass
class Basis:
    """The gap basis: the families of FAMILIES with p_max >= 0, in that order; the coefficients of an Order's gap are
    stacked family by family (slices). p_max is the spec it was built from (an int: the lambda = 1/6 family only)."""
    p_max: object
    families: list

    @property
    def size(self):
        return sum(f.p_max + 1 for f in self.families)

    @property
    def slices(self):
        out, start = [], 0
        for f in self.families:
            out.append((f, slice(start, start + f.p_max + 1)))
            start += f.p_max + 1
        return out

    def values(self, t):
        """u_p(t) of every family, stacked, shape (size,) + shape(t), for abs(t) < 1."""
        t = np.asarray(t, dtype=float)
        return np.concatenate([np.array([eval_gegenbauer(2 * p, f.lam, t) for p in range(f.p_max + 1)])
                               * (1.0 - t ** 2) ** (f.lam - 0.5) for f in self.families])

    def project(self, cyl, func):
        """int over the gap of u_p(t(z)) func(s) dz, s = z + h = g abs(t), stacked, by each family's Gauss-Jacobi rule;
        func maps an array of s to an array (..., len(s))."""
        g = cyl.gap
        return np.concatenate([0.5 * g * (f.poly * f.weights) @ np.asarray(func(g * np.abs(f.nodes))).T
                               for f in self.families])


GRAM_CUTOFF = 1e-10   # eigenvalues of the Gram matrix of a multi-family basis below this times the largest are dropped
GRAM_NODES = 200      # Gauss-Legendre nodes in sigma, 1 - abs(t) = sigma**3, for the Gram matrix


def _gram(basis):
    """G[p, q] = int_-1^1 u_p u_q dt for the stacked basis. With 1 - abs(t) = sigma**3 every product u_p u_q dt is a
    polynomial in sigma times a smooth factor (the exponents 3 lam - 3/2 are integers), so Gauss-Legendre in sigma."""
    sigma, w = _legendre_rule(0.0, 1.0, GRAM_NODES)
    t = 1.0 - sigma ** 3
    one_minus_t2 = sigma ** 3 * (2.0 - sigma ** 3)                 # 1 - t**2 without cancellation
    u = np.concatenate([np.array([eval_gegenbauer(2 * p, f.lam, t) for p in range(f.p_max + 1)])
                        * one_minus_t2 ** (f.lam - 0.5) for f in basis.families])
    return 2.0 * (u * (3.0 * sigma ** 2 * w)) @ u.T             # even: twice the integral over 0 < t < 1


def _gram_transform(basis, cutoff=None):
    """T, size x K: the columns span the basis with orthonormal combinations, T^T G T = I, keeping the eigenvalues of G
    above cutoff times the largest; and the kept eigenvalue range (for the report)."""
    cutoff = GRAM_CUTOFF if cutoff is None else cutoff
    values, vectors = np.linalg.eigh(_gram(basis))
    keep = values > cutoff * values[-1]
    return vectors[:, keep] / np.sqrt(values[keep]), (values[keep][0] / values[-1], int(np.sum(~keep)))


def basis_spec(p_max):
    """The p_max of each family in FAMILIES: an int P is (P, -1, -1), the single family of Tasks 1 to 4; a tuple gives the
    three (-1: the family is absent)."""
    if isinstance(p_max, (int, np.integer)):
        return (int(p_max), -1, -1)
    spec = tuple(int(p) for p in p_max)
    if len(spec) != len(FAMILIES) or spec[0] < 0:
        raise ValueError("p_max: an int or a tuple of %d ints with the first >= 0, got %r" % (len(FAMILIES), p_max))
    return spec


def gegenbauer_basis(p_max, nodes=GAUSS_JACOBI_NODES):
    """The gap basis for p_max (basis_spec), each family with a Gauss-Jacobi rule of `nodes` points for
    alpha = beta = lam - 1/2."""
    families = []
    for lam, p in zip(FAMILIES, basis_spec(p_max)):
        if p < 0:
            continue
        t, w = roots_jacobi(nodes, lam - 0.5, lam - 0.5)
        poly = np.array([eval_gegenbauer(2 * q, lam, t) for q in range(p + 1)])
        families.append(Family(lam, p, t, w, poly))
    return Basis(p_max, families)


def _gegenbauer_constant(p, lam=LAMBDA):
    """c_p = pi 2**(1 - lam) Gamma(2p + 2 lam)/((2p)! Gamma(lam)), elementwise in p (lam = 1/6: pi 2**(5/6) ...)."""
    p = np.asarray(p, dtype=float)
    return math.pi * 2.0 ** (1.0 - lam) * np.exp(gammaln(2 * p + 2 * lam) - gammaln(2 * p + 1) - gammaln(lam))


def _cos_integral(p_max, a, lam=LAMBDA):
    """M[p, j] = int_-1^1 (1 - t**2)**(lam - 1/2) C_2p^lam(t) cos(a_j t) dt = c_p (-1)**p J_(2p+lam)(a) a**(-lam), a_j >= 0."""
    p = np.arange(p_max + 1)[:, None]
    a = np.asarray(a, dtype=float)[None, :]
    safe = np.where(a > 0, a, 1.0)
    out = _gegenbauer_constant(p, lam) * (-1.0) ** p * jv(2 * p + lam, safe) * safe ** (-lam)
    at_zero = np.where(p == 0, _gegenbauer_constant(0, lam) * 2.0 ** (-lam) / math.gamma(1 + lam), 0.0)
    return np.where(a > 0, out, at_zero)


def inner_exterior(cyl, k, kn, basis):
    """E[p, n] = int_-h^-d u_p(t(z)) Z_n(z) dz, column 0 the propagating Z_0, columns 1 ... len(kn) the evanescent Z_n;
    shape (basis.size, len(kn) + 1), the families stacked."""
    g, h = cyl.gap, cyl.depth
    a = k * g
    ratio = 2.0 * math.exp(a - k * h) / (1.0 + math.exp(-2.0 * k * h))   # exp(a)/cosh(kh), without overflow
    kn = np.asarray(kn, dtype=float)
    rows = []
    for f in basis.families:
        p = np.arange(f.p_max + 1)
        propagating = _gegenbauer_constant(p, f.lam) * ive(2 * p + f.lam, a) * a ** (-f.lam) * ratio
        evanescent = _cos_integral(f.p_max, kn * g, f.lam) / np.cos(kn * h)[None, :]
        rows.append(np.column_stack([propagating, evanescent]))
    return 0.5 * g * np.concatenate(rows)


def inner_interior(cyl, n_max, basis):
    """F[p, n] = int_-h^-d u_p(t(z)) cos(lambda_n (z + h)) dz, lambda_n = n pi/g, n = 0 ... n_max;
    shape (basis.size, n_max + 1)."""
    lam_a = np.arange(n_max + 1) * math.pi
    return 0.5 * cyl.gap * np.concatenate([_cos_integral(f.p_max, lam_a, f.lam) for f in basis.families])


def _gauss_jacobi_exterior(cyl, k, kn, basis):
    """inner_exterior by the Gauss-Jacobi rules of the basis on the reflected interval, z + h = g abs(t). For the checks."""
    h = cyl.depth
    return basis.project(cyl, lambda s: np.array([np.cosh(k * s) / math.cosh(k * h)]
                                                 + [np.cos(q * s) / math.cos(q * h) for q in kn]))


def _gauss_jacobi_interior(cyl, n_max, basis):
    """inner_interior by the Gauss-Jacobi rules of the basis. For the checks."""
    return basis.project(cyl, lambda s: np.cos(np.arange(n_max + 1)[:, None] * math.pi * s[None, :] / cyl.gap))


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
    Gap, r = a: dphi_m/dr = sum_p gap[p] u_p(t). Wall, r = a, -d < z < 0: dphi_m/dr = wall_velocity[0] + wall_velocity[1] z,
    the body condition the order was solved with (zero for diffraction). The coefficients are amplitudes at r = a; the
    ratios are evaluated with hankel1, kve and ive and never overflow."""
    m: int
    incident: complex
    scattered: complex
    evanescent: np.ndarray
    interior: np.ndarray
    gap: np.ndarray
    particular: object = None
    wall_velocity: tuple = (0.0, 0.0)


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


def _inner_products(cyl, k, kn, basis):
    """(inner_exterior, inner_interior): the same for every order, so the solvers compute them once."""
    return inner_exterior(cyl, k, kn, basis), inner_interior(cyl, len(kn), basis)


MATRIX_TAIL = True   # add the sums over n > N of the leading asymptotic terms to the matching system (_matching_tail)


def _matching_tail(cyl, basis, n_modes):
    """(T, w): the sums over n > N of the non-oscillating leading terms of the matching system, by the Hurwitz zeta
    function. From the Hankel expansion, E_pn ~ (g/2) sqrt(2/pi) c_p (k_n g)**(-beta_p) cos(k_n g - chi_p)/cos(k_n h)
    and F_pn ~ (g/2) sqrt(2/pi) c_p (n pi)**(-beta_p) (-1)**n cos(chi_p), beta = lam + 1/2, chi = (lam/2 + 1/4) pi, with
    N_n kappa_n ~ -h k_n/(2 cos(k_n h)**2), nu_n mu_n ~ n pi/2, k_n ~ n pi/h, and the wall projection of a wall velocity
    V at the corner ~ -V sin(k_n g)/(k_n cos(k_n h)). T[p, q] is the tail of sum E E/(N kappa) - sum F F/(nu mu), w[p]
    that of sum E wall/(N kappa) per unit V. What is left falls one power of n faster (the O(1/n) terms of the expansions
    and the oscillating products, cos(2 k_n g - ...))."""
    g, h = cyl.gap, cyl.depth
    lam = np.concatenate([np.full(f.p_max + 1, f.lam) for f in basis.families])
    c = np.concatenate([_gegenbauer_constant(np.arange(f.p_max + 1), f.lam) for f in basis.families])
    beta, chi = lam + 0.5, (0.5 * lam + 0.25) * math.pi
    s = beta[:, None] + beta[None, :] + 1.0
    hurwitz = zeta(s, n_modes + 1.0)
    cc = np.outer(c, c)
    ext = -(g ** 2 / (2.0 * math.pi * h)) * cc * np.cos(chi[:, None] - chi[None, :])         * g ** (1.0 - s) * (h / math.pi) ** s * hurwitz
    intr = (g ** 2 / math.pi) * cc * np.outer(np.cos(chi), np.cos(chi)) * math.pi ** (-s) * hurwitz
    w = (g / (2.0 * h)) * math.sqrt(2.0 / math.pi) * c * np.sin(chi) * g ** (-beta) * (h / math.pi) ** (beta + 2.0)         * zeta(beta + 2.0, n_modes + 1.0)
    return ext - intr, w


def _solve_order(cyl, m, k, kn, basis, incident=0.0, wall=None, particular=None, inner=None, wall_velocity=(0.0, 0.0)):
    """One order of the matching; see the module docstring. A complex system of size basis.size (plus one for m = 0
    with several families). inner: _inner_products, if already computed."""
    a, g = cyl.radius, cyl.gap
    n_modes = len(kn)
    big_e, big_f = _inner_products(cyl, k, kn, basis) if inner is None else inner   # (size, N+1) each
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
    pot_p = np.zeros(basis.size) if particular is None else np.asarray(particular.potential_gap)

    ext = big_e / (norms * kappa)                           # E[p, n]/(N_n kappa_n)
    live = mu != 0.0                                        # m = 0, n = 0: no radial velocity, constant unknown
    intr = np.where(live, big_f / np.where(live, nu_int * mu, 1.0), 0.0)
    matrix = ext @ big_e.T - intr @ big_f.T
    if MATRIX_TAIL:
        tail, tail_wall = _matching_tail(cyl, basis, n_modes)
        matrix = matrix + tail
    # incident: E_p0 inc (J_m - J_m' H_m/H_m') = E_p0 inc 2i/(pi ka H_m'(ka)), by the Wronskian
    incident_at_wall = incident * 2j / (math.pi * ka * h1vp(m, ka))
    rhs = pot_p - big_e[:, 0] * incident_at_wall - ext @ wall - intr @ flux_p
    if MATRIX_TAIL:
        rhs = rhs - tail_wall * (wall_velocity[0] - cyl.draft * wall_velocity[1])
    if len(basis.families) > 1:
        # the families are not orthogonal and together nearly dependent: solve in the orthonormal combinations of
        # _gram_transform, gap = T c. For m = 0 the zero-flux condition sum_q gap[q] F[q, 0] = flux_p[0] borders the
        # system, with the constant interior[0] as one more unknown: [[M, -F0], [F0, 0]] (each family's p = 0 carries flux)
        trans = _gram_transform(basis)[0]
        reduced, rhs_r = trans.T @ matrix @ trans, trans.T @ rhs
        if m == 0:
            size = trans.shape[1]
            f0 = trans.T @ big_f[:, 0]
            bordered = np.zeros((size + 1, size + 1), complex)
            bordered[:size, :size] = reduced
            bordered[:size, size] = -f0
            bordered[size, :size] = f0
            sol = np.linalg.solve(bordered, np.concatenate([rhs_r, [flux_p[0]]]))
            gap, const = trans @ sol[:size], sol[size]
        else:
            gap, const = trans @ np.linalg.solve(reduced, rhs_r), None
    elif m == 0:
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
    inner = _inner_products(cyl, k, kn, basis)
    amplitude = -make_mesh.G / omega                       # phi_I = -(g A/omega) Z_0 e^(i k x), A = 1
    solved = {}
    for m in _order_list(orders):
        inc = amplitude * (1.0 if m == 0 else 2.0) * 1j ** m   # e^(ikx) = sum eps_m i^m J_m(kr) cos(m theta)
        solved[m] = _solve_order(cyl, m, k, kn, basis, incident=inc, inner=inner)
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
        self.potential_gap = basis.project(cyl, lambda s: self.fields(np.full(s.shape, a), s - h)[0])
        self.flux_gap = self.flux(n_modes)
        if mode == 3:
            self.bottom = (g ** 2 * a ** 2 / 2.0 - a ** 4 / 8.0) / (2.0 * g)
        else:
            self.bottom = -(g ** 2 * a ** 4 / 4.0 - a ** 6 / 24.0) / (2.0 * g)

    def flux(self, n_modes):
        """flux_gap[n], n = 0 ... n_modes, in closed form."""
        a, g = self.cyl.radius, self.cyl.gap
        lam = np.arange(1, n_modes + 1) * math.pi / g
        one = np.concatenate([[g], np.zeros(n_modes)])                                # int_0^g cos(lambda_n s) ds
        s2 = np.concatenate([[g ** 3 / 3.0], 2.0 * g * (-1.0) ** np.arange(1, n_modes + 1) / lam ** 2])   # s**2 cos
        if self.mode == 3:
            return -a / (2.0 * g) * one
        return -(s2 - 0.75 * a ** 2 * one) / (2.0 * g)

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
    velocity = {1: (1.0, 0.0), 3: (0.0, 0.0), 5: (0.0, 1.0)}[mode]   # dphi/dr on the wall, as `wall` above
    order = _solve_order(cyl, m, k, kn, basis, wall=wall, particular=particular, wall_velocity=velocity)
    order.wall_velocity = velocity
    return Solution(cyl, omega, k, kn, p_max, {m: order})


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


def body_matrices(cyl):
    """(M, C), 6 x 6, about the origin (the waterplane centre at the still waterline): make_mesh.hydrostatics's formulas
    for a uniform solid prism at its draft (steps = 1, exact for a prism), applied to the exact circle, area pi a**2 and
    second moments int x**2 dA = int y**2 dA = pi a**4/4. Mass rho pi a**2 d, z_G = z_B = -d/2, M15 = m z_G, M24 = -m z_G,
    M44 = M55 = rho d (pi a**4/4 + pi a**2 d**2/3), M66 = rho d pi a**4/2; C33 = rho g pi a**2,
    C44 = C55 = rho g (pi a**4/4 + V z_B) - m g z_G = rho g pi a**4/4. SI units, as in Hydrostatic.in."""
    a, d = cyl.radius, cyl.draft
    rho, g = make_mesh.RHO, make_mesh.G
    area, second = math.pi * a ** 2, math.pi * a ** 4 / 4.0
    vol = area * d
    mass = rho * vol
    zg = zb = -0.5 * d
    z2 = d ** 2 / 3.0                                  # (1/d) int_-d^0 z**2 dz
    m = np.zeros((6, 6))
    m[0, 0] = m[1, 1] = m[2, 2] = mass
    m[0, 4] = m[4, 0] = mass * zg
    m[1, 3] = m[3, 1] = -mass * zg
    m[3, 3] = m[4, 4] = rho * (second + area * z2) * d
    m[5, 5] = rho * 2.0 * second * d
    c = np.zeros((6, 6))
    c[2, 2] = rho * g * area
    c[3, 3] = c[4, 4] = rho * g * (second + vol * zb) - mass * g * zg
    return m, c


MOTION_MODES = (1, 3, 5)
_MODE_INDEX = [j - 1 for j in MOTION_MODES]


def _same_problem(name, cyl, omega, sol):
    if sol.cyl != cyl or sol.omega != omega:
        raise ValueError("%s: a solution was solved for another cylinder or frequency" % name)


def motions(cyl, omega, diff, rads, motion="free"):
    """xi_j, j = 1, 3, 5: the complex motion amplitudes per unit wave amplitude, time factor exp(-i omega t), rotations in
    rad per m (L = 1), about the origin, from [-omega**2 (M + A) - i omega B + C] xi = X with the dimensional A, B
    (added_mass_damping times rho and rho omega), X = i omega rho int phi n dS of the diffraction solution diff and
    (M, C) = body_matrices. motion "fixed" gives zeros (rads is not used). This is HAMS's equation of motion
    (SolveMotion in src/SolveMotion.f90, no external damping or stiffness) in the module's time convention; hams_motion
    gives the complex value Motion.4 writes, and the modulus abs(xi_j) is Motion.4's MOD."""
    if motion == "fixed":
        return {j: 0j for j in MOTION_MODES}
    if motion != "free":
        raise ValueError("motions: motion is 'fixed' or 'free', got %r" % (motion,))
    _same_problem("motions", cyl, omega, diff)
    if set(rads) != set(MOTION_MODES):
        raise ValueError("motions: rads needs the modes 1, 3 and 5")
    a, b = added_mass_damping(cyl, omega, rads)
    integrals = _body_integrals(cyl, diff)
    force = np.array([1j * omega * make_mesh.RHO * integrals[j] for j in MOTION_MODES])
    mass, restoring = body_matrices(cyl)
    rho = make_mesh.RHO
    lhs = -omega ** 2 * (mass + rho * a) - 1j * omega * (rho * omega * b) + restoring
    xi = np.linalg.solve(lhs[np.ix_(_MODE_INDEX, _MODE_INDEX)], force)
    return {j: complex(x) for j, x in zip(MOTION_MODES, xi)}


def long_wave_tolerance(cyl, omega, added_mass):
    """{3: tol, 5: tol} of check_long_wave_motions: omega**2 (M_jj + A_jj)/C_jj, the size of the inertia term against the
    restoring one, j = 3, 5. added_mass is added_mass_damping's A (over rho). The O(omega**2) correction to xi_3/eta_I(0)
    and xi_5/(-d eta_I/dx(0)) is that inertia term less the depth decay of the incident pressure over the draft and the
    body's own volume, terms of the same size that partly cancel (module docstring)."""
    mass, restoring = body_matrices(cyl)
    rho = make_mesh.RHO
    return {j: omega ** 2 * (mass[j - 1, j - 1] + rho * added_mass[j - 1, j - 1]) / restoring[j - 1, j - 1] for j in (3, 5)}


def hams_motion(xi):
    """Motion.4's NREL + i NIMG of a motion amplitude xi (exp(-i omega t), per unit wave amplitude): conj(xi/eta_I(0)),
    eta_I(0) = -i, the convention of ExcForce.3 (module docstring)."""
    return complex(np.conj(xi / -1j))


class _Combination:
    """A linear combination sum c_i p_i of particular solutions (the hook of Solution), for total_solution."""

    def __init__(self, terms):
        self.terms = terms
        self.potential_gap = sum(c * np.asarray(p.potential_gap) for c, p in terms)
        self.flux_gap = sum(c * np.asarray(p.flux_gap) for c, p in terms)
        self.bottom = sum(c * p.bottom for c, p in terms)

    def flux(self, n_modes):
        return sum(c * p.flux(n_modes) for c, p in self.terms)

    def fields(self, r, z):
        parts = [(c, p.fields(r, z)) for c, p in self.terms]
        return tuple(sum(c * f[i] for c, f in parts) for i in range(3))


def total_solution(diff, rads, xi):
    """phi_I + phi_D + sum_j (-i omega xi_j) phi_j: the total first-order potential of the free body, per unit wave
    amplitude, as a Solution with the orders of diff. The radiation potentials are per unit velocity and the velocity is
    -i omega xi in exp(-i omega t). xi from motions; if every xi_j is zero (fixed), rads may be None and diff is returned.
    rads must have the cylinder, frequency, N and P of diff."""
    if any(o.particular is not None for o in diff.orders.values()):
        # the combination below keeps only the radiation particular solutions; diffraction has none (solve_diffraction)
        raise ValueError("total_solution: diff has a particular solution, which would be dropped")
    velocity = {j: -1j * diff.omega * xi[j] for j in MOTION_MODES if xi.get(j, 0.0) != 0.0}
    if not velocity:
        return diff
    for j in velocity:
        rad = rads[j]
        _same_problem("total_solution", diff.cyl, diff.omega, rad)
        if rad.p_max != diff.p_max or not np.array_equal(rad.kn, diff.kn):
            raise ValueError("total_solution: rads[%s] has another N or P than diff" % j)
    orders = {}
    for m, o in diff.orders.items():
        parts = [(v, rads[j].orders[m]) for j, v in velocity.items() if RADIATION_ORDER[j] == m]
        if not parts:
            orders[m] = o
            continue
        particular = [(v, r.particular) for v, r in parts if r.particular is not None]
        orders[m] = Order(m, o.incident, o.scattered + sum(v * r.scattered for v, r in parts),
                          o.evanescent + sum(v * r.evanescent for v, r in parts),
                          o.interior + sum(v * r.interior for v, r in parts),
                          o.gap + sum(v * r.gap for v, r in parts),
                          _Combination(particular) if particular else None,
                          tuple(o.wall_velocity[i] + sum(v * r.wall_velocity[i] for v, r in parts) for i in (0, 1)))
    return Solution(diff.cyl, diff.omega, diff.k, diff.kn, diff.p_max, orders)


def group_ratio(k, depth):
    """C_g/C_p = (1 + 2 k h/sinh(2 k h))/2, safe at large k h (as group_ratio in src/FarFieldDrift.f90)."""
    x = 2.0 * k * depth
    return 0.5 * (1.0 + 2.0 * x * math.exp(-x) / (1.0 - math.exp(-2.0 * x)))


def far_field_fx(cyl, omega, total):
    """Mean surge drift F_x/(rho g A**2) by Maruo's formula, (1/2)(C_g/C_p) int_0^2pi abs(a)**2 (1 - cos theta) dtheta,
    a(theta) = sum_m a_m e^(i m theta) from far_field_amplitude of the total solution (the scattered and radiated waves),
    heading 0. Positive along +x, the wave direction, as MeanDriftMomentum.8 (module docstring). With a_-m = a_m the
    integral is 2 pi [sum_m abs(a_m)**2 - Re sum_m a_m conj(a_(m+1))], m = -M ... M."""
    _same_problem("far_field_fx", cyl, omega, total)
    a = far_field_amplitude(total)
    full = np.concatenate([a[:0:-1], a])                     # a_m, m = -M ... M
    energy = np.sum(np.abs(full) ** 2)
    cross = np.real(np.sum(full[:-1] * np.conj(full[1:])))
    return float(0.5 * group_ratio(total.k, cyl.depth) * 2.0 * math.pi * (energy - cross))


# Rule of the corner zone: "cube", Gauss-Legendre in (rho/L)**(1/3), the default, or "jacobi", Gauss-Jacobi for the weight
# rho**(-2/3) (module docstring, "Near field", for why the default is not the brief's Gauss-Jacobi rule)
CORNER_RULE = "cube"
CORNER_EXPONENT = -2.0 / 3.0   # abs(grad phi)**2 ~ rho**(-2/3) at the 270 degree corner
NEAR_THETA = 72                # trapezoidal points in theta: exact for trigonometric polynomials of degree < 72
NEAR_QUADRATURE = dict(n_corner=12, n_far=8)


@dataclass
class Face:
    """One face of the wetted surface on a meridian. Nodes (r, z); w, the weights of int f dS = int_0^2pi sum(w f) dtheta
    (a dz on the wall, r dr on the bottom); normal = (n_r, n_z), the unit normal out of the fluid (into the body), constant
    on the face; corner = distance of each node from the corner (a, -d)."""
    name: str
    r: np.ndarray
    z: np.ndarray
    w: np.ndarray
    normal: tuple
    corner: np.ndarray


@dataclass
class Quadrature:
    """The wetted surface: the wall r = a, -d < z < 0, the bottom z = -d, 0 < r < a, and the trapezoidal rule in theta,
    theta_k = 2 pi k/n, weight 2 pi/n, which is exact for the products of the near-field terms when n > 2 M + 4."""
    wall: Face
    bottom: Face
    theta: np.ndarray
    theta_weight: float

    @property
    def faces(self):
        return (self.wall, self.bottom)


def _corner_rule(length, n):
    """(rho, W): int_0^length f(rho) drho ~ sum(W f(rho)), Gauss-Jacobi for the weight rho**(-2/3) (exact for f =
    rho**(-2/3) times a polynomial of degree < 2 n). scipy's weight (1 - x)**0 (1 + x)**beta, rho = length (1 + x)/2."""
    if CORNER_RULE == "cube":
        s, w = _legendre_rule(0.0, 1.0, n)
        return length * s ** 3, 3.0 * length * s ** 2 * w
    x, w = roots_jacobi(n, 0.0, CORNER_EXPONENT)
    return 0.5 * length * (1.0 + x), 0.5 * length * w * (1.0 + x) ** (-CORNER_EXPONENT)


def _legendre_rule(lo, hi, n):
    """Gauss-Legendre nodes and weights on [lo, hi]."""
    x, w = np.polynomial.legendre.leggauss(n)
    return lo + 0.5 * (hi - lo) * (1.0 + x), 0.5 * (hi - lo) * w


def body_quadrature(cyl, n_corner, n_far, n_theta=NEAR_THETA):
    """Quadrature of the wetted surface. On each face, Gauss-Jacobi with the weight rho**(-2/3) in the distance rho from
    the corner, n_corner nodes, within min(gap, d) of it on the wall and min(gap, a/2) on the bottom, and Gauss-Legendre with
    n_far nodes beyond (on the default cylinder the corner zone is the whole wall, d = 10 < g = 90, so the wall has no
    Gauss-Legendre part)."""
    a, d, g = cyl.radius, cyl.draft, cyl.gap
    faces = []
    for name, extent, zone in (("wall", d, min(g, d)), ("bottom", a, min(g, 0.5 * a))):
        rho, w = _corner_rule(zone, n_corner)
        if extent > zone:
            far, wf = _legendre_rule(zone, extent, n_far)
            rho, w = np.concatenate([rho, far]), np.concatenate([w, wf])
        if name == "wall":
            faces.append(Face(name, np.full(rho.shape, a), -d + rho, a * w, (-1.0, 0.0), rho))
        else:
            faces.append(Face(name, a - rho, np.full(rho.shape, -d), (a - rho) * w, (0.0, 1.0), rho))
    theta = 2.0 * math.pi * np.arange(n_theta) / n_theta
    return Quadrature(faces[0], faces[1], theta, 2.0 * math.pi / n_theta)


NEAR_FIELD_TERMS = ("I", "Ia", "II", "III", "IV", "V")   # the order of MeanDriftTerms.txt (T_WATERLINE ... T_HYDROSTATIC)


@dataclass
class _SurfaceState:
    """One first-order state on the quadrature, per unit wave amplitude, exp(-i omega t). Arrays over (node, theta):
    per face, phi, grad (3, ...), the displacement disp (3, ...) = xi + alpha x (x - xr); at the waterline r = a, z = 0,
    phi_wl and disp3_wl over theta. xi = translation of xr, alpha = rotation, 3-vectors."""
    omega: float
    xi: np.ndarray
    alpha: np.ndarray
    faces: list
    phi_wl: np.ndarray
    disp3_wl: np.ndarray


def _rigid(xi, xr):
    """(xi', alpha): the translation of the point xr and the rotation of the module's motions xi (about the origin),
    xi' = xi + alpha x xr, so that the displacement xi + alpha x x = xi' + alpha x (x - xr)."""
    alpha = np.array([0.0, xi.get(5, 0.0), 0.0], complex)
    trans = np.array([xi.get(1, 0.0), 0.0, xi.get(3, 0.0)], complex)
    return trans + np.cross(alpha, xr), alpha


NEAR_SERIES_MODES = 2 ** 15   # modes of the extended series on the body (_body_series), at least the N of the solution
NEAR_WINDOW_START = 0.25      # the smooth window of the wall series falls from 1 to 0 over n/NEAR_SERIES_MODES in [0.25, 1]
BOTTOM_DECAY = 40.0           # bottom series: modes with lambda_n (a - r) > 40 are dropped (exp(-40) of the first)


def _extended_roots(sol, count):
    """k_n, n = 1 ... count: sol.kn, then y = n pi - atan(nu h/y), y = k_n h, by fixed-point iteration (its derivative
    nu h/(y**2 + (nu h)**2) is below 1e-6 for the n beyond any N used here, so eight steps reach round-off)."""
    h = sol.cyl.depth
    c = sol.omega ** 2 / make_mesh.G * h
    n = np.arange(len(sol.kn) + 1, count + 1)
    y = (n - 0.5) * math.pi
    for _ in range(8):
        y = n * math.pi - np.arctan(c / y)
    return np.concatenate([sol.kn, y / h])


def _smooth_window(count, start):
    """w_n, n = 1 ... count: 1 up to n = start count, then a C-infinity step down to 0 at n = count."""
    x = np.clip((np.arange(1, count + 1) - start * count) / ((1.0 - start) * count), 0.0, 1.0)
    inside = (x > 0.0) & (x < 1.0)
    e0 = np.where(inside, np.exp(-1.0 / np.where(inside, x, 1.0)), 0.0)
    e1 = np.where(inside, np.exp(-1.0 / np.where(inside, 1.0 - x, 1.0)), 0.0)
    return np.where(x <= 0.0, 1.0, np.where(x >= 1.0, 0.0, e1 / np.where(inside, e0 + e1, 1.0)))


POLYLOG_TERMS = 40   # terms of the zeta series of _polylog; |mu| <= 0.2 pi here, so they fall like 0.1**k


def _polylog(order, mu):
    """Li_s(e^mu) = sum_(n >= 1) n**(-s) e^(n mu), Re mu <= 0, 0 < |mu| < 2 pi, elementwise. Non-integer s: the zeta series
    Li_s(e^mu) = Gamma(1 - s) (-mu)**(s - 1) + sum_k zeta(s - k) mu**k/k! (DLMF 25.12.12); s = 1: -log(1 - e^mu)."""
    mu = np.asarray(mu, complex)
    if order == 1:
        return -np.log(-np.expm1(mu))
    out = gamma(1.0 - order) * (-mu) ** (order - 1.0)
    term = np.ones_like(mu)
    for k in range(POLYLOG_TERMS):
        out = out + zeta(order - k) * term
        term = term * mu / (k + 1)
    return out


def _unit_phase(phi):
    """i (phi reduced to (-pi, pi]): the argument of _polylog for e^(i phi)."""
    return 1j * (np.mod(np.asarray(phi) + math.pi, 2.0 * math.pi) - math.pi)


def _bottom_keep(quad, cyl, count):
    """Modes kept per bottom node: lambda_n (a - r) <= BOTTOM_DECAY, at most count."""
    rho = np.maximum(quad.bottom.corner, 1e-300)
    return np.minimum(count, np.ceil(BOTTOM_DECAY * cyl.gap / (math.pi * rho))).astype(int)


def _corner_models(quad, cyl, count, window, lam):
    """Per node, the exact sum less the windowed (wall) or cut (bottom) sum of the model tails of the family lam of
    _body_series, beta = lam + 1/2, chi = (lam/2 + 1/4) pi: on the wall "dz" = sum (n pi g/h)**(-beta)
    [sin(n pi sigma/h - chi) + sin(n pi rho/h + chi)]/2 (sigma = s + g), "phi" the same with h/(n pi) and cos; on the
    bottom "bottom_dr" = sqrt(a/r) sum (n pi)**(-beta) exp(-n pi rho/g) and "bottom_phi" with g/(n pi). lam = 1/2 gives
    "jump" only, sum [cos(n pi rho/h) - cos(n pi sigma/h)]/(n pi): its other tails converge absolutely (n**(-2)) or vanish
    (cos chi = 0). Module docstring, "Near field"."""
    h, g, a = cyl.depth, cyl.gap, cyl.radius
    n = np.arange(1, count + 1)
    rho = quad.wall.corner
    sigma = quad.wall.z + h + g                       # s + g
    near, far = np.outer(rho, n) * (math.pi / h), np.outer(sigma, n) * (math.pi / h)
    pn, pf = _unit_phase(rho * math.pi / h), _unit_phase(sigma * math.pi / h)
    if lam == 0.5:
        windowed_jump = (np.cos(near) - np.cos(far)) @ (window / (n * math.pi))
        exact_jump = (np.real(_polylog(1, pn)) - np.real(_polylog(1, pf))) / math.pi
        return dict(jump=exact_jump - windowed_jump)
    beta, chi = lam + 0.5, (0.5 * lam + 0.25) * math.pi
    amp = (n * math.pi * g / h) ** (-beta) * window
    windowed_dz = 0.5 * (np.sin(far - chi) + np.sin(near + chi)) @ amp
    windowed_phi = 0.5 * (np.cos(near + chi) + np.cos(far - chi)) @ (amp * h / (n * math.pi))
    scale = (math.pi * g / h) ** (-beta)
    e = np.exp(1j * chi)
    exact_dz = 0.5 * scale * (np.imag(np.conj(e) * _polylog(beta, pf)) + np.imag(e * _polylog(beta, pn)))
    exact_phi = 0.5 * scale * h / math.pi * (np.real(e * _polylog(beta + 1.0, pn))
                                             + np.real(np.conj(e) * _polylog(beta + 1.0, pf)))
    rho_b = quad.bottom.corner
    tail = np.exp(-np.outer(rho_b, n) * (math.pi / g)) * (n * math.pi) ** (-beta)
    keep = _bottom_keep(quad, cyl, count)
    cut_dr = np.array([np.sum(tail[i, :q]) for i, q in enumerate(keep)])
    cut_phi = np.array([np.sum(tail[i, :q] * g / (n[:q] * math.pi)) for i, q in enumerate(keep)])
    mu_b = -rho_b * math.pi / g
    exact_dr = math.pi ** (-beta) * np.real(_polylog(beta, mu_b))
    exact_bphi = math.pi ** (-beta) * g / math.pi * np.real(_polylog(beta + 1.0, mu_b))
    root = np.sqrt(a / quad.bottom.r)
    return dict(dz=exact_dz - windowed_dz, phi=exact_phi - windowed_phi,
                bottom_dr=root * (exact_dr - cut_dr), bottom_phi=root * (exact_bphi - cut_phi))


def _body_series(sol, quad, count=None):
    """Per order m: (phi_m, d phi_m/dz) on the wall nodes, (phi_m, d phi_m/dr) on the bottom nodes, phi_m at the
    waterline, of the extended, corner-corrected series (module docstring, "Near field"). The gap velocity
    sum_p gap[p] u_p and the wall velocity of each order define the exterior and interior coefficients for every n in
    closed form (as in _solve_order); the first N are the solution's own, and the series is carried to n = count. The
    wall and waterline sums are summed with _smooth_window, the bottom sums cut where exp(-lambda_n (a - r)) is below
    exp(-BOTTOM_DECAY); then the exact sum less the windowed (or cut) sum of each model tail is added (_corner_models)."""
    cyl = sol.cyl
    a, d, h, g = cyl.radius, cyl.draft, cyl.depth, cyl.gap
    count = max(NEAR_SERIES_MODES, len(sol.kn)) if count is None else count
    if count < len(sol.kn):
        raise ValueError("_body_series: count %d is below the solution's N = %d" % (count, len(sol.kn)))
    n_sol = len(sol.kn)
    kn = _extended_roots(sol, count)
    basis = gegenbauer_basis(sol.p_max)
    big_e = inner_exterior(cyl, sol.k, kn, basis)[:, 1:]
    big_f = inner_interior(cyl, count, basis)[:, 1:]
    norms = _mode_norms(sol.k, kn, h)[1:]
    i0, i1 = _wall_moments(cyl, sol.k, kn)
    window = _smooth_window(count, NEAR_WINDOW_START)
    models = {f.lam: _corner_models(quad, cyl, count, window, f.lam) for f in basis.families}
    if 0.5 not in models:                                        # the jump of the wall velocity needs it
        models[0.5] = _corner_models(quad, cyl, count, window, 0.5)
    s_wall = quad.wall.z + h
    secant = 1.0 / np.cos(kn * h)
    phase = np.outer(s_wall, kn)
    cos_wall = np.cos(phase) * (secant * window)
    sin_wall = -np.sin(phase) * (kn * secant * window)
    z0, dz0 = _cosh_ratio(sol.k, s_wall, h)
    dz0 = sol.k * dz0
    lam = np.arange(1, count + 1) * math.pi / g
    sign = (-1.0) ** np.arange(1, count + 1)                     # cos(lambda_n g), the bottom z = -d
    rb = quad.bottom.r
    keep = _bottom_keep(quad, cyl, count)
    top = max(sol.orders) + 2                                    # ive(k, .) for k = 0 ... max order + 1, shared by the orders
    bessel_a = [ive(k, lam * a) for k in range(top)]
    bessel_r = [[ive(k, lam[:q] * r) for k in range(top)] for r, q in zip(rb, keep)]
    growth = [np.exp(lam[:q] * (r - a)) for r, q in zip(rb, keep)]
    out = {}
    for m, o in sol.orders.items():
        x = kn * a
        kappa = -0.5 * kn * (kve(m - 1, x) + kve(m + 1, x)) / kve(m, x)
        wall = o.wall_velocity[0] * i0[1:] + o.wall_velocity[1] * i1[1:]
        b = (o.gap @ big_e + wall) / (norms * kappa)
        y = lam * a
        mu = 0.5 * lam * (ive(m - 1, y) + ive(m + 1, y)) / ive(m, y)
        flux = o.particular.flux(count)[1:] if o.particular is not None else 0.0
        c = (o.gap @ big_f - flux) / (0.5 * g * mu)
        for mine, theirs, name in ((b[:n_sol], o.evanescent, "exterior"), (c[:n_sol], o.interior[1:], "interior")):
            if n_sol and np.max(np.abs(mine - theirs)) > 1e-8 * np.max(np.abs(theirs)):
                raise ValueError("_body_series: the %s coefficients of order %d differ from the solution's" % (name, m))
        # tail strengths S = sum_p gap[p] c_p per family; the lam = 1/2 family's corner value sum_p gap[p] (C_2p^(1/2)(1)
        # = 1) is the regular gap velocity at the corner, so the jump of the radial velocity there is v_corner minus it
        strength = {f.lam: o.gap[sl] @ _gegenbauer_constant(np.arange(f.p_max + 1), f.lam) for f, sl in basis.slices}
        v_corner = o.wall_velocity[0] - d * o.wall_velocity[1]   # wall velocity at the corner
        v_corner = v_corner - 0.5 * math.sqrt(2.0 / math.pi) * strength.get(0.5, 0.0)
        radial = o.incident * jv(m, sol.k * a) + o.scattered     # the propagating amplitude at r = a
        wall_phi = z0 * radial + cos_wall @ b
        wall_dz = dz0 * radial + sin_wall @ b - v_corner * models[0.5]["jump"]
        bottom_model_phi = bottom_model_dr = 0.0
        for lam_f, s_f in strength.items():
            if lam_f == 0.5:
                continue
            corner = (g / h) * math.sqrt(2.0 / math.pi) * s_f
            wall_phi = wall_phi - corner * models[lam_f]["phi"]
            wall_dz = wall_dz + corner * models[lam_f]["dz"]
            k_c = math.sqrt(2.0 / math.pi) * s_f * math.cos((0.5 * lam_f + 0.25) * math.pi)
            bottom_model_phi = bottom_model_phi + k_c * models[lam_f]["bottom_phi"]
            bottom_model_dr = bottom_model_dr + k_c * models[lam_f]["bottom_dr"]
        waterline = radial + np.sum(b * window)                  # Z_n(0) = 1
        bottom_phi = o.interior[0] * (rb / a) ** m + 0j
        bottom_dr = o.interior[0] * m * rb ** max(m - 1, 0) / a ** m + 0j
        for i in range(rb.size):
            q = keep[i]
            at_r, at_a = bessel_r[i], bessel_a
            grow = growth[i] / at_a[m][:q]
            cs = c[:q] * sign[:q]
            bottom_phi[i] += np.sum(cs * at_r[m] * grow)
            bottom_dr[i] += np.sum(cs * 0.5 * lam[:q] * (at_r[abs(m - 1)] + at_r[m + 1]) * grow)
        bottom_phi = bottom_phi + bottom_model_phi
        bottom_dr = bottom_dr + bottom_model_dr
        if o.particular is not None:
            pp, pr, _ = o.particular.fields(rb, np.full(rb.shape, -d))
            bottom_phi, bottom_dr = bottom_phi + pp, bottom_dr + pr
        out[m] = dict(wall=(wall_phi, wall_dz), bottom=(bottom_phi, bottom_dr), waterline=waterline)
    return out


def _surface_state(cyl, omega, sol, xi, quad, xr):
    """phi and grad phi of sol on the quadrature. Normal derivative: the body condition, d phi/dn = -i omega disp.n (n out
    of the fluid), on the wall d phi/dr = -d phi/dn, on the bottom d phi/dz = d phi/dn. Tangential derivatives and phi:
    the extended series of _body_series, exterior on the wall and at the waterline, interior on the bottom."""
    _same_problem("near field", cyl, omega, sol)
    trans, alpha = _rigid(xi, xr)
    ct, st = np.cos(quad.theta), np.sin(quad.theta)
    series = _body_series(sol, quad)
    faces = []
    for face in quad.faces:
        phi = np.zeros((face.r.size, ct.size), complex)
        g_t = phi.copy()
        g_s = phi.copy()                                       # the series' tangential derivative on the meridian
        for m in sol.orders:
            pm, dm = series[m][face.name]
            cm, sm = np.cos(m * quad.theta), np.sin(m * quad.theta)
            phi += np.outer(pm, cm)
            g_t += np.outer(-m * pm / face.r, sm)
            g_s += np.outer(dm, cm)
        x = np.array([np.outer(face.r, ct), np.outer(face.r, st), np.outer(face.z, np.ones_like(ct))])
        disp = trans[:, None, None] + np.cross(alpha, x - np.asarray(xr)[:, None, None], axis=0)
        n_r, n_z = face.normal
        normal = np.array([n_r * ct[None, :] + 0 * x[0], n_r * st[None, :] + 0 * x[0], n_z + 0 * x[0]])
        g_n = -1j * omega * np.sum(disp * normal, axis=0)      # d phi/dn, the body condition
        if face.name == "wall":
            g_r, g_z = -g_n, g_s                               # n = -e_r
        else:
            g_r, g_z = g_s, g_n                                # n = +e_z
        grad = np.array([ct * g_r - st * g_t, st * g_r + ct * g_t, g_z])
        faces.append(dict(face=face, phi=phi, grad=grad, disp=disp, x=x - np.asarray(xr)[:, None, None],
                          normal=normal))
    a = cyl.radius
    phi_wl = sum(series[m]["waterline"] * np.cos(m * quad.theta) for m in sol.orders)
    x_wl = np.array([a * ct, a * st, 0.0 * ct]) - np.asarray(xr)[:, None]
    disp3_wl = trans[2] + np.cross(alpha, x_wl, axis=0)[2]
    return _SurfaceState(omega, trans, alpha, faces, phi_wl, disp3_wl)


def _pair(fi, gi, fj, gj):
    """The kernel T- of spec 5.4 without c/2: (fi conj(gj) + conj(fj) gi)/4, elementwise; for i = j, Re(f conj(g))/2."""
    return 0.25 * (fi * np.conj(gj) + np.conj(fj) * gi)


def _pair_cross(ai, bi, aj, bj):
    """(a x b) under the kernel: component l = sum eps_lmn pair(a_m, b_n), 3-vectors."""
    out = np.zeros(3, complex)
    for l, m, n in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        out[l] = _pair(ai[m], bi[n], aj[m], bj[n]) - _pair(ai[n], bi[m], aj[n], bj[m])
    return out


def _load(vector, x):
    """The force and moment densities (6, ...) about xr of a force density (3, ...) at the points x - xr (3, ...); the
    caller integrates them."""
    return np.concatenate([vector, np.cross(x, vector, axis=0)])


def near_field_terms(cyl, omega_i, sol_i, xi_i, omega_j, sol_j, xi_j, quad, xr=(0.0, 0.0, 0.0)):
    """T-_ij of each near-field term (spec 5.3 and 5.4 of docs/superpowers/specs/2026-11-04-second-order-waves-design.md),
    over rho g, per unit wave amplitudes, as (6,) complex arrays (F_x, F_y, F_z, M_x, M_y, M_z, moments about xr), components
    1, 3, 5 filled; for i = j the real part is the mean load over rho g A**2. xi_i, xi_j: the motions (module convention,
    about the origin); xr in the x-z plane. See the module docstring for the terms as coded."""
    xr = np.asarray(xr, dtype=float)
    if xr[1] != 0.0:
        raise ValueError("near_field_terms: xr must lie in the x-z plane (F_y, M_x, M_z are not computed)")
    orders = max(max(sol_i.orders), max(sol_j.orders))
    if quad.theta.size <= 2 * orders + 4:
        raise ValueError("near_field_terms: %d theta points are too few for order %d" % (quad.theta.size, orders))
    si = _surface_state(cyl, omega_i, sol_i, xi_i, quad, xr)
    same = sol_j is sol_i and omega_j == omega_i and xi_j == xi_i
    sj = si if same else _surface_state(cyl, omega_j, sol_j, xi_j, quad, xr)
    grav = make_mesh.G
    dth = quad.theta_weight
    terms = {key: np.zeros(6, complex) for key in NEAR_FIELD_TERMS}

    # I: (1/2) int_WL r' [P(zeta, zeta) - 2 P(zeta, X3)] dl, zeta = (i omega/g) phi, r' = n' or (x - xr) x n', dl = a dtheta
    a = cyl.radius
    ct, st = np.cos(quad.theta), np.sin(quad.theta)
    n_wl = np.array([-ct, -st, 0.0 * ct])                     # unit horizontal normal into the body
    x_wl = np.array([a * ct, a * st, 0.0 * ct]) - xr[:, None]
    zi, zj = 1j * omega_i / grav * si.phi_wl, 1j * omega_j / grav * sj.phi_wl
    bracket = _pair(zi, zi, zj, zj) - 2.0 * _pair(zi, si.disp3_wl, zj, sj.disp3_wl)
    terms["I"] = 0.5 * np.sum(_load(n_wl * bracket, x_wl), axis=1) * a * dth
    # Ia: (1/2) int_WL k n_z/sqrt(1 - n_z**2) P(zeta - X3, zeta - X3) dl is zero: the wall is vertical at the waterline

    gi_t = [-1j * omega_i * si_f["grad"] for si_f in si.faces]   # grad phi_t, phi_t = -i omega phi
    gj_t = [-1j * omega_j * sj_f["grad"] for sj_f in sj.faces]
    big_g = [np.zeros(3, complex), np.zeros(3, complex)]          # G/(rho g) = -(1/g) int n phi_t dS, state i and j
    big_h = [np.zeros(3, complex), np.zeros(3, complex)]          # H/(rho g) = -(1/g) int (x - xr) x n phi_t dS
    for f, (fi, fj) in enumerate(zip(si.faces, sj.faces)):
        w = fi["face"].w[:, None] * dth
        n, x = fi["normal"], fi["x"]
        # II: -(1/g) int r (1/2) sum_k P(d_k phi, d_k phi) dS
        q = 0.5 * np.sum(_pair(fi["grad"], fi["grad"], fj["grad"], fj["grad"]), axis=0)
        terms["II"] += -np.sum(_load(n * q, x) * w, axis=(1, 2)) / grav
        # III: -(1/g) int r sum_k P(X_k, d_k phi_t) dS
        q = np.sum(_pair(fi["disp"], gi_t[f], fj["disp"], gj_t[f]), axis=0)
        terms["III"] += -np.sum(_load(n * q, x) * w, axis=(1, 2)) / grav
        for s, (state, omega) in enumerate(((fi, omega_i), (fj, omega_j))):
            dens = _load(n * (-1j * omega * state["phi"]), x)
            total = -np.sum(dens * w, axis=(1, 2)) / grav
            big_g[s] += total[:3]
            big_h[s] += total[3:]
    # IV: force alpha x G; moment alpha x H + xi x G (xi the translation of xr)
    terms["IV"][:3] = _pair_cross(si.alpha, big_g[0], sj.alpha, big_g[1])
    terms["IV"][3:] = (_pair_cross(si.alpha, big_h[0], sj.alpha, big_h[1])
                       + _pair_cross(si.xi, big_g[0], sj.xi, big_g[1]))
    terms["V"] = _hydrostatic_term(cyl, si, sj, xr)
    for key in terms:                                          # components 2, 4, 6 vanish by symmetry
        terms[key][[1, 3, 5]] = 0.0
    return terms


def _hydrostatic_term(cyl, si, sj, xr):
    """V: F_S and M_S of Task 7 of docs/superpowers/plans/2026-11-04-mean-wave-drift.md (WAMIT eq. 15.56, 15.57), over
    rho g, under the kernel, with xi = translation of xr, alpha = rotation, waterplane and buoyancy quantities relative to
    xr: area A, centre (x_f, y_f), L_ij = int (x - xr)_i (x - xr)_j dA, volume V, centre of buoyancy (x_b, y_b),
    Z_o = z of xr. F_S,x = F_S,y = 0; only F_S,z and M_S,y are filled."""
    a, d = cyl.radius, cyl.draft
    area = math.pi * a ** 2
    xf, yf, zo = -xr[0], -xr[1], xr[2]
    l11 = math.pi * a ** 4 / 4.0 + area * xf ** 2
    l22 = math.pi * a ** 4 / 4.0 + area * yf ** 2
    l12 = area * xf * yf
    vol, xb = area * d, -xr[0]

    x1, x3 = (si.xi[0], sj.xi[0]), (si.xi[2], sj.xi[2])
    al1, al2, al3 = ((si.alpha[k], sj.alpha[k]) for k in range(3))

    def pp(u, v):   # P(u, v) with u from state i, j: (u_i, u_j), (v_i, v_j)
        return _pair(u[0], v[0], u[1], v[1])

    out = np.zeros(6, complex)
    out[2] = -area * (pp(al1, al3) * xf + pp(al2, al3) * yf + 0.5 * (pp(al1, al1) + pp(al2, al2)) * zo)
    out[4] = (-pp(x1, al2) * area * xf + (pp(x1, al1) - pp(x3, al3)) * area * yf
              + (0.5 * pp(al1, al1) + 1.5 * pp(al2, al2)) * area * xf * zo - pp(al1, al2) * area * yf * zo
              + pp(x1, x3) * area - pp(x3, al2) * area * zo
              + pp(al1, al3) * (l11 - l22) + 2.0 * pp(al2, al3) * l12 + 0.5 * (pp(al2, al2) + pp(al3, al3)) * vol * xb)
    return out


def near_field(cyl, omega, total, xi, quad, xr=(0.0, 0.0, 0.0)):
    """(terms, total): the diagonal of near_field_terms, real, over rho g A**2; total is their sum, (6,), the normalisation
    of MeanDriftPressure.9 (conj(F)/(rho g REFL**k), REFL = 1, write_drift_row in src/SecondOrderOutput.f90; F is real)."""
    terms = near_field_terms(cyl, omega, total, xi, omega, total, xi, quad, xr)
    terms = {key: value.real for key, value in terms.items()}
    return terms, sum(terms.values())


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


def _quad_gap(cyl, p, mode, lam=LAMBDA):
    """int over -h < z < -d of u_p(t(z)) mode(z) dz, u_p of the family lam, by QUADPACK on the unreflected gap, the corner
    weight (-d - z)**(lam - 1/2) handled by weight="alg"; (1 - t)**(lam - 1/2) = g**(1/2 - lam) (-d - z)**(lam - 1/2), and
    (1 + t)**(lam - 1/2) is smooth on the gap."""
    g, h, d = cyl.gap, cyl.depth, cyl.draft
    alpha = lam - 0.5

    def f(z):
        t = (z + h) / g
        return eval_gegenbauer(2 * p, lam, t) * (1.0 + t) ** alpha * mode(z)

    # epsabs: entries with p large against a slowly varying mode are near zero by cancellation, and a smaller absolute
    # tolerance is below their rounding floor (QUADPACK warns). The matrices' largest entries are above 60.
    with warnings.catch_warnings():
        # QUADPACK reports round-off on a few near-zero entries of the lam = 1/2 and 5/6 families; the closed form and
        # the Gauss-Jacobi sums agree with these values to 1e-13 of the largest entry (check_inner_products)
        warnings.simplefilter("ignore", IntegrationWarning)
        value, _ = quad(f, -h, -d, weight="alg", wvar=(0.0, alpha), epsabs=1e-12, epsrel=1e-12, limit=400)
    return g ** (-alpha) * value


def _relative(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


def _relative_families(mine, ref, basis):
    """_relative per family of the basis (rows), the largest."""
    return max(_relative(mine[sl], ref[sl]) for _, sl in basis.slices)


def check_inner_products():
    """inner_exterior, inner_interior and the Gauss-Jacobi sums against quad, every family of FAMILIES with p <= CHECK_P;
    see the module docstring."""
    cyl = Cylinder()
    h = cyl.depth
    basis = gegenbauer_basis((CHECK_P,) * len(FAMILIES))

    def reference(modes):
        return np.concatenate([np.array([[_quad_gap(cyl, p, m, f.lam) for m in modes] for p in range(f.p_max + 1)])
                               for f in basis.families])

    _relative = functools.partial(_relative_families, basis=basis)
    worst = {"exterior": 0.0, "interior": 0.0, "exterior GJ": 0.0, "interior GJ": 0.0}
    for omega in drift_convergence.WIDE_FREQUENCIES:
        k, kn = wave_numbers(omega, h, CHECK_N)
        modes = [lambda z: math.cosh(k * (z + h)) / math.cosh(k * h)]
        modes += [lambda z, q=q: math.cos(q * (z + h)) / math.cos(q * h) for q in kn]
        ref = reference(modes)
        worst["exterior"] = max(worst["exterior"], _relative(inner_exterior(cyl, k, kn, basis), ref))
        worst["exterior GJ"] = max(worst["exterior GJ"], _relative(_gauss_jacobi_exterior(cyl, k, kn, basis), ref))
    lam = np.arange(CHECK_N + 1) * math.pi / cyl.gap
    modes = [lambda z, q=q: math.cos(q * (z + h)) for q in lam]
    ref = reference(modes)
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
    see the module docstring. At LIGHT_RESOLUTION: the extrapolated error (5.5e-5) and the ratio are the same at N = 500,
    P = 10 and N = 8000, P = 20, so the finer resolution only costs time."""
    h = Cylinder().depth
    res = LIGHT_RESOLUTION
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
    """Fixed body: scattered energy flux against the optical theorem; see the module docstring. At LIGHT_RESOLUTION: the
    identity holds order by order at any truncation (1e-14 at N = 500 and 8000), so the finer resolution only costs
    time."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
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
    particular solution enters, Green's identities held for the truncated series exactly: surge's Haskind relation held
    to round-off at any N and P (1e-16 at N = 50, P = 5). MATRIX_TAIL (Matching tail) adds the tails of the sums to the
    solve but not to the forces, so it breaks that exactness: surge's Haskind residual is now 7e-4 at N = 50, P = 5, and
    5e-8, 5e-9 at N = 16000, 64000 (gap basis (20, 4, 4), 0.70 rad/s). It still tests the Haskind formula, the far-field
    normalisation and the sign and conjugation of excitation and added_mass_damping, now after extrapolation in N. Pitch and heave have
    the interior particular solution, which is not in the span of the N truncated interior modes; the residual of the
    identities is the tail of its cosine projection (coefficients n**-2) against that of the corner-singular gap velocity
    (n**-2/3), so it falls as N**(-5/3): measured exponent 1.68 to 1.73 for N = 1000 ... 8000, 1.3e-7 for A15 = A51 still at
    N = 32000. One Richardson step with p = 5/3 removes that term; what remains falls about as N**-2.9 (2e-8 at N = 4000,
    8000). So A15 = A51, B15 = B51 and Haskind for modes 3 and 5 test the particular solutions through their gap
    projections (potential_gap, flux_gap) and the pitch wall term. They do not test _Particular.bottom: the bottom
    integral of phi_p is real and enters only the diagonal Re A33 and Re A55, which no identity here contains;
    check_particular_bottom tests it against quad. They do not bound the error of A, B or X themselves, which converge as
    about N**-1.6 (Convergence)."""
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
    held to round-off at any truncation before MATRIX_TAIL (as check_reciprocity says; now it converges in N like the
    others), so for it the check tests the flux formula and the normalisation and sign of B. Heave and pitch have the residual
    N**(-5/3) of the particular solution and test it. It does not test the convergence of B."""
    def error(v):
        return np.max(np.abs(np.diag(v["b"])[[0, 2, 4]] - v["flux"]) / np.abs(v["flux"]))

    worst = raw = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        q, hi = _first_order_extrapolated(omega)
        worst, raw = max(worst, error(q)), max(raw, error(hi))
    return worst < 1e-6, ("extrapolated in N: max |B_jj - flux|/|flux| %.1e, j = 1, 3, 5; at N = %d itself %.1e"
                          % (worst, 2 * CHECK_RESOLUTION["n_modes"], raw))


def check_particular_bottom():
    """_Particular.bottom, the closed-form int_0^a phi_p(r, -d) r**(m+1) dr, equals scipy.integrate.quad of
    fields(r, -d)[0] r**(m+1) to 1e-12 relative, heave and pitch, default cylinder and a gap of 1e-3 h. It is the only
    test of the bottom term: it is real, so it enters Re A33 and Re A55 only, which no identity check sees (module
    docstring). It does not test fields itself; reciprocity and radiation_energy do, through the gap projections."""
    worst = 0.0
    h = Cylinder().depth
    basis = gegenbauer_basis(2)
    for cyl in (Cylinder(), Cylinder(draft=h - BOTTOM_GAP * h)):
        for mode in (3, 5):
            p = _Particular(cyl, mode, basis, 4)
            m = RADIATION_ORDER[mode]
            ref, _ = quad(lambda r: float(p.fields(r, -cyl.draft)[0]) * r ** (m + 1), 0.0, cyl.radius,
                          epsabs=0.0, epsrel=1e-13)
            worst = max(worst, abs(p.bottom - ref) / abs(ref))
    return worst < 1e-12, "max |bottom - quad|/|quad| %.1e, heave and pitch, gaps 90 m and 0.1 m" % worst


def check_bottom_limit_drift():
    """Gap 1e-3 h, fixed: far_field_fx against drift_convergence.maccamy_fuchs_fx, after extrapolation in the gap,
    2 F(g) - F(2g), to 1e-3; see the module docstring."""
    h = Cylinder().depth
    res = LIGHT_RESOLUTION
    worst = worst_raw = 0.0
    ratios, skipped = [], []
    for omega in drift_convergence.WIDE_FREQUENCIES:
        f = []
        for gap in (BOTTOM_GAP * h, 2.0 * BOTTOM_GAP * h):
            cyl = Cylinder(draft=h - gap)
            diff = solve_diffraction(cyl, omega, res["orders"], res["n_modes"], res["p_max"])
            f.append(far_field_fx(cyl, omega, diff))
        ref = drift_convergence.maccamy_fuchs_fx(omega, Cylinder().radius, h)
        worst = max(worst, abs(2.0 * f[0] - f[1] - ref) / abs(ref))
        worst_raw = max(worst_raw, abs(f[0] - ref) / abs(ref))
        if abs(f[0] - ref) > LINEARITY_FLOOR * abs(ref):
            ratios.append(abs(f[1] - ref) / abs(f[0] - ref))
        else:
            skipped.append(omega)
    linear = bool(ratios) and all(LINEARITY_RATIO[0] <= q <= LINEARITY_RATIO[1] for q in ratios)
    return worst < 1e-3 and linear, (
        "extrapolated to g = 0: max |Fx - MF|/|MF| %.1e; at g = 1e-3 h itself %.1e; |F(2g) - MF|/|F(g) - MF| %s%s"
        % (worst, worst_raw, "%.3f ... %.3f" % (min(ratios), max(ratios)) if ratios else "not tested anywhere",
           ", not tested at %s rad/s (below the round-off floor)" % ", ".join("%.2f" % w for w in skipped)
           if skipped else ""))


def check_free_energy():
    """Free body, default cylinder, LIGHT_RESOLUTION, every WIDE_FREQUENCIES: the optical-theorem term of the total
    outgoing wave equals its energy flux, 2 pi sum_m abs(a_m)**2, to 1e-10 of the fixed body's flux; see the module
    docstring."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
    worst = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        diff = solve_diffraction(cyl, omega, res["orders"], res["n_modes"], res["p_max"])
        rads = {j: solve_radiation(cyl, omega, j, res["n_modes"], res["p_max"]) for j in MOTION_MODES}
        total = total_solution(diff, rads, motions(cyl, omega, diff, rads))
        a = far_field_amplitude(total)
        full = np.concatenate([a[:0:-1], a])
        flux = 2.0 * math.pi * np.sum(np.abs(full) ** 2)
        optical = -math.sqrt(8.0 * math.pi / total.k) * np.real(np.exp(0.75j * math.pi) * np.sum(full))
        fixed = 2.0 * math.pi * np.sum(np.abs(far_field_amplitude(diff)) ** 2 * np.where(np.arange(len(a)), 2.0, 1.0))
        worst = max(worst, abs(flux - optical) / fixed)
    return worst < 1e-10, "max |flux - optical|/(flux of the fixed body) %.1e, free body" % worst


BODY_MATRICES_SIDES = 4096


def check_body_matrices():
    """body_matrices against make_mesh.hydrostatics on a BODY_MATRICES_SIDES-sided polygon of the same radius, 1e-5;
    see the module docstring."""
    cyl = Cylinder()
    theta = 2.0 * math.pi * np.arange(BODY_MATRICES_SIDES) / BODY_MATRICES_SIDES
    polygon = [(cyl.radius * math.cos(t), cyl.radius * math.sin(t)) for t in theta]
    hydro = make_mesh.hydrostatics(polygon, cyl.draft)
    mass, restoring = body_matrices(cyl)
    worst = 0.0
    for mine, ref in ((mass, np.array(hydro["mass"])), (restoring, np.array(hydro["restoring"]))):
        nonzero = ref != 0.0
        if np.any((mine != 0.0) != nonzero):
            return False, "the zero pattern differs from make_mesh.hydrostatics"
        worst = max(worst, np.max(np.abs(mine[nonzero] - ref[nonzero]) / np.abs(ref[nonzero])))
    return worst < 1e-5, "max relative difference per nonzero entry %.1e (zero pattern equal)" % worst


def _first_order_case(cyl, omega, motion, n_modes, p_max, orders):
    """(total, xi) of the default problem at one resolution: diffraction, and for "free" the radiation solutions, the
    motions and the total potential."""
    diff = solve_diffraction(cyl, omega, orders, n_modes, p_max)
    if motion == "fixed":
        return diff, motions(cyl, omega, diff, None, "fixed")
    rads = {j: solve_radiation(cyl, omega, j, n_modes, p_max) for j in MOTION_MODES}
    xi = motions(cyl, omega, diff, rads)
    return total_solution(diff, rads, xi), xi


# Resolution of fx_identity: the three-family gap basis (20, 4, 4) and N = 64000, where the free-body error at 0.70 and
# 0.80 rad/s (4.0e-4 and 4.1e-4) has settled to about 3e-5 (module docstring, "Gap basis")
NEAR_CHECK_RESOLUTION = dict(n_modes=64000, p_max=(20, 4, 4), orders=15)


def _fx_pair(cyl, omega, motion, res=None):
    """(near F_x, far F_x) over rho g A**2 of the default problem at NEAR_CHECK_RESOLUTION."""
    res = NEAR_CHECK_RESOLUTION if res is None else res
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    total, xi = _first_order_case(cyl, omega, motion, res["n_modes"], res["p_max"], res["orders"])
    return near_field(cyl, omega, total, xi, quad)[1][0], far_field_fx(cyl, omega, total)


FX_FLOOR = 1e-4   # fx_identity: below FX_FLOOR times the largest abs(far F_x) of the motion, the error is over that floor


def check_fx_identity():
    """Near-field F_x (sum of the terms) against far-field F_x, fixed and free, every WIDE_FREQUENCIES, at
    NEAR_CHECK_RESOLUTION, to 5e-4 of max(abs(far F_x), FX_FLOOR scale), scale the largest abs(far F_x) of the motion
    (the floor of drift_convergence.percent; it acts only at free 0.40 rad/s, far F_x = 3e-6); see the module docstring."""
    cyl = Cylinder()
    worst, rows = 0.0, []
    for motion in ("fixed", "free"):
        pairs = [(omega, _fx_pair(cyl, omega, motion)) for omega in drift_convergence.WIDE_FREQUENCIES]
        scale = max(abs(far) for _, (_, far) in pairs)
        for omega, (near, far) in pairs:
            floored = abs(far) < FX_FLOOR * scale
            err = abs(near - far) / max(abs(far), FX_FLOOR * scale)
            worst = max(worst, err)
            rows.append("%s %.2f %.1e%s" % (motion, omega, err, " (floor)" if floored else ""))
    return worst < 5e-4, "max |near - far|/|far| %.1e; per case: %s" % (worst, ", ".join(rows))


def check_reference_point():
    """M_y about xr = (0, 0, -3) against M_y about the origin: M_y(xr) = M_y(0) - (xr x F)_y, fixed and free, at 0.70 and
    1.40 rad/s, to 1e-10 of max(abs(M_y(0)), abs((xr x F)_y)); see the module docstring."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    xr = np.array([0.0, 0.0, -3.0])
    worst = 0.0
    for motion in ("fixed", "free"):
        for omega in (0.70, 1.40):
            total, xi = _first_order_case(cyl, omega, motion, res["n_modes"], res["p_max"], res["orders"])
            f0 = near_field(cyl, omega, total, xi, quad, (0.0, 0.0, 0.0))[1]
            f1 = near_field(cyl, omega, total, xi, quad, tuple(xr))[1]
            shift = np.cross(xr, f0[:3])[1]
            worst = max(worst, abs(f1[4] - (f0[4] - shift)) / max(abs(f0[4]), abs(shift)),
                        abs(f1[0] - f0[0]) / abs(f0[0]))
    return worst < 1e-10, "max |M_y(xr) - M_y(0) + (xr x F)_y|/scale, and |F_x(xr) - F_x(0)|/|F_x| %.1e" % worst


def check_pair_symmetry():
    """T-_ij = conj(T-_ji) for omega_i = 0.70, omega_j = 0.95, free body, every term and component, to 1e-12 of the
    largest abs(T-_ij) of the term; see the module docstring."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    (si, xi_i), (sj, xi_j) = (_first_order_case(cyl, w, "free", res["n_modes"], res["p_max"], res["orders"])
                              for w in (0.70, 0.95))
    tij = near_field_terms(cyl, 0.70, si, xi_i, 0.95, sj, xi_j, quad)
    tji = near_field_terms(cyl, 0.95, sj, xi_j, 0.70, si, xi_i, quad)
    worst = 0.0
    for key in NEAR_FIELD_TERMS:
        scale = np.max(np.abs(tij[key]))
        if scale > 0.0:
            worst = max(worst, np.max(np.abs(tij[key] - np.conj(tji[key]))) / scale)
    return worst < 1e-12, "max |T_ij - conj(T_ji)|/max |T_ij| over the terms %.1e" % worst


LONG_WAVE_OMEGA = 0.05


def check_long_wave_motions():
    """Free body at omega = LONG_WAVE_OMEGA: heave follows the incident elevation, xi_3/eta_I(0) -> 1, and pitch the
    incident slope, xi_5/(-d eta_I/dx(0)) -> 1, within LONG_WAVE_TOLERANCE; see the module docstring."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
    omega = LONG_WAVE_OMEGA
    diff = solve_diffraction(cyl, omega, (0, 1), res["n_modes"], res["p_max"])
    rads = {j: solve_radiation(cyl, omega, j, res["n_modes"], res["p_max"]) for j in MOTION_MODES}
    xi = motions(cyl, omega, diff, rads)
    eta, slope = -1j, diff.k          # eta_I = -i e^(i k x): eta_I(0) = -i, d eta_I/dx(0) = k
    heave = abs(xi[3] / eta - 1.0)
    pitch = abs(xi[5] / -slope - 1.0)
    tol = long_wave_tolerance(cyl, omega, added_mass_damping(cyl, omega, rads)[0])
    ok = heave < tol[3] and pitch < tol[5]
    return ok, ("|xi_3/eta_I(0) - 1| %.2e (tolerance %.2e), |xi_5/(-d eta_I/dx(0)) - 1| %.2e (tolerance %.2e)"
                % (heave, tol[3], pitch, tol[5]))


CHECKS = {"roots": check_roots, "inner_products": check_inner_products, "bottom_limit_x1": check_bottom_limit_x1,
          "matching": check_matching, "scattering_energy": check_scattering_energy, "reciprocity": check_reciprocity,
          "radiation_energy": check_radiation_energy, "particular_bottom": check_particular_bottom,
          "bottom_limit_drift": check_bottom_limit_drift, "body_matrices": check_body_matrices,
          "free_energy": check_free_energy, "long_wave_motions": check_long_wave_motions,
          "fx_identity": check_fx_identity, "reference_point": check_reference_point,
          "pair_symmetry": check_pair_symmetry}


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
