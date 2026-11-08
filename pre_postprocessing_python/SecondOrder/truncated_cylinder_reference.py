"""Semi-analytical reference for a truncated vertical cylinder, by eigenfunction matching, to validate HAMS-MREL's mean drift.

Problem (the method of Garrett 1971 and Yeung 1981, with a corner basis; see Sources at the end, and the section
"Semi-analytical reference: sharp truncated cylinder" of docs/second_order_validation.md for the results). A vertical
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
(-d - z)**(-1/3). Reflecting the gap about the sea bed, where dphi/dz = 0, gives t = (z + h)/g in [-1, 1],
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
size P + 1 in c_p for the single family (a Galerkin method; for the three families of "Gap basis", in the kept Gram
combinations):

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
heading 0, with a(theta) the outgoing elevation amplitude of far_field_amplitude, eta_scattered ~ a(theta) e^(i k r)/sqrt(r)
as k r -> infinity, in exp(-i omega t) (Maruo 1960; Newman 1967). The general momentum form,
(1/2)(C_g/C_p)[X - int abs(a)**2 cos theta dtheta] with the optical term
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
in t, with their own p_max; an int P is the single lam = 1/6 family, unchanged. Each family's inner products are
Gegenbauer's integral with its lam (c_p = pi 2**(1 - lam) Gamma(2p + 2 lam)/((2p)! Gamma(lam))), and check_inner_products
tests all three against quad. The families are not orthogonal and together nearly dependent (the Gram matrix
int u_p u_q dt of (20, 4, 4) has 8 eigenvalues below 1e-10 of the largest; _gram, exact in sigma, 1 - abs(t) = sigma**3),
so _solve_order solves the Galerkin system in the orthonormal combinations of the eigenvectors above GRAM_CUTOFF = 1e-10
(_gram_transform), gap = T c; for m = 0 the zero-flux condition borders the system, since each family's p = 0 carries
flux. Without that reduction the coefficients reach 3e4 and the near field is erratic in N and P.

Matching tail (_matching_tail, MATRIX_TAIL). The matrix sums sum_n E E/(N kappa) - F F/(nu mu) of the corner family
converge only like N**(-4/3) (E_pn ~ n**(-2/3)), and the wall projection like N**(-5/3); the Hankel expansion gives
their non-oscillating leading tails in closed form (Hurwitz zeta), checked against brute-force sums from N = 2000 to 2e5
to 1e-3 of the tail (check_matching_tail: 1.1e-3 for the matrix, 1.8e-3 for the wall, N = 2000 to 64000). Adding them
makes the near-field identity converge like N**(-2) instead of N**(-1.3) and leaves the single-family results
N-converged at N = 4000 (fixed body, 0.70 rad/s, 6.21e-4 at N = 4000 and 6.25e-4 at N = 32000). It applies to every
solve, single family included, but only where the leading terms are valid (_tail_rows): the Hankel expansion needs
k_N g >> nu**2, nu = 2p + lam, and the omitted oscillating products are about 1/(k_N g) of the tail. At the default
cylinder that holds for every row at every resolution used here; at the 0.1 and 0.2 m gaps of the bottom-limit checks
(N = 1000, k_N g = pi and 2 pi) it holds for none, and those solves have no tail.

Convergence. The exterior and interior series of dphi/dr at r = a are expansions of a velocity that is singular at the
corner and jumps there, so their pointwise value converges slowly in N, about N**(-2/3) a few metres from the corner.
Without the matching tail phi, the excitation forces, the added mass and the damping converged about as N**(-1.6), and N
had to grow with P (N >= 100 P or so). With it they converge as N**(-5/3), for the single family and for (20, 4, 4)
alike: X_1, A_11 and B_11 at 0.70 rad/s change by 1.1e-8, 2.2e-8 and 1.1e-8 (relative) from N = 32000 to 64000,
measured exponent 1.67 (1.5 to 2 for X_5, A_55, B_55). N does not grow the linear system (basis.size unknowns, or the
fewer kept Gram combinations), only the inner products, so N in the tens of thousands is cheap: a diffraction solve
(16 orders) and the three radiation solves at N = 64000, (20, 4, 4), take about 3 s, most of it the inner products,
which _modes computes once per (omega, N, gap basis) for all of them.

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

Near field. near_field_terms writes the pressure-integration terms of the WAMIT v7.5 user manual, section 15.8,
eqs. (15.54) to (15.57), restated below as coded (the six terms docs/input_files.md lists for HAMS-MREL), for a pair of
first-order states i, j under the difference-frequency kernel P(f, g) = (f_i conj(g_j) + conj(f_j) g_i)/4 (for i = j,
Re(f conj(g))/2, the time average of the product of two harmonic quantities), over rho g, n the unit normal out of the
fluid (into the body), x relative to xr, xi' = xi + alpha x xr the translation of xr and alpha = (0, xi_5, 0) the
rotation (the module's xi are about the origin), X = xi' + alpha x (x - xr), phi_t = -i omega phi:

    I    (1/2) int_WL r' [P(zeta, zeta) - 2 P(zeta, X_3)] a dtheta, zeta = (i omega/g) phi(a, theta, 0), r' = n' or
         (x - xr) x n', n' = (-cos theta, -sin theta, 0) the unit horizontal normal into the body,
    Ia   zero: the wall is vertical at the waterline (n_z = 0),
    II   -(1/g) int_S r (1/2) sum_k P(d_k phi, d_k phi) dS,  r = n or (x - xr) x n,
    III  -(1/g) int_S r sum_k P(X_k, d_k phi_t) dS,
    IV   force P(alpha x G), moment P(alpha x H) + P(xi' x G), G = -(1/g) int n phi_t dS, H = -(1/g) int (x - xr) x n phi_t dS,
         the first-order dynamic pressure only,
    V    F_S,z and M_S,y of WAMIT eqs. (15.56) and (15.57), the formula HAMS-MREL uses, with A = pi a**2, (x_f, y_f) =
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

Self-checks: python -I truncated_cylinder_reference.py check [NAME ...] runs the named checks, or without NAME all of them
but uncertainty and compare_signs (SLOW_CHECKS: they are run by name, and by the reference subcommand).
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
  sums, the forces only the N modes); all converge in N, and the extrapolated residuals are 2e-9. The same at one
  frequency for the three-family basis (THREE_FAMILY_CASE: 0.70 rad/s, (20, 4, 4), N = 8000 and 16000): 3.7e-10
  extrapolated, 3.8e-7 and 7.5e-8 at N = 16000 itself. check_reciprocity's docstring says what each part tests.
- radiation_energy: same resolution and extrapolation: B11, B33 and B55 from the body pressure equal
  radiated_energy_damping, the energy flux at infinity, to 1e-6, and the same for THREE_FAMILY_CASE (3.7e-10
  extrapolated). It tests the flux formula and the sign and normalisation of B, and for heave and pitch the particular
  solutions. Neither radiation check bounds the error of A, B or X themselves (convergence does, for B11).
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
  wave. The same for THREE_FAMILY_CASE at N = 16000 (9e-16). It tests the sign and phase of X against the radiation
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
  the near field: it tests every term of F_x, the quadrature, and through the free body the motions, A and the
  radiation potentials near the body. Of the body velocities it tests what F_x sees: n_x is zero on the bottom, so F_x
  takes II and III from the wall only (the body condition and the wall series of phi and d phi/dz with their corner
  tails), and from the bottom only phi, through G_z in the rotation term IV. It does not test the bottom's tangential
  velocity d phi/dr (the interior series and its corner tail), which enters only F_z and M_y; corner_limit does. The
  rows are the reference values (DEFAULT_RESOLUTION, cached in _reference_row). Prescribed single
  motions at 0.70 rad/s, (20, 4, 4), N = 64000, against the momentum far field (X - int abs(a)**2 cos(theta) dtheta,
  valid without the lossless assumption): surge 1.8e-6, heave 5.6e-8, pitch -1.6e-4 absolute (F 1.92, 0.15, 0.80); the
  pitch residual falls slowly with the basis ((10, 1, 1): -2.0e-4) and is most of the free-body error at 0.70 rad/s.
  It does not test F_z or M_y; reference_point tests M_y's moment arms, corner_limit the bottom velocity.
- reference_point: xr = (0, 0, -3) against the origin, fixed and free, 0.70 and 1.40 rad/s, LIGHT_RESOLUTION:
  M_y(xr) = M_y(0) - (xr x F)_y and F_x unchanged, to 1e-10. The motions are re-expressed about xr (xi' = xi + alpha x xr),
  so it tests that every moment arm, the IV moment (alpha x H + xi' x G) and M_S transform as a rigid shift. F_z is not
  compared: F_S,z depends on Z_o (Near field).
- pair_symmetry: free body, omega_i = 0.70, omega_j = 0.95: T-_ij = conj(T-_ji) for every term and component, and
  T+_ij = T+_ji under the sum-frequency kernel _pair_sum, (f_i g_j + f_j g_i)/4, each to 1e-12. It tests that each
  state's frequency, motion and conjugation enter the kernels symmetrically, not the values; T+ is otherwise unused
  here (only i = j, the mean drift, is computed and checked).
- polylog: _polylog against direct sums, 1e-10, for s = 1 and the orders of the corner tails, real mu from -0.05 to -60
  and complex mu off the unit circle; and that it raises for Re mu > 0. Measured 1.4e-14. Its docstring says what it
  does not cover (the unit circle itself).
- matching_tail: at 0.70 rad/s, orders 0 and 1, (20, 4, 4), the closed-form tail of _matching_tail summed from N = 2000
  to 64000 (the tail at 2000 less the tail at 64000) equals the brute-force sums of the matrix and of the wall projection
  over those modes to 5e-3 of the largest (measured 1.1e-3 and 1.8e-3: what is left is the O(1/n) part the model omits);
  _tail_rows keeps every row at the default cylinder at LIGHT, CHECK, the matching check's and the default resolutions,
  and none at a 0.1 m gap, N = 1000, P = 10, where the tail is zero. It covers the tail's coefficients, phases and
  powers, and where it applies, not its effect on the solution (the N convergence of the other checks shows that).
- corner_limit: fixed and free, 0.70 and 1.40 rad/s, orders 0 ... 3, CORNER_LIMIT_RESOLUTION ((20, 4, 4), N = 16000):
  rho**(1/3) d phi_m/dz on the wall and rho**(1/3) d phi_m/dr on the bottom, from _body_series at rho = 1e-3, 1e-4,
  1e-5 m, each fitted with c0 + c1 s + c2 s**2 in s = rho**(1/3) (the wedge expansion), have the same corner limit c0,
  to 1e-3 relative: the 270 degree wedge's leading term rho**(2/3) cos(2 alpha/3) gives the same tangential velocity
  rho**(-1/3) on both faces, in this sign convention. Measured 1e-4 to 1.9e-4 (4e-6 to 1e-4 at N = 64000, 2.7e-3 at
  N = 4000); the raw values at rho = 1e-3 differ by up to 0.12 (free, m = 0), from the c1 s term, which differs between
  the faces, so the fit is needed. Injected errors: the bottom corner tail dropped gives 1.0,
  scaled by 1.01 gives 1.1e-2 (both caught); dropping the lam = 5/6 bottom tail changes the limit by 3e-5 only (it is
  of order rho**(1/3) in rho**(1/3) u_t), so that tail is not tested. It covers the bottom series d phi/dr, its corner
  tail and the particular solution's velocity, which only F_z and M_y see; it does not test them away from the corner.
- convergence and uncertainty: see "Reference values".
- compare_fixture: compare on a temporary RESULTS written in HAMS's formats from LIGHT_RESOLUTION reference rows at 0.40,
  0.70 and 1.60 rad/s, passed through the reference CSV: every difference 0 (to 1e-9) for a case with all files, heading
  30 rows ignored; "missing" terms and control surface for a case without MeanDriftTerms.txt and
  MeanDriftControlSurface.7; the 1.60 rows "irregular" and no other; free F_x at 0.40 rad/s empty (below the floor); an
  a_sharp_ folder that is not a study case and one without output/Wamit_format skipped; and irregular_frequencies equal
  to the 1.548 and 1.939 of drift_convergence.IRREGULAR_FREQUENCIES. Every complex value is written with the phase
  FIXTURE_PHASE = 0.6 rad, Re = MOD cos(phase), Im = MOD sin(phase), so a reader that takes the wrong column fails
  (measured: MOD instead of Re for the drift files, Re instead of MOD for ExcForce.3, Im instead of Re for
  MeanDriftTerms.txt, each FAIL). The LIGHT_RESOLUTION rows carry no stated uncertainty (every *_unc nan), and compare's
  q_unc comes from the row: nan everywhere except FIXTURE_UNCERTAINTY, set on one row before compare. It tests the
  readers, the indices, the floor, the flags and the uncertainty passing, not HAMS.
- compare_signs (in SLOW_CHECKS): see "Comparison with HAMS". Without the study's results it prints SKIP and counts as
  neither passed nor failed.

Reference values. python -I truncated_cylinder_reference.py reference [--out CSV] writes, at DEFAULT_RESOLUTION
(N = 64000, gap basis (20, 4, 4), orders 0 ... 15) and every WIDE_FREQUENCIES, the fixed and the free rows of reference
(A, B, abs(X), abs(xi), F_x far and near, F_z, M_y and the near-field terms, in HAMS's normalisations, listed in
reference's docstring), then the convergence table, then runs every check; exit 1 if one fails.

- convergence (the check): for F_x far, F_z, M_y (near field, about the origin) and B11 at 0.70, 0.95 and 1.40 rad/s
  (CONVERGENCE_FREQUENCIES), fixed and free, the values at the three resolutions of CONVERGENCE_LADDER, N = 32000 with
  (10, 2, 2), N = 64000 with (20, 4, 4) (the default) and N = 128000 with (30, 6, 6), and at the default with orders
  0 ... 30, are finite, and the default's spread against the finest and against 2M is below 1e-4 of its value
  (CONVERGENCE_SPREAD; for a near-zero value, of CONVERGENCE_FLOOR = 1e-4 of the largest of its quantity and motion,
  which applies nowhere here). One row exceeds 1e-4 and is a named known exceedance (KNOWN_UNCERTAINTY): free F_z at
  0.70 rad/s, -0.6128 from II -29.94, III +42.58, IV -13.25, spread 2.6e-4 (1.6e-4 absolute, all of it in II on the
  bottom); it passes against its stated uncertainty, 3e-3, and would fail if the spread grew past that. Every other
  spread is below 1.6e-5 (F_x far below 1.2e-7, B11 below 4e-9); the 2M spread is zero to round-off everywhere (the
  far-field amplitude of order 16 is 1e-25 of order 0 at 1.40 rad/s).
- DEFAULT_RESOLUTION is the coarsest that meets both gates: N = 32000 with (10, 2, 2) as the default would fail the
  spread against N = 64000 with (20, 4, 4) (fixed F_z at 1.40 rad/s, 1.3e-4), and the free body needs the three
  families and N = 64000 for fx_identity (Gap basis; at N = 32000 (10, 2, 2) passes it at 0.70 rad/s only because it is
  not converged in N: -3.4e-4, and -5.1e-4 at N = 128000).
- Scope. DEFAULT_RESOLUTION and the stated uncertainties below are measured for the default hull only (Cylinder():
  a 10, d 10, h 100, h/a = 10) at the frequencies of WIDE_FREQUENCIES. The N needed grows with h/a (measured, free F_z at
  0.90 rad/s, gap basis (20, 4, 4), N = 4000 against 16000: 4.9 % for a 10, d 10, h 500; 0.8 % for a 5, d 8, h 40), so
  uncertainty() returns nan (not measured) for any other cylinder, resolution or frequency (omega matched to
  WIDE_FREQUENCIES to 1e-9 relative), and a new geometry needs the convergence and fx_identity measurements repeated
  for it (convergence_table takes the cylinder; fx_identity's near against far F_x) before its values carry one.
- Uncertainty, stated per row (uncertainty(): KNOWN_UNCERTAINTY, else DEFAULT_UNCERTAINTY = 1e-4 relative, within the
  scope above; the CSV columns Fx_far_unc, Fx_near_unc, Fz_unc, My_unc). Evidence for the default 1e-4:
  - Fixed body: spread of the default against N = 128000 with (30, 6, 6), measured at every WIDE_FREQUENCIES for F_z
    and M_y (outside the check, once): F_z 6e-7 to 6.6e-6 up to 1.20 rad/s, 1.6e-5 at 1.40, 3.1e-5 at 1.60, 4.8e-5 at
    1.80 (F_z falls to -0.025 there); M_y below 1.5e-6. fx_identity below 2.4e-6. F_x far and B11 far below.
  - Free body away from 0.60 to 0.80 rad/s: the convergence spreads at 0.95 and 1.40 (below 4e-6), fx_identity below
    2e-4 (7e-6 at 0.95 rad/s). Not measured for F_z and M_y at 0.40, 1.00 ... 1.20, 1.60, 1.80 beyond fx_identity.
  - Free body at 0.60 to 0.80 rad/s: near-field F_x settles about 4.2e-4 of far F_x from it as N grows, systematically,
    pitch only (a rho**2 log rho corner term, which no power family represents, is a plausible cause, Gap basis), and
    does not fall with the gap basis: Fx_near is stated at 5e-4 there (far F_x changes by less than 1.2e-7 over the
    resolutions of convergence; its stated uncertainty stays 1e-4). F_z and M_y are not
    covered by that identity, and their spreads understate their error. check_uncertainty (drift_uncertainty) fits F_z
    and M_y against the F_x residual r over the gap bases (10, 2, 2), (20, 4, 4), (30, 6, 6) at N = 128000:
      r = -2.56e-4, -1.68e-4, -2.09e-4 (0.60 rad/s); -5.11e-4, -4.17e-4, -4.46e-4 (0.70); +4.94e-4, +4.24e-4,
      +4.26e-4 (0.80);
    slopes dF_z/dr = -0.47, -3.8, +24 and dM_y/dr = 0.03, 0.82, 114; the linear fit's largest residual is 2e-7, 3e-7,
    8e-5 for F_z. The estimated error is the distance of the default value to the fit's value at r = 0: F_z 2.4e-5,
    2.5e-3, 3.3e-3 relative, M_y 3.0e-5, 4.1e-5, 3.2e-5 relative. So KNOWN_UNCERTAINTY states F_z 3e-3 at 0.70 and 4e-3 at
    0.80 rad/s; the rest stays within 1e-4. The lever arm is short: the residuals span only 0.9e-4 (0.70 rad/s) and
    0.7e-4 (0.80) and the extrapolation to r = 0 is 5 to 6 times that span, and it assumes the error vanishes where
    the residual does, which nothing here tests (the residual does not fall with the basis). It is an estimate of
    the size, not a bound.
  - Free body at 1.40 to 1.80 rad/s: near-field F_x is -1.11e-4, -1.29e-4, -1.26e-4 from far F_x (fx_identity), so
    KNOWN_UNCERTAINTY states 2e-4 for it there. F_z and M_y there have no separate estimate (check_uncertainty covers
    0.60 to 0.80 rad/s only); they carry DEFAULT_UNCERTAINTY.
- Run time, laptop, measured once each (November 2026): check (the 20 checks without SLOW_CHECKS) 5 min 4 s, most of
  it convergence (146 s) and fx_identity (93 s, the 22 reference rows). reference --out 7 min 11 s for all 22 checks: the
  rows 92 s, convergence 135 s, uncertainty 137 s (the rows at N = 128000 that no other check computes); compare_signs
  takes 0 s there, since the rows are cached, and alone about the 92 s of the rows (an estimate, not measured). compare
  93 s computing the rows, under 1 s with --reference. A row pair (fixed and free) at the default costs about
  6 s: they share the solves (_solved), the inner products and radial ratios (_modes) and the solution-independent
  part of the body series (_series_context).

Comparison with HAMS. python -I truncated_cylinder_reference.py compare RESULTS [--out CSV] [--reference CSV] reads every
folder of RESULTS whose name starts with a_sharp_ (the sharp hull of group A of drift_convergence.study_cases: fixed and
free, differencing and source, multiples 1 to 4, 448 to 7168 panels), output/Wamit_format/, heading 0 only, and prints per
case and omega the percentage difference (HAMS - reference)/abs(reference) of A and B (11, 15, 33, 55; 13 and 35 are zero
by symmetry), abs(X), abs(RAO) (both motions against the free body: HAMS writes Motion.4 for a fixed body too), F_x near
and far, F_z, M_y, each near-field term, and the control-surface F_x of MeanDriftControlSurface.7 against the far-field
reference, "missing" where a file is absent. Only moduli and real drift values are compared, so the factor i between
HAMS's phases and this module's (eta_I(0) = 1 against -i) plays no part. Each drift percentage is followed by the
stated uncertainty the reference row carries (its *_unc column, computed or read from the CSV; [n/m], not measured,
where it is nan: rows that are not the default hull at DEFAULT_RESOLUTION), and the header prints the rows'
resolution. A HAMS frequency is matched to a reference row within a relative FREQUENCY_MATCH = 1e-6, else the row is
flagged "no reference"; the files' frequencies are matched to each other the same way. Each block is headed by the
panels and the area error of the waterline, a ntheta-sided inscribed polygon (polygon_area_error: -6.4e-3 at 32 sides, -4.0e-4 at 128). A percentage is empty where
abs(reference) is below drift_convergence.percent's floor, 1e-4 of the largest abs(reference) of the quantity and motion
over the reference rows (for a term, of the term and its total). Rows within 5 % of an irregular frequency of the hull
(irregular_frequencies: 1.5483 and 1.9394 rad/s, the only two below the top of the scan, which compare asserts; the 1.60
rows) are marked "irregular", not counted as HAMS errors. An a_sharp_ folder that is not a case of
drift_convergence.study_cases, or that has no output/Wamit_format, is skipped with a printed note. The
reference rows are computed at DEFAULT_RESOLUTION (cached; 93 s) or read from a CSV of the reference subcommand (under 1 s).
compare_signs, at the finest source mesh (7168 panels), tests the sign of every near-field term and of far-field F_x
above the floor, and of the near-field totals where they are at least SIGN_CANCELLATION = 5 % of their largest term: the
free F_z at 0.70 and 0.80 rad/s (1.4 and 1.0 % of their terms) and the free F_x at 0.60 and 0.70 rad/s are near
cancellations whose sign HAMS's term errors of a few % can flip. Measured once against the study's results (not
reproducible from the repository): every tested sign agrees (210 values); HAMS's free F_z at 0.70 rad/s is +0.059 against
-0.613, and its Richardson extrapolation from 448, 1792 and 7168 panels (source) is -0.65. compare_signs is in
SLOW_CHECKS: it reads data outside the repository and computes the 22 reference rows.

Sources.
- Garrett, C.J.R. (1971). Wave forces on a circular dock. J. Fluid Mech. 46(1), 129-139. doi:10.1017/S0022112071000430.
  Eigenfunction matching for a truncated cylinder (diffraction).
- Yeung, R.W. (1981). Added mass and damping of a vertical cylinder in finite-depth waters. Appl. Ocean Res. 3(3),
  119-133. doi:10.1016/0141-1187(81)90101-2. The same matching for the radiation problems.
- Mavrakos, S.A. (1988). The vertical drift force and pitch moment on axisymmetric bodies in regular waves. Appl. Ocean
  Res. 10(4), 207-218. doi:10.1016/S0141-1187(88)80005-1. Mean F_z and M_y of axisymmetric bodies by matched eigenfunctions.
- MacCamy, R.C. and Fuchs, R.A. (1954). Wave forces on piles: a diffraction theory. Technical Memorandum 69, Beach
  Erosion Board, US Army Corps of Engineers. The bottom-mounted limit.
- Maruo, H. (1960). The drift of a body floating on waves. J. Ship Res. 4, 1-10. Far-field F_x.
- Newman, J.N. (1967). The drift force and moment on ships in waves. J. Ship Res. 11(1), 51-60.
  doi:10.5957/jsr.1967.11.1.51. Far-field drift force and moment.
- WAMIT, Inc. WAMIT User Manual, Version 7.5, section 15.8, eqs. (15.54) to (15.57): mean drift by pressure integration.
- DLMF 18.17.17 (Gegenbauer's integral) and 25.12.12 (polylogarithm), https://dlmf.nist.gov.
"""
import argparse
import csv
import functools
import math
import os
import sys
import time
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


@dataclass(frozen=True)
class Cylinder:
    """Truncated vertical cylinder, lengths in m. Defaults: the study's sharp hull (HULLS["sharp"] in drift_convergence.py).
    Frozen, so that it can key the caches of the solves (_modes, _solved)."""
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


@dataclass(eq=False)
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


@functools.lru_cache(maxsize=None)
def _gram_transform(basis, cutoff=None):
    """T, size x K: the columns span the basis with orthonormal combinations, T^T G T = I, keeping the eigenvalues of G
    above cutoff times the largest; and the kept eigenvalue range (for the report). Cached per basis (gegenbauer_basis
    returns one object per spec)."""
    cutoff = GRAM_CUTOFF if cutoff is None else cutoff
    values, vectors = np.linalg.eigh(_gram(basis))
    keep = values > cutoff * values[-1]
    return vectors[:, keep] / np.sqrt(values[keep]), (values[keep][0] / values[-1], int(np.sum(~keep)))


def basis_spec(p_max):
    """The p_max of each family in FAMILIES: an int P is (P, -1, -1), the single lam = 1/6 family; a tuple gives the
    three (-1: the family is absent)."""
    if isinstance(p_max, (int, np.integer)):
        return (int(p_max), -1, -1)
    spec = tuple(int(p) for p in p_max)
    if len(spec) != len(FAMILIES) or spec[0] < 0:
        raise ValueError("p_max: an int or a tuple of %d ints with the first >= 0, got %r" % (len(FAMILIES), p_max))
    return spec


def _spec_key(p_max):
    """p_max as a hashable cache key: an int, or a tuple of ints."""
    return int(p_max) if isinstance(p_max, (int, np.integer)) else tuple(int(p) for p in p_max)


def gegenbauer_basis(p_max, nodes=GAUSS_JACOBI_NODES):
    """The gap basis for p_max (basis_spec), each family with a Gauss-Jacobi rule of `nodes` points for
    alpha = beta = lam - 1/2. One object per (p_max, nodes): the arrays are not to be changed."""
    return _gegenbauer_basis(_spec_key(p_max), nodes)


@functools.lru_cache(maxsize=None)
def _gegenbauer_basis(p_max, nodes):
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


def _radial_ratios(cyl, m, k, kn):
    """(kappa, mu), n = 0 ... len(kn): the radial log derivatives at r = a of the exterior modes (H_m'/H_m for n = 0,
    K_m'/K_m) and of the interior modes (m/a for n = 0, I_m'/I_m), by scaled Bessel functions."""
    a, g = cyl.radius, cyl.gap
    kappa = np.empty(len(kn) + 1, complex)
    kappa[0] = k * h1vp(m, k * a) / hankel1(m, k * a)
    x = kn * a
    kappa[1:] = -0.5 * kn * (kve(m - 1, x) + kve(m + 1, x)) / kve(m, x)
    lam = np.arange(1, len(kn) + 1) * math.pi / g
    mu = np.empty(len(kn) + 1)
    mu[0] = m / a
    y = lam * a
    mu[1:] = 0.5 * lam * (ive(m - 1, y) + ive(m + 1, y)) / ive(m, y)
    return kappa, mu


@dataclass(eq=False)
class _Modes:
    """The order-independent part of the matching at one (cylinder, omega, N, p_max): the roots, the basis, the inner
    products, the mode norms and the wall moments, and per order the radial ratios (computed once per order). Shared by
    solve_diffraction, solve_radiation and _body_series through the cache of _modes."""
    cyl: Cylinder
    omega: float
    k: float
    kn: np.ndarray
    basis: Basis
    inner: tuple
    norms: np.ndarray
    moments: tuple
    _ratios: dict

    def ratios(self, m):
        if m not in self._ratios:
            self._ratios[m] = _radial_ratios(self.cyl, m, self.k, self.kn)
        return self._ratios[m]


@functools.lru_cache(maxsize=2)
def _interior_cached(cyl, n_modes, p_max):
    """inner_interior, which does not depend on omega, cached for the frequency loops."""
    return inner_interior(cyl, n_modes, gegenbauer_basis(p_max))


@functools.lru_cache(maxsize=2)
def _modes_cached(cyl, omega, n_modes, p_max):
    k, kn = wave_numbers(omega, cyl.depth, n_modes)
    basis = gegenbauer_basis(p_max)
    inner = (inner_exterior(cyl, k, kn, basis), _interior_cached(cyl, n_modes, p_max))
    return _Modes(cyl, omega, k, kn, basis, inner, _mode_norms(k, kn, cyl.depth), _wall_moments(cyl, k, kn), {})


def _modes(cyl, omega, n_modes, p_max):
    """_Modes of (cyl, omega, N, p_max), cached (the last two): the fixed and the free body at one omega, and the
    diffraction and radiation solves, share one."""
    return _modes_cached(cyl, float(omega), int(n_modes), _spec_key(p_max))


MATRIX_TAIL = True   # add the sums over n > N of the leading asymptotic terms to the matching system (_matching_tail)
TAIL_ORDER_RATIO = 0.25     # _tail_rows: the tail applies to the u_p with nu_p**2/(2 x_N) <= this, nu_p = 2p + lam, and
TAIL_MIN_ARGUMENT = 100.0   # only where x_N = (N + 1) pi g/h, the smallest Bessel argument of the tail, is at least this


def _tail_rows(cyl, basis, n_modes):
    """Boolean per basis function: where the leading terms of _matching_tail are valid. The Hankel expansion of
    J_nu(x) needs x >> nu**2 (its first correction is (4 nu**2 - 1)/(8 x)), and the oscillating products the tail
    omits are about 1/x of it; x_N = (N + 1) pi g/h is the smallest argument k_n g of the tail. Measured against
    brute-force sums (check_matching_tail): 1e-3 of the tail at the default gap, N = 2000, (20, 4, 4), where every row
    has nu**2/(2 x_N) <= 0.14; at a 0.1 m gap, N = 1000, P = 10, x_N = pi and the leading terms are 10 to 16 % off
    even for p = 0, and wrong in order of magnitude for nu > x_N. So the tail applies to the rows and columns with
    nu_p**2/(2 x_N) <= TAIL_ORDER_RATIO, and nowhere if x_N < TAIL_MIN_ARGUMENT: at every resolution of the default
    cylinder used here it applies everywhere, at the bottom-limit gaps (0.1 and 0.2 m, N = 1000) nowhere."""
    x = (n_modes + 1) * math.pi * cyl.gap / cyl.depth
    nu = np.concatenate([2.0 * np.arange(f.p_max + 1) + f.lam for f in basis.families])
    if x < TAIL_MIN_ARGUMENT:
        return np.zeros(nu.shape, bool)
    return nu ** 2 / (2.0 * x) <= TAIL_ORDER_RATIO


def _matching_tail(cyl, basis, n_modes):
    """(T, w): the sums over n > N of the non-oscillating leading terms of the matching system, by the Hurwitz zeta
    function. From the Hankel expansion, E_pn ~ (g/2) sqrt(2/pi) c_p (k_n g)**(-beta_p) cos(k_n g - chi_p)/cos(k_n h)
    and F_pn ~ (g/2) sqrt(2/pi) c_p (n pi)**(-beta_p) (-1)**n cos(chi_p), beta = lam + 1/2, chi = (lam/2 + 1/4) pi, with
    N_n kappa_n ~ -h k_n/(2 cos(k_n h)**2), nu_n mu_n ~ n pi/2, k_n ~ n pi/h, and the wall projection of a wall velocity
    V at the corner ~ -V sin(k_n g)/(k_n cos(k_n h)). T[p, q] is the tail of sum E E/(N kappa) - sum F F/(nu mu), w[p]
    that of sum E wall/(N kappa) per unit V. What is left falls one power of n faster (the O(1/n) terms of the expansions
    and the oscillating products, cos(2 k_n g - ...)). Rows, columns and entries of w where they are not valid
    (_tail_rows) are zero."""
    g, h = cyl.gap, cyl.depth
    lam = np.concatenate([np.full(f.p_max + 1, f.lam) for f in basis.families])
    c = np.concatenate([_gegenbauer_constant(np.arange(f.p_max + 1), f.lam) for f in basis.families])
    beta, chi = lam + 0.5, (0.5 * lam + 0.25) * math.pi
    s = beta[:, None] + beta[None, :] + 1.0
    hurwitz = zeta(s, n_modes + 1.0)
    cc = np.outer(c, c)
    ext = (-(g ** 2 / (2.0 * math.pi * h)) * cc * np.cos(chi[:, None] - chi[None, :])
           * g ** (1.0 - s) * (h / math.pi) ** s * hurwitz)
    intr = (g ** 2 / math.pi) * cc * np.outer(np.cos(chi), np.cos(chi)) * math.pi ** (-s) * hurwitz
    w = ((g / (2.0 * h)) * math.sqrt(2.0 / math.pi) * c * np.sin(chi) * g ** (-beta) * (h / math.pi) ** (beta + 2.0)
         * zeta(beta + 2.0, n_modes + 1.0))
    valid = _tail_rows(cyl, basis, n_modes)
    return (ext - intr) * np.outer(valid, valid), w * valid


def _solve_order(cyl, m, k, kn, basis, incident=0.0, wall=None, particular=None, inner=None, wall_velocity=(0.0, 0.0),
                 modes=None):
    """One order of the matching; see the module docstring. A complex system of the basis.size unknowns, or of the kept
    Gram combinations with several families (plus one for m = 0). inner: _inner_products, if already computed; modes: the
    _Modes of (cyl, omega, N, p_max), which supplies inner, the norms and the radial ratios."""
    a, g = cyl.radius, cyl.gap
    n_modes = len(kn)
    if modes is not None:
        big_e, big_f = modes.inner
        norms = modes.norms
        kappa, mu = modes.ratios(m)
    else:
        big_e, big_f = _inner_products(cyl, k, kn, basis) if inner is None else inner   # (size, N+1) each
        norms = _mode_norms(k, kn, cyl.depth)
        kappa, mu = _radial_ratios(cyl, m, k, kn)            # Z_n' radial log derivatives at r = a
    ka = k * a
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
    modes = _modes(cyl, omega, n_modes, p_max)
    k, kn = modes.k, modes.kn
    amplitude = -make_mesh.G / omega                       # phi_I = -(g A/omega) Z_0 e^(i k x), A = 1
    solved = {}
    for m in _order_list(orders):
        inc = amplitude * (1.0 if m == 0 else 2.0) * 1j ** m   # e^(ikx) = sum eps_m i^m J_m(kr) cos(m theta)
        solved[m] = _solve_order(cyl, m, k, kn, modes.basis, incident=inc, modes=modes)
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
    modes = _modes(cyl, omega, n_modes, p_max)
    k, kn, basis = modes.k, modes.kn, modes.basis
    m = RADIATION_ORDER[mode]
    i0, i1 = modes.moments
    wall = {1: i0, 3: None, 5: i1}[mode]
    particular = None if mode == 1 else _Particular(cyl, mode, basis, n_modes)
    velocity = {1: (1.0, 0.0), 3: (0.0, 0.0), 5: (0.0, 1.0)}[mode]   # dphi/dr on the wall, as `wall` above
    order = _solve_order(cyl, m, k, kn, basis, wall=wall, particular=particular, wall_velocity=velocity, modes=modes)
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
# rho**(-2/3) (module docstring, "Quadrature", for why the default is not Gauss-Jacobi)
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


@dataclass(eq=False)
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


def _corner_rule(length, n, rule=None):
    """(rho, W): int_0^length f(rho) drho ~ sum(W f(rho)), n nodes. rule (default CORNER_RULE): "cube", Gauss-Legendre
    in s = (rho/length)**(1/3), rho = length s**3, W = 3 length s**2 w, exact for f = rho**(-2/3) times a polynomial in
    rho**(1/3) of degree < 2 n (module docstring, "Quadrature"); or "jacobi", Gauss-Jacobi for the weight rho**(-2/3),
    exact for f = rho**(-2/3) times a polynomial in rho of degree < 2 n (scipy's weight (1 - x)**0 (1 + x)**beta,
    rho = length (1 + x)/2)."""
    rule = CORNER_RULE if rule is None else rule
    if rule == "cube":
        s, w = _legendre_rule(0.0, 1.0, n)
        return length * s ** 3, 3.0 * length * s ** 2 * w
    x, w = roots_jacobi(n, 0.0, CORNER_EXPONENT)
    return 0.5 * length * (1.0 + x), 0.5 * length * w * (1.0 + x) ** (-CORNER_EXPONENT)


def _legendre_rule(lo, hi, n):
    """Gauss-Legendre nodes and weights on [lo, hi]."""
    x, w = np.polynomial.legendre.leggauss(n)
    return lo + 0.5 * (hi - lo) * (1.0 + x), 0.5 * (hi - lo) * w


def body_quadrature(cyl, n_corner, n_far, n_theta=NEAR_THETA):
    """Quadrature of the wetted surface (_body_quadrature), one object per arguments and CORNER_RULE, so that the
    cache of _series_context sees the same quadrature."""
    return _body_quadrature(cyl, int(n_corner), int(n_far), int(n_theta), CORNER_RULE)


@functools.lru_cache(maxsize=None)
def _body_quadrature(cyl, n_corner, n_far, n_theta, rule):
    """Quadrature of the wetted surface. On each face, the corner rule `rule` (_corner_rule; the default "cube",
    Gauss-Legendre in (rho/L)**(1/3)) in the distance rho from the corner, n_corner nodes, within min(gap, d) of it on the
    wall and min(gap, a/2) on the bottom, and Gauss-Legendre in rho with n_far nodes beyond (on the default cylinder the
    corner zone is the whole wall, d = 10 < g = 90, so the wall has no Gauss-Legendre part)."""
    a, d, g = cyl.radius, cyl.draft, cyl.gap
    faces = []
    for name, extent, zone in (("wall", d, min(g, d)), ("bottom", a, min(g, 0.5 * a))):
        rho, w = _corner_rule(zone, n_corner, rule)
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


def _extended_roots(omega, h, kn, count):
    """k_n, n = 1 ... count: the solution's kn, then y = n pi - atan(nu h/y), y = k_n h, by fixed-point iteration (its
    derivative nu h/(y**2 + (nu h)**2) is below 1e-6 for the n beyond any N used here, so eight steps reach round-off)."""
    c = omega ** 2 / make_mesh.G * h
    n = np.arange(len(kn) + 1, count + 1)
    y = (n - 0.5) * math.pi
    for _ in range(8):
        y = n * math.pi - np.arctan(c / y)
    return np.concatenate([kn, y / h])


def _smooth_window(count, start):
    """w_n, n = 1 ... count: 1 up to n = start count, then a C-infinity step down to 0 at n = count."""
    x = np.clip((np.arange(1, count + 1) - start * count) / ((1.0 - start) * count), 0.0, 1.0)
    inside = (x > 0.0) & (x < 1.0)
    e0 = np.where(inside, np.exp(-1.0 / np.where(inside, x, 1.0)), 0.0)
    e1 = np.where(inside, np.exp(-1.0 / np.where(inside, 1.0 - x, 1.0)), 0.0)
    return np.where(x <= 0.0, 1.0, np.where(x >= 1.0, 0.0, e1 / np.where(inside, e0 + e1, 1.0)))


POLYLOG_TERMS = 60    # terms of the zeta series of _polylog; they fall like (|mu|/(2 pi))**k, below 0.53**k where used
POLYLOG_DIRECT = 40   # terms of the direct sum of _polylog where Re mu < -1: the remainder is below e**(-40) of the sum
POLYLOG_RADIUS = 3.5  # the zeta series is used for |mu| <= this (with Re mu >= -1, |Im mu| <= pi: |mu| <= 3.3)


def _polylog(order, mu):
    """Li_s(e^mu) = sum_(n >= 1) n**(-s) e^(n mu), Re mu <= 0, mu != 0, elementwise. Where Re mu < -1, any s: the
    direct sum (POLYLOG_DIRECT terms). Elsewhere, s = 1: -log(1 - e^mu) (by expm1); s not an integer: the zeta series
    Li_s(e^mu) = Gamma(1 - s) (-mu)**(s - 1) + sum_k zeta(s - k) mu**k/k! (DLMF 25.12.12), which converges only for
    |mu| < 2 pi and is used for |mu| <= POLYLOG_RADIUS (ValueError beyond, and for Re mu > 0). Without the direct branch
    the series diverged silently: Li_(2/3)(e**(-12)) gave -6.4e9 for 6.1e-6 (the bottom tails reach
    mu = -pi a/(2 g), -12 at a 1.3 m gap). check_polylog tests both branches."""
    mu = np.asarray(mu, complex)
    if np.any(mu.real > 0.0):
        raise ValueError("_polylog: Re mu > 0, outside the domain of the series")
    direct = mu.real < -1.0
    near = np.where(direct, -0.5 + 0j, mu)                     # a harmless point where the direct sum is used
    if order == 1:
        out = -np.log(-np.expm1(near))                         # expm1 keeps 1 - e^mu accurate near mu = 0
    else:
        if np.any(np.abs(mu[~direct]) > POLYLOG_RADIUS):
            raise ValueError("_polylog: |mu| > %g with Re mu >= -1, where the zeta series of %d terms is not accurate"
                             % (POLYLOG_RADIUS, POLYLOG_TERMS))
        out = gamma(1.0 - order) * (-near) ** (order - 1.0)
        term = np.ones_like(near)
        for k in range(POLYLOG_TERMS):
            out = out + zeta(order - k) * term
            term = term * near / (k + 1)
    if np.any(direct):
        n = np.arange(1, POLYLOG_DIRECT + 1).reshape((-1,) + (1,) * mu.ndim)
        far = np.where(direct, mu, -2.0 + 0j)
        out = np.where(direct, np.sum(n ** (-order) * np.exp(n * far), axis=0), out)
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


@dataclass(eq=False)
class _SeriesContext:
    """The solution-independent part of _body_series at one (cylinder, omega, N, p_max, count, quadrature, orders): the
    extended roots, inner products, norms and wall moments, the window, the corner models and the Bessel functions on
    the bottom nodes. The fixed and the free body at one omega share it (_series_context)."""
    k: float
    kn: np.ndarray
    basis: Basis
    big_e: np.ndarray
    big_f: np.ndarray
    norms: np.ndarray
    moments: tuple
    window: np.ndarray
    models: dict
    cos_wall: np.ndarray
    sin_wall: np.ndarray
    z0: np.ndarray
    dz0: np.ndarray
    lam: np.ndarray
    sign: np.ndarray
    keep: np.ndarray
    bessel_a: list
    bessel_r: list
    growth: list
    modes: object                 # the _Modes of the solve when count = N (its radial ratios are shared), else None
    _ratios: dict

    def ratios(self, cyl, m):
        """(kappa_n, mu_n), n = 1 ... count."""
        if self.modes is not None:
            kappa, mu = self.modes.ratios(m)
            return kappa[1:], mu[1:]
        if m not in self._ratios:
            kappa, mu = _radial_ratios(cyl, m, self.k, self.kn)
            self._ratios[m] = (kappa[1:], mu[1:])
        return self._ratios[m]


@functools.lru_cache(maxsize=1)
def _series_context(cyl, omega, k, n_sol, p_max, count, quad, top):
    a, h, g = cyl.radius, cyl.depth, cyl.gap
    modes = _modes(cyl, omega, n_sol, p_max) if count == n_sol else None
    if modes is not None:
        kn, basis = modes.kn, modes.basis
        big_e, big_f = modes.inner[0][:, 1:], modes.inner[1][:, 1:]
        norms, (i0, i1) = modes.norms[1:], modes.moments
    else:
        kn = _extended_roots(omega, h, wave_numbers(omega, h, n_sol)[1], count)
        basis = gegenbauer_basis(p_max)
        big_e = inner_exterior(cyl, k, kn, basis)[:, 1:]
        big_f = inner_interior(cyl, count, basis)[:, 1:]
        norms = _mode_norms(k, kn, h)[1:]
        i0, i1 = _wall_moments(cyl, k, kn)
    window = _smooth_window(count, NEAR_WINDOW_START)
    models = {f.lam: _corner_models(quad, cyl, count, window, f.lam) for f in basis.families}
    if 0.5 not in models:                                        # the jump of the wall velocity needs it
        models[0.5] = _corner_models(quad, cyl, count, window, 0.5)
    s_wall = quad.wall.z + h
    secant = 1.0 / np.cos(kn * h)
    phase = np.outer(s_wall, kn)
    cos_wall = np.cos(phase) * (secant * window)
    sin_wall = -np.sin(phase) * (kn * secant * window)
    z0, dz0 = _cosh_ratio(k, s_wall, h)
    lam = np.arange(1, count + 1) * math.pi / g
    sign = (-1.0) ** np.arange(1, count + 1)                     # cos(lambda_n g), the bottom z = -d
    rb = quad.bottom.r
    keep = _bottom_keep(quad, cyl, count)
    bessel_a = [ive(q, lam * a) for q in range(top)]             # ive(q, .) for q = 0 ... max order + 1
    bessel_r = [[ive(q, lam[:n] * r) for q in range(top)] for r, n in zip(rb, keep)]
    growth = [np.exp(lam[:n] * (r - a)) for r, n in zip(rb, keep)]
    return _SeriesContext(k, kn, basis, big_e, big_f, norms, (i0, i1), window, models, cos_wall, sin_wall, z0, k * dz0,
                          lam, sign, keep, bessel_a, bessel_r, growth, modes, {})


def _body_series(sol, quad, count=None):
    """Per order m: (phi_m, d phi_m/dz) on the wall nodes, (phi_m, d phi_m/dr) on the bottom nodes, phi_m at the
    waterline, of the extended, corner-corrected series (module docstring, "Near field"). The gap velocity
    sum_p gap[p] u_p and the wall velocity of each order define the exterior and interior coefficients for every n in
    closed form (as in _solve_order); the first N are the solution's own, and the series is carried to n = count. The
    wall and waterline sums are summed with _smooth_window, the bottom sums cut where exp(-lambda_n (a - r)) is below
    exp(-BOTTOM_DECAY); then the exact sum less the windowed (or cut) sum of each model tail is added (_corner_models).
    What does not depend on the solution is cached (_series_context), so the fixed and the free body share it."""
    cyl = sol.cyl
    a, d, h, g = cyl.radius, cyl.draft, cyl.depth, cyl.gap
    count = max(NEAR_SERIES_MODES, len(sol.kn)) if count is None else count
    if count < len(sol.kn):
        raise ValueError("_body_series: count %d is below the solution's N = %d" % (count, len(sol.kn)))
    n_sol = len(sol.kn)
    ctx = _series_context(cyl, float(sol.omega), float(sol.k), n_sol, _spec_key(sol.p_max), count, quad,
                          max(sol.orders) + 2)
    basis, big_e, big_f, norms, models = ctx.basis, ctx.big_e, ctx.big_f, ctx.norms, ctx.models
    i0, i1 = ctx.moments
    cos_wall, sin_wall, z0, dz0, lam, sign = ctx.cos_wall, ctx.sin_wall, ctx.z0, ctx.dz0, ctx.lam, ctx.sign
    rb, keep, bessel_a, bessel_r, growth = quad.bottom.r, ctx.keep, ctx.bessel_a, ctx.bessel_r, ctx.growth
    window = ctx.window
    out = {}
    for m, o in sol.orders.items():
        kappa, mu = ctx.ratios(cyl, m)
        wall = o.wall_velocity[0] * i0[1:] + o.wall_velocity[1] * i1[1:]
        b = (o.gap @ big_e + wall) / (norms * kappa)
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
    """The difference-frequency kernel T-: (fi conj(gj) + conj(fj) gi)/4, elementwise; for i = j, Re(f conj(g))/2,
    the time average of the product of two harmonic quantities f and g."""
    return 0.25 * (fi * np.conj(gj) + np.conj(fj) * gi)


def _pair_sum(fi, gi, fj, gj):
    """The sum-frequency kernel T+: (fi gj + fj gi)/4, elementwise, with the normalisation of _pair: the products
    f g exp(-i (omega_i + omega_j) t), no conjugate, so symmetric in i and j. Used by check_pair_symmetry and qtf_table;
    the values of qtf_table are compared with HAMS-MREL and Nemoh by qtf_compare.py."""
    return 0.25 * (fi * gj + fj * gi)


def _pair_cross(ai, bi, aj, bj, kernel=_pair):
    """(a x b) under the kernel: component l = sum eps_lmn pair(a_m, b_n), 3-vectors."""
    out = np.zeros(3, complex)
    for l, m, n in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        out[l] = kernel(ai[m], bi[n], aj[m], bj[n]) - kernel(ai[n], bi[m], aj[n], bj[m])
    return out


def _load(vector, x):
    """The force and moment densities (6, ...) about xr of a force density (3, ...) at the points x - xr (3, ...); the
    caller integrates them."""
    return np.concatenate([vector, np.cross(x, vector, axis=0)])


def near_field_terms(cyl, omega_i, sol_i, xi_i, omega_j, sol_j, xi_j, quad, xr=(0.0, 0.0, 0.0),
                     kernel=_pair):
    """T-_ij of each near-field term (WAMIT v7.5 user manual, eqs. 15.54 to 15.57, under the kernel _pair; T+_ij with
    kernel=_pair_sum),
    over rho g, per unit wave amplitudes, as (6,) complex arrays (F_x, F_y, F_z, M_x, M_y, M_z, moments about xr), components
    1, 3, 5 filled; for i = j the real part is the mean load over rho g A**2. xi_i, xi_j: the motions (module convention,
    about the origin); xr in the x-z plane. See the module docstring for the terms as coded."""
    xr = np.asarray(xr, dtype=float)
    si = surface_state(cyl, omega_i, sol_i, xi_i, quad, xr)
    same = sol_j is sol_i and omega_j == omega_i and xi_j == xi_i
    sj = si if same else surface_state(cyl, omega_j, sol_j, xi_j, quad, xr)
    return near_field_terms_states(cyl, omega_i, si, omega_j, sj, quad, xr, kernel)


def surface_state(cyl, omega, sol, xi, quad, xr=(0.0, 0.0, 0.0)):
    """The fields of one first-order solution on the quadrature (_surface_state), after the checks of near_field_terms:
    xr in the x-z plane, enough theta points for the orders of sol. The QTF table keeps these, not the solutions."""
    xr = np.asarray(xr, dtype=float)
    if xr[1] != 0.0:
        raise ValueError("near_field_terms: xr must lie in the x-z plane (F_y, M_x, M_z are not computed)")
    if quad.theta.size <= 2 * max(sol.orders) + 4:
        raise ValueError("near_field_terms: %d theta points are too few for order %d" % (quad.theta.size, max(sol.orders)))
    return _surface_state(cyl, omega, sol, xi, quad, xr)


def near_field_terms_states(cyl, omega_i, si, omega_j, sj, quad, xr=(0.0, 0.0, 0.0), kernel=_pair):
    """near_field_terms from the surface states si, sj of surface_state (same quad and xr)."""
    xr = np.asarray(xr, dtype=float)
    grav = make_mesh.G
    dth = quad.theta_weight
    terms = {key: np.zeros(6, complex) for key in NEAR_FIELD_TERMS}

    # I: (1/2) int_WL r' [P(zeta, zeta) - 2 P(zeta, X3)] dl, zeta = (i omega/g) phi, r' = n' or (x - xr) x n', dl = a dtheta
    a = cyl.radius
    ct, st = np.cos(quad.theta), np.sin(quad.theta)
    n_wl = np.array([-ct, -st, 0.0 * ct])                     # unit horizontal normal into the body
    x_wl = np.array([a * ct, a * st, 0.0 * ct]) - xr[:, None]
    zi, zj = 1j * omega_i / grav * si.phi_wl, 1j * omega_j / grav * sj.phi_wl
    bracket = kernel(zi, zi, zj, zj) - 2.0 * kernel(zi, si.disp3_wl, zj, sj.disp3_wl)
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
        q = 0.5 * np.sum(kernel(fi["grad"], fi["grad"], fj["grad"], fj["grad"]), axis=0)
        terms["II"] += -np.sum(_load(n * q, x) * w, axis=(1, 2)) / grav
        # III: -(1/g) int r sum_k P(X_k, d_k phi_t) dS
        q = np.sum(kernel(fi["disp"], gi_t[f], fj["disp"], gj_t[f]), axis=0)
        terms["III"] += -np.sum(_load(n * q, x) * w, axis=(1, 2)) / grav
        for s, (state, omega) in enumerate(((fi, omega_i), (fj, omega_j))):
            dens = _load(n * (-1j * omega * state["phi"]), x)
            total = -np.sum(dens * w, axis=(1, 2)) / grav
            big_g[s] += total[:3]
            big_h[s] += total[3:]
    # IV: force alpha x G; moment alpha x H + xi x G (xi the translation of xr)
    terms["IV"][:3] = _pair_cross(si.alpha, big_g[0], sj.alpha, big_g[1], kernel)
    terms["IV"][3:] = (_pair_cross(si.alpha, big_h[0], sj.alpha, big_h[1], kernel)
                       + _pair_cross(si.xi, big_g[0], sj.xi, big_g[1], kernel))
    terms["V"] = _hydrostatic_term(cyl, si, sj, xr, kernel)
    for key in terms:                                          # components 2, 4, 6 vanish by symmetry
        terms[key][[1, 3, 5]] = 0.0
    return terms


def _hydrostatic_term(cyl, si, sj, xr, kernel=_pair):
    """V: F_S and M_S of the WAMIT v7.5 user manual, eqs. (15.56) and (15.57), the formula HAMS-MREL uses, over
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
        return kernel(u[0], v[0], u[1], v[1])

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
    linear = bool(ratios) and all(LINEARITY_RATIO[0] <= q <= LINEARITY_RATIO[1] for q in ratios)
    ok = worst_mod < 1e-3 and worst_complex < 1e-3 and linear
    return ok, ("extrapolated to g = 0: max ||X1| - |MF||/|MF| %.1e, max |X1 - MF|/|MF| %.1e; at g = 1e-3 h itself "
                "max ||X1| - |MF||/|MF| %.1e; |X(2g) - MF|/|X(g) - MF| %s%s"
                % (worst_mod, worst_complex, worst_raw,
                   "%.3f ... %.3f" % (min(ratios), max(ratios)) if ratios else "not tested anywhere",
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


# The three-family case of the radiation and energy checks: one frequency, gap basis (20, 4, 4), N and 2 N
THREE_FAMILY_CASE = dict(omega=0.70, n_modes=8000, p_max=(20, 4, 4))


def _first_order_extrapolated(omega, n=None, p_max=None):
    """_first_order at N and 2 N (default: CHECK_RESOLUTION's), extrapolated in N: (2**p Q(2N) - Q(N))/(2**p - 1),
    p = RICHARDSON_EXPONENT. Also returns the values at 2 N."""
    n = CHECK_RESOLUTION["n_modes"] if n is None else n
    p_max = CHECK_RESOLUTION["p_max"] if p_max is None else p_max
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
    normalisation and the sign and conjugation of excitation and added_mass_damping, now after extrapolation in N.
    Pitch and heave have the interior particular solution, which is not in the span of the N truncated interior modes;
    the residual of the identities is the tail of its cosine projection (coefficients n**-2) against that of the corner-singular gap velocity
    (n**-2/3), so it falls as N**(-5/3): measured exponent 1.68 to 1.73 for N = 1000 ... 8000, 1.3e-7 for A15 = A51 still at
    N = 32000. One Richardson step with p = 5/3 removes that term; what remains falls about as N**-2.9 (2e-8 at N = 4000,
    8000). So A15 = A51, B15 = B51 and Haskind for modes 3 and 5 test the particular solutions through their gap
    projections (potential_gap, flux_gap) and the pitch wall term. They do not test _Particular.bottom: the bottom
    integral of phi_p is real and enters only the diagonal Re A33 and Re A55, which no identity here contains;
    check_particular_bottom tests it against quad. They do not bound the error of A, B or X themselves, which converge as
    N**(-5/3) with MATRIX_TAIL (about N**(-1.6) without it; Convergence)."""
    worst_sym = worst_haskind = raw_sym = raw_haskind = 0.0
    for omega in drift_convergence.WIDE_FREQUENCIES:
        q, hi = _first_order_extrapolated(omega)
        sym, haskind = _reciprocity_errors(q)
        worst_sym, worst_haskind = max(worst_sym, sym), max(worst_haskind, haskind)
        sym, haskind = _reciprocity_errors(hi)
        raw_sym, raw_haskind = max(raw_sym, sym), max(raw_haskind, haskind)
    case = THREE_FAMILY_CASE
    q, hi = _first_order_extrapolated(case["omega"], case["n_modes"], case["p_max"])
    three_sym, three_haskind = _reciprocity_errors(q)
    ok = max(worst_sym, three_sym) < 1e-8 and max(worst_haskind, three_haskind) < 1e-6
    return ok, ("extrapolated in N: max |A15 - A51|/|A15|, |B15 - B51|/|B15| %.1e, max |X Haskind - X diffraction|/|X| "
                "%.1e, modes 1, 3, 5; at N = %d itself %.1e and %.1e; gap basis %s at %.2f rad/s, N = %d and %d, "
                "extrapolated: %.1e and %.1e (at N = %d itself %.1e and %.1e)"
                % (worst_sym, worst_haskind, 2 * CHECK_RESOLUTION["n_modes"], raw_sym, raw_haskind, case["p_max"],
                   case["omega"], case["n_modes"], 2 * case["n_modes"], three_sym, three_haskind,
                   2 * case["n_modes"], *_reciprocity_errors(hi)))


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
    case = THREE_FAMILY_CASE
    q, hi = _first_order_extrapolated(case["omega"], case["n_modes"], case["p_max"])
    three, three_raw = error(q), error(hi)
    return max(worst, three) < 1e-6, ("extrapolated in N: max |B_jj - flux|/|flux| %.1e, j = 1, 3, 5; at N = %d itself "
                                      "%.1e; gap basis %s at %.2f rad/s, N = %d and %d, extrapolated: %.1e (at N = %d "
                                      "itself %.1e)" % (worst, 2 * CHECK_RESOLUTION["n_modes"], raw, case["p_max"],
                                                        case["omega"], case["n_modes"], 2 * case["n_modes"], three,
                                                        2 * case["n_modes"], three_raw))


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
    cases = [(omega, res["n_modes"], res["p_max"]) for omega in drift_convergence.WIDE_FREQUENCIES]
    three = THREE_FAMILY_CASE
    cases.append((three["omega"], 2 * three["n_modes"], three["p_max"]))
    errors = []
    for omega, n_modes, p_max in cases:
        diff = solve_diffraction(cyl, omega, res["orders"], n_modes, p_max)
        rads = {j: solve_radiation(cyl, omega, j, n_modes, p_max) for j in MOTION_MODES}
        total = total_solution(diff, rads, motions(cyl, omega, diff, rads))
        a = far_field_amplitude(total)
        full = np.concatenate([a[:0:-1], a])
        flux = 2.0 * math.pi * np.sum(np.abs(full) ** 2)
        optical = -math.sqrt(8.0 * math.pi / total.k) * np.real(np.exp(0.75j * math.pi) * np.sum(full))
        fixed = 2.0 * math.pi * np.sum(np.abs(far_field_amplitude(diff)) ** 2 * np.where(np.arange(len(a)), 2.0, 1.0))
        errors.append(abs(flux - optical) / fixed)
    return max(errors) < 1e-10, ("max |flux - optical|/(flux of the fixed body) %.1e, free body; gap basis %s at %.2f "
                                 "rad/s, N = %d: %.1e" % (max(errors[:-1]), three["p_max"], three["omega"],
                                                          2 * three["n_modes"], errors[-1]))


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


def _orders_key(orders):
    """orders as a hashable cache key: an int M, or a tuple of orders."""
    return int(orders) if isinstance(orders, (int, np.integer)) else tuple(int(m) for m in orders)


@functools.lru_cache(maxsize=1)
def _solved(cyl, omega, n_modes, p_max, orders):
    """(diff, rads): the diffraction solution and the radiation solutions of modes 1, 3, 5 at one resolution, cached (the
    last one), so that the fixed and the free body at one omega share them. Keys from _spec_key and _orders_key."""
    diff = solve_diffraction(cyl, omega, orders, n_modes, p_max)
    return diff, {j: solve_radiation(cyl, omega, j, n_modes, p_max) for j in MOTION_MODES}


def _first_order_case(cyl, omega, motion, n_modes, p_max, orders):
    """(total, xi) of the default problem at one resolution: diffraction, and for "free" the radiation solutions, the
    motions and the total potential (_solved)."""
    diff, rads = _solved(cyl, float(omega), int(n_modes), _spec_key(p_max), _orders_key(orders))
    if motion == "fixed":
        return diff, motions(cyl, omega, diff, None, "fixed")
    xi = motions(cyl, omega, diff, rads)
    return total_solution(diff, rads, xi), xi


# The resolution of the reference values (reference): N exterior and interior modes, the gap basis (P_1/6, P_1/2, P_5/6)
# (Gap basis), azimuthal orders 0 ... M. The coarsest that meets the gate of fx_identity (the free body needs the three
# families and N = 64000) and the 1e-4 spread of convergence (module docstring, "Reference values").
DEFAULT_RESOLUTION = dict(n_modes=64000, p_max=(20, 4, 4), orders=15)
# Resolution of fx_identity: the default, so that its rows are the reference values (cached in _reference_row)
NEAR_CHECK_RESOLUTION = DEFAULT_RESOLUTION
# The (N, gap basis) ladder of convergence, increasing, with DEFAULT_RESOLUTION in it; the orders are tested by doubling M
CONVERGENCE_LADDER = (dict(n_modes=32000, p_max=(10, 2, 2), orders=15), DEFAULT_RESOLUTION,
                      dict(n_modes=128000, p_max=(30, 6, 6), orders=15))
# convergence: below CONVERGENCE_FLOOR times the largest abs(value) of the quantity and motion over CONVERGENCE_FREQUENCIES,
# a value counts as near zero and its spreads are measured against that floor (the floor of drift_convergence.percent);
# it acts nowhere in the present table
CONVERGENCE_FLOOR = 1e-4
# The stated relative uncertainty of the reference values: DEFAULT_UNCERTAINTY (the spread gate of convergence) except
# where KNOWN_UNCERTAINTY names a larger one, (motion, quantity, omega) -> relative uncertainty. The free body at 0.60 to
# 0.80 rad/s carries the systematic residual of the F_x identity (module docstring, "Reference values"): near-field F_x
# there is +-5e-4, and F_z at 0.70 and 0.80 rad/s is estimated at 2.5e-3 and 3.3e-3 (check_uncertainty reproduces the
# estimates and fails if one exceeds its entry). Free near-field F_x at 1.40, 1.60 and 1.80 rad/s is 2e-4: its residual
# against far F_x (which changes by less than 1.2e-7 over the resolutions) is -1.11e-4, -1.29e-4, -1.26e-4 (fx_identity);
# free F_z and M_y there have no separate estimate (check_uncertainty covers 0.60 to 0.80 rad/s only). convergence passes
# a row named here against this value instead of the gate, and prints it as a known exceedance. Every other row of the
# default hull at DEFAULT_RESOLUTION and WIDE_FREQUENCIES carries DEFAULT_UNCERTAINTY; any other row nan (not measured).
DEFAULT_UNCERTAINTY = 1e-4
KNOWN_UNCERTAINTY = {("free", "Fx_near", 0.60): 5e-4, ("free", "Fx_near", 0.70): 5e-4, ("free", "Fx_near", 0.80): 5e-4,
                     ("free", "Fx_near", 1.40): 2e-4, ("free", "Fx_near", 1.60): 2e-4, ("free", "Fx_near", 1.80): 2e-4,
                     ("free", "Fz", 0.70): 3e-3, ("free", "Fz", 0.80): 4e-3}
UNCERTAIN_QUANTITIES = ("Fx_far", "Fx_near", "Fz", "My")   # the quantities with an uncertainty column in the CSV
UNCERTAINTY_FREQUENCY_MATCH = 1e-9   # uncertainty: omega is a WIDE_FREQUENCIES entry if within this relative distance


def _is_default_resolution(n_modes, p_max, orders):
    res = DEFAULT_RESOLUTION
    return ((int(n_modes), _spec_key(p_max), _orders_key(orders))
            == (res["n_modes"], _spec_key(res["p_max"]), _orders_key(res["orders"])))


def uncertainty(motion, quantity, omega, cyl, n_modes, p_max, orders):
    """The stated relative uncertainty of a reference value: KNOWN_UNCERTAINTY, else DEFAULT_UNCERTAINTY, for the
    default hull (cyl == Cylinder(), h/a = 10) at DEFAULT_RESOLUTION and an omega of WIDE_FREQUENCIES (to a relative
    UNCERTAINTY_FREQUENCY_MATCH), the only rows the evidence of "Reference values" measures; else nan (not measured).
    The N needed grows with h/a (measured, free F_z at 0.90 rad/s, gap basis (20, 4, 4), between N = 4000 and 16000:
    a 10, d 10, h 500 changes by 4.9 %; a 5, d 8, h 40 by 0.8 %), so a new geometry needs convergence and fx_identity
    re-run before any uncertainty is stated for it."""
    if cyl != Cylinder() or not _is_default_resolution(n_modes, p_max, orders):
        return float("nan")
    for w in drift_convergence.WIDE_FREQUENCIES:
        if abs(float(omega) - w) <= UNCERTAINTY_FREQUENCY_MATCH * w:
            return KNOWN_UNCERTAINTY.get((motion, quantity, w), DEFAULT_UNCERTAINTY)
    return float("nan")


FREE_BODY_UNCERTAINTY = (
    "Stated uncertainty (KNOWN_UNCERTAINTY; module docstring, \"Reference values\"): relative %g everywhere except the "
    "free body at 0.60 to 0.80 rad/s, where near-field F_x settles about 4.2e-4 of far F_x from it (pitch only, "
    "systematic), and free near-field F_x at 1.40 to 1.80 rad/s, 1.1e-4 to 1.3e-4 from far F_x (fx_identity): %s. "
    "Free F_z and M_y at 1.40 to 1.80 rad/s have no separate estimate. The estimates come from check_uncertainty (F_z and M_y against the F_x residual across gap bases); "
    "the CSV carries them per row (columns *_unc)."
    % (DEFAULT_UNCERTAINTY, ", ".join("%s %s %.2f: %g" % (k[0], k[1], k[2], v) for k, v in sorted(KNOWN_UNCERTAINTY.items()))))

REFERENCE_PAIRS = (11, 13, 15, 33, 35, 55)


@functools.lru_cache(maxsize=None)
def _reference_row(cyl, omega, motion, n_modes, p_max, orders):
    """One row of reference (keys: _spec_key, _orders_key), cached: fx_identity, convergence and the reference
    subcommand share the rows at DEFAULT_RESOLUTION."""
    diff, rads = _solved(cyl, omega, n_modes, p_max, orders)
    a, b = added_mass_damping(cyl, omega, rads)
    x = excitation(cyl, diff)
    xi = motions(cyl, omega, diff, rads, motion)
    total = total_solution(diff, rads, xi)
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    terms, near = near_field(cyl, omega, total, xi, quad)
    row = dict(omega=omega, motion=motion, n_modes=n_modes, p_max=p_max, orders=orders)
    for ij in REFERENCE_PAIRS:
        i, j = divmod(ij, 10)
        row["A%d" % ij], row["B%d" % ij] = float(a[i - 1, j - 1]), float(b[i - 1, j - 1])
    for j in MOTION_MODES:
        row["X%d" % j] = abs(x[j])
    for j in MOTION_MODES:
        row["RAO%d" % j] = abs(xi[j]) if motion == "free" else float("nan")
    row["Fx_far"] = far_field_fx(cyl, omega, total)
    row["Fx_near"], row["Fz"], row["My"] = float(near[0]), float(near[2]), float(near[4])
    for name, index in (("Fx", 0), ("Fz", 2), ("My", 4)):
        for key in NEAR_FIELD_TERMS:
            row["%s_%s" % (name, key)] = float(terms[key][index])
    for quantity in UNCERTAIN_QUANTITIES:
        row["%s_unc" % quantity] = uncertainty(motion, quantity, omega, cyl, n_modes, p_max, orders)
    return row


def reference(cyl, frequencies, motion, n_modes, p_max, orders):
    """One row per omega of frequencies, motion "fixed" or "free", at the resolution (n_modes, p_max, orders): omega,
    motion and the resolution; A and B (AmssDamp.1's normalisation, A/rho and B/(rho omega), keys A11, A13, A15, A33,
    A35, A55 and B..); abs(X) (ExcForce.3, over rho g A; X1, X3, X5); abs(xi) (Motion.4's MOD, rotations in rad/m;
    RAO1, RAO3, RAO5, nan for the fixed body); Fx_far (MeanDriftMomentum.8), Fx_near, Fz, My (MeanDriftPressure.9, over
    rho g A**2, moments in m about the origin) and each near-field term of them, Fx_I ... Fx_V, Fz_.., My_..
    (MeanDriftTerms.txt); and Fx_far_unc, Fx_near_unc, Fz_unc, My_unc, the stated relative uncertainty of each
    (uncertainty: KNOWN_UNCERTAINTY or DEFAULT_UNCERTAINTY for the default hull at DEFAULT_RESOLUTION and
    WIDE_FREQUENCIES, else nan, not measured). p_max is the key of _spec_key, a tuple for the three
    families; the CSV writes it space-separated ("20 4 4"). See the module docstring, "Reference values"."""
    if motion not in ("fixed", "free"):
        raise ValueError("reference: motion is 'fixed' or 'free', got %r" % (motion,))
    return [dict(_reference_row(cyl, float(omega), motion, int(n_modes), _spec_key(p_max), _orders_key(orders)))
            for omega in frequencies]


FX_FLOOR = 1e-4   # fx_identity: below FX_FLOOR times the largest abs(far F_x) of the motion, the error is over that floor


def check_fx_identity():
    """Near-field F_x (sum of the terms) against far-field F_x, fixed and free, every WIDE_FREQUENCIES, at
    NEAR_CHECK_RESOLUTION, to 5e-4 of max(abs(far F_x), FX_FLOOR scale), scale the largest abs(far F_x) of the motion
    (the floor of drift_convergence.percent; it acts only at free 0.40 rad/s, far F_x = 3e-6); see the module docstring."""
    cyl = Cylinder()
    res = NEAR_CHECK_RESOLUTION
    rows = {}
    for omega in drift_convergence.WIDE_FREQUENCIES:          # omega outermost: fixed and free share the solves
        for motion in ("fixed", "free"):
            rows[motion, omega] = reference(cyl, [omega], motion, res["n_modes"], res["p_max"], res["orders"])[0]
    worst, text = 0.0, []
    for motion in ("fixed", "free"):
        scale = max(abs(rows[motion, w]["Fx_far"]) for w in drift_convergence.WIDE_FREQUENCIES)
        for omega in drift_convergence.WIDE_FREQUENCIES:
            near, far = rows[motion, omega]["Fx_near"], rows[motion, omega]["Fx_far"]
            floored = abs(far) < FX_FLOOR * scale
            err = abs(near - far) / max(abs(far), FX_FLOOR * scale)
            worst = max(worst, err)
            text.append("%s %.2f %.1e%s" % (motion, omega, err, " (floor)" if floored else ""))
    return worst < 5e-4, "max |near - far|/|far| %.1e; per case: %s" % (worst, ", ".join(text))


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
    """T-_ij = conj(T-_ji) and T+_ij = T+_ji (kernel _pair_sum) for omega_i = 0.70, omega_j = 0.95, free body, every
    term and component, to 1e-12 of the largest abs(T_ij) of the term; see the module docstring."""
    cyl = Cylinder()
    res = LIGHT_RESOLUTION
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    (si, xi_i), (sj, xi_j) = (_first_order_case(cyl, w, "free", res["n_modes"], res["p_max"], res["orders"])
                              for w in (0.70, 0.95))
    tij = near_field_terms(cyl, 0.70, si, xi_i, 0.95, sj, xi_j, quad)
    tji = near_field_terms(cyl, 0.95, sj, xi_j, 0.70, si, xi_i, quad)
    sij = near_field_terms(cyl, 0.70, si, xi_i, 0.95, sj, xi_j, quad, kernel=_pair_sum)
    sji = near_field_terms(cyl, 0.95, sj, xi_j, 0.70, si, xi_i, quad, kernel=_pair_sum)
    worst = worst_sum = 0.0
    for key in NEAR_FIELD_TERMS:
        scale, scale_sum = np.max(np.abs(tij[key])), np.max(np.abs(sij[key]))
        if scale > 0.0:
            worst = max(worst, np.max(np.abs(tij[key] - np.conj(tji[key]))) / scale)
        if scale_sum > 0.0:
            worst_sum = max(worst_sum, np.max(np.abs(sij[key] - sji[key])) / scale_sum)
    return max(worst, worst_sum) < 1e-12, ("max |T-_ij - conj(T-_ji)|/max |T-_ij| over the terms %.1e; "
                                           "max |T+_ij - T+_ji|/max |T+_ij| %.1e" % (worst, worst_sum))


QTF_FREQUENCIES = tuple(round(0.40 + 0.05 * i, 2) for i in range(21))   # the grid of the Nemoh runs (NEMOH_FREQUENCY_RANGE)
QTF_COLUMNS = ("motion", "omega_i", "omega_j", "kernel", "Fx_re", "Fx_im", "Fz_re", "Fz_im", "My_re", "My_im")


def qtf_table(cyl, frequencies, motion, n_modes, p_max, orders, xr=(0.0, 0.0, 0.0)):
    """Quadratic QTF of the near-field terms for every pair omega_i <= omega_j of frequencies, both kernels: rows of
    QTF_COLUMNS with kernel "diff" (T-_ij, _pair) or "sum" (T+_ij, _pair_sum), F/(rho g), M_y in m about xr, per unit
    amplitudes, in the internal e^{-i omega t} convention (the conjugates of QTF_Quadratic.10d/.10s). The other half of
    the square follows: T-_ji = conj(T-_ij), T+_ji = T+_ij. One first-order solve per frequency; only the surface states
    are kept.

    Phase of the wave. The solutions here take HAMS's VINP as the incident potential, so their wave has elevation
    eta_I = -i A e^(i k x) (module docstring, "Incident wave"), while the QTFs of HAMS-MREL and Nemoh are per unit
    amplitude of a wave with eta_I = A e^(i k x) at the origin (spec 5.2: the physical potential is i VINP). Every
    first-order quantity here is therefore -i times that of the unit wave. The factor cancels in T- (one factor
    conjugated) and in the mean drift, but gives (-i)**2 = -1 in T+, so the sum rows are multiplied by -1."""
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    states = []
    for omega in frequencies:
        total, xi = _first_order_case(cyl, omega, motion, n_modes, p_max, orders)
        states.append(surface_state(cyl, omega, total, xi, quad, xr))
    rows = []
    for i, wi in enumerate(frequencies):
        for j in range(i, len(frequencies)):
            wj = frequencies[j]
            for name, kernel in (("diff", _pair), ("sum", _pair_sum)):
                terms = near_field_terms_states(cyl, wi, states[i], wj, states[j], quad, xr, kernel)
                total = sum(terms.values()) * (-1.0 if name == "sum" else 1.0)   # the phase of the wave, above
                rows.append(dict(motion=motion, omega_i=wi, omega_j=wj, kernel=name,
                                 Fx_re=total[0].real, Fx_im=total[0].imag, Fz_re=total[2].real, Fz_im=total[2].imag,
                                 My_re=total[4].real, My_im=total[4].imag))
    return rows


def check_qtf_diagonal():
    """The diagonal of qtf_table (kernel diff) is the mean drift of near_field, fixed and free at 0.70 and 0.95 rad/s,
    LIGHT_RESOLUTION, to 1e-12. Both go through near_field_terms_states, so this checks the table's bookkeeping only; the
    sign of its sum rows is checked against HAMS-MREL and Nemoh by qtf_compare.py (docs/second_order_validation.md)."""
    cyl, res = Cylinder(), LIGHT_RESOLUTION
    quad = body_quadrature(cyl, NEAR_QUADRATURE["n_corner"], NEAR_QUADRATURE["n_far"])
    worst = 0.0
    for motion in ("fixed", "free"):
        rows = qtf_table(cyl, (0.70, 0.95), motion, res["n_modes"], res["p_max"], res["orders"])
        for omega in (0.70, 0.95):
            total, xi = _first_order_case(cyl, omega, motion, res["n_modes"], res["p_max"], res["orders"])
            mean = near_field(cyl, omega, total, xi, quad)[1]
            row = [r for r in rows if r["kernel"] == "diff" and r["omega_i"] == omega and r["omega_j"] == omega][0]
            got = np.array([row["Fx_re"], row["Fz_re"], row["My_re"]])
            worst = max(worst, np.max(np.abs(got - mean[[0, 2, 4]])) / np.max(np.abs(mean)),
                        abs(row["Fx_im"]) / np.max(np.abs(mean)))
    return worst < 1e-12, "max |diagonal - mean drift|/max |mean| %.1e" % worst


def run_qtf(out, motions=("fixed", "free"), frequencies=QTF_FREQUENCIES, resolution=None):
    """The qtf subcommand: qtf_table at DEFAULT_RESOLUTION on QTF_FREQUENCIES for each motion, written to out."""
    res = resolution or DEFAULT_RESOLUTION
    cyl, start = Cylinder(), time.time()
    rows = []
    for motion in motions:
        rows += qtf_table(cyl, frequencies, motion, res["n_modes"], res["p_max"], res["orders"])
        print("qtf %s: %d rows, %.0f s" % (motion, len(rows), time.time() - start), flush=True)
    with open(out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(QTF_COLUMNS)
        for row in rows:
            writer.writerow([row[c] if isinstance(row[c], str) else "%.10g" % row[c] for c in QTF_COLUMNS])
    print("wrote %d rows to %s" % (len(rows), out))
    return 0


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


POLYLOG_ORDERS = (1.0, 2.0 / 3.0, 7.0 / 6.0, 4.0 / 3.0, 5.0 / 3.0, 11.0 / 6.0, 13.0 / 6.0)
POLYLOG_POINTS = (-0.05, -0.5, -1.0, -1.5, -3.0, -12.0, -60.0, -1.0 + 3.0j, -0.05 + 0.1j, -0.05 - 1.0j, -0.05 + 3.1j,
                  -2.0 + 3.1j)


def check_polylog():
    """_polylog(s, mu) against the direct sum sum_(n=1)^(4000) n**(-s) e^(n mu), to 1e-10 relative, for the orders of
    _corner_models (1 and lam + 1/2, lam + 3/2 of the three families) at real mu from -0.05 to -60 (the bottom tails,
    mu = -pi rho/g, reach -pi a/(2 g), large for a small gap) and complex mu with Re mu = -0.05 ... -2 and Im mu up to
    3.1 (Re mu <= -0.05, so that 4000 terms reach round-off); and that it raises outside its domain (Re mu > 0). It tests
    both branches of _polylog (zeta series and direct sum) and the switch between them, not the unit circle Re mu = 0,
    where the direct sum does not converge absolutely; the wall tails there are tested through the wall series
    (fx_identity) and corner_limit."""
    n = np.arange(1, 4001)
    worst = 0.0
    for order in POLYLOG_ORDERS:
        for mu in POLYLOG_POINTS:
            direct = np.sum(n ** (-order) * np.exp(n * mu))
            worst = max(worst, abs(complex(_polylog(order, mu)) - direct) / abs(direct))
    try:
        _polylog(2.0 / 3.0, 0.5)
        raised = False
    except ValueError:
        raised = True
    return worst < 1e-10 and raised, ("max |Li_s(e^mu) - direct sum|/|direct sum| %.1e over %d orders and %d points; "
                                      "raises for Re mu > 0: %s" % (worst, len(POLYLOG_ORDERS), len(POLYLOG_POINTS),
                                                                   raised))


MATCHING_TAIL_CASE = dict(omega=0.70, n_modes=2000, p_max=(20, 4, 4), brute=64000)
MATCHING_TAIL_TOLERANCE = 5e-3   # measured 1.1e-3 (matrix) and 1.8e-3 (wall)


def check_matching_tail():
    """_matching_tail against brute-force sums, and where it applies (_tail_rows); see the module docstring."""
    cyl = Cylinder()
    case = MATCHING_TAIL_CASE
    n, big = case["n_modes"], case["brute"]
    basis = gegenbauer_basis(case["p_max"])
    k, kn = wave_numbers(case["omega"], cyl.depth, big)
    big_e, big_f = _inner_products(cyl, k, kn, basis)
    norms = _mode_norms(k, kn, cyl.depth)
    i0, _ = _wall_moments(cyl, k, kn)
    tail_n, wall_n = _matching_tail(cyl, basis, n)
    tail_big, wall_big = _matching_tail(cyl, basis, big)
    model, model_wall = tail_n - tail_big, wall_n - wall_big          # the model's sum over N < n <= brute
    sl = slice(n + 1, big + 1)
    worst = worst_wall = 0.0
    for m in (0, 1):
        kappa, mu = _radial_ratios(cyl, m, k, kn)
        ext = big_e[:, sl] / (norms[sl] * kappa[sl])
        brute = ext @ big_e[:, sl].T - (big_f[:, sl] / (0.5 * cyl.gap * mu[sl])) @ big_f[:, sl].T
        brute_wall = ext @ i0[sl]                                     # a unit wall velocity: the moment i0
        worst = max(worst, np.max(np.abs(model - brute)) / np.max(np.abs(brute)))
        worst_wall = max(worst_wall, np.max(np.abs(model_wall - brute_wall)) / np.max(np.abs(brute_wall)))
    small = Cylinder(draft=cyl.depth - 0.1)
    rows_small = _tail_rows(small, gegenbauer_basis(10), 1000)
    tail_small, wall_small = _matching_tail(small, gegenbauer_basis(10), 1000)
    restricted = not rows_small.any() and not np.any(tail_small) and not np.any(wall_small)
    matching = dict(n_modes=MATCHING_MODES_PER_P * max(MATCHING_P), p_max=max(MATCHING_P))
    everywhere = all(_tail_rows(cyl, gegenbauer_basis(res["p_max"]), res["n_modes"]).all()
                     for res in (LIGHT_RESOLUTION, CHECK_RESOLUTION, matching, DEFAULT_RESOLUTION))
    ok = worst < MATCHING_TAIL_TOLERANCE and worst_wall < MATCHING_TAIL_TOLERANCE and restricted and everywhere
    return ok, ("gap basis %s, N = %d, sums to n = %d, orders 0 and 1: max |model - brute|/max |brute| %.1e (matrix), "
                "%.1e (wall); gap 0.1 m, N = 1000, P = 10: no row kept %s; default gap at LIGHT, CHECK, matching and "
                "default resolutions: every row kept %s" % (case["p_max"], n, big, worst, worst_wall, restricted,
                                                            everywhere))


CORNER_LIMIT_RHO = (1e-3, 1e-4, 1e-5)    # distances from the corner of the probe nodes, m
CORNER_LIMIT_RESOLUTION = dict(n_modes=16000, p_max=(20, 4, 4), orders=3)
CORNER_LIMIT_TOLERANCE = 1e-3            # measured 1.9e-4 (fixed and free, 0.70 and 1.40 rad/s, orders 0 ... 3)


def _probe_quadrature(cyl, rho):
    """A Quadrature whose wall and bottom nodes lie at the distances rho from the corner, with zero weights: for
    evaluating the body series there, not for integrating."""
    a, d = cyl.radius, cyl.draft
    rho = np.asarray(rho, dtype=float)
    wall = Face("wall", np.full(rho.shape, a), -d + rho, np.zeros_like(rho), (-1.0, 0.0), rho)
    bottom = Face("bottom", a - rho, np.full(rho.shape, -d), np.zeros_like(rho), (0.0, 1.0), rho)
    return Quadrature(wall, bottom, np.zeros(1), 2.0 * math.pi)


def check_corner_limit():
    """rho**(1/3) times the tangential velocity on the wall (d phi_m/dz, exterior series) and on the bottom (d phi_m/dr,
    interior series), each extrapolated to the corner, agree; see the module docstring."""
    cyl = Cylinder()
    rho = np.array(CORNER_LIMIT_RHO)
    s = rho ** (1.0 / 3.0)
    fit = np.vander(s, 3, increasing=True)                # c0 + c1 s + c2 s**2, s = rho**(1/3)
    quad = _probe_quadrature(cyl, rho)
    res = CORNER_LIMIT_RESOLUTION
    worst, cases = 0.0, 0
    for omega in (0.70, 1.40):
        for motion in ("fixed", "free"):
            total, _ = _first_order_case(cyl, omega, motion, res["n_modes"], res["p_max"], res["orders"])
            for m, values in _body_series(total, quad).items():
                wall = np.linalg.lstsq(fit, values["wall"][1] * s, rcond=None)[0][0]
                bottom = np.linalg.lstsq(fit, values["bottom"][1] * s, rcond=None)[0][0]
                worst = max(worst, abs(wall - bottom) / abs(wall))
                cases += 1
    return worst < CORNER_LIMIT_TOLERANCE, ("max |c_wall - c_bottom|/|c_wall| %.1e over %d cases (corner limits of "
                                           "rho**(1/3) u_t, fitted at rho = %s m)"
                                           % (worst, cases, ", ".join("%g" % r for r in CORNER_LIMIT_RHO)))


CONVERGENCE_FREQUENCIES = (0.70, 0.95, 1.40)
CONVERGENCE_QUANTITIES = ("Fx_far", "Fz", "My", "B11")
CONVERGENCE_SPREAD = 1e-4


def convergence_orders():
    """DEFAULT_RESOLUTION with twice its orders."""
    return dict(DEFAULT_RESOLUTION, orders=2 * DEFAULT_RESOLUTION["orders"])


def convergence_table(cyl=None):
    """For CONVERGENCE_QUANTITIES at CONVERGENCE_FREQUENCIES, fixed and free: the values at the resolutions of
    CONVERGENCE_LADDER and at DEFAULT_RESOLUTION with twice its orders, and the spreads of the default against the
    finest of the ladder and against 2M, each over max(abs(default value), floor), floor = CONVERGENCE_FLOOR times the
    largest abs(default value) of the quantity and motion over the frequencies (near-zero values only); "allowed" is the
    larger of CONVERGENCE_SPREAD and the default row's stated uncertainty (its *_unc, uncertainty). A list of dicts,
    one per (motion, quantity, omega)."""
    cyl = Cylinder() if cyl is None else cyl
    resolutions = list(CONVERGENCE_LADDER) + [convergence_orders()]
    values = {}
    for omega in CONVERGENCE_FREQUENCIES:                 # omega outermost: the solves are cached per omega
        for i, res in enumerate(resolutions):
            for motion in ("fixed", "free"):
                values[omega, motion, i] = _reference_row(cyl, omega, motion, res["n_modes"], _spec_key(res["p_max"]),
                                                          _orders_key(res["orders"]))
    default, finest = CONVERGENCE_LADDER.index(DEFAULT_RESOLUTION), len(CONVERGENCE_LADDER) - 1
    out = []
    for motion in ("fixed", "free"):
        for quantity in CONVERGENCE_QUANTITIES:
            scale = max(abs(values[w, motion, default][quantity]) for w in CONVERGENCE_FREQUENCIES)
            for omega in CONVERGENCE_FREQUENCIES:
                v = [values[omega, motion, i][quantity] for i in range(len(resolutions))]
                ref, floor = v[default], CONVERGENCE_FLOOR * scale
                spread_n = abs(v[finest] - ref) / max(abs(ref), floor)
                spread_m = abs(v[-1] - ref) / max(abs(ref), floor)
                stated = values[omega, motion, default].get("%s_unc" % quantity, float("nan"))
                allowed = CONVERGENCE_SPREAD if math.isnan(stated) else max(CONVERGENCE_SPREAD, stated)
                out.append(dict(omega=omega, motion=motion, quantity=quantity, values=v, spread_n=spread_n,
                                spread_m=spread_m, floored=abs(ref) < floor, allowed=allowed,
                                known=max(spread_n, spread_m) >= CONVERGENCE_SPREAD and allowed > CONVERGENCE_SPREAD,
                                finite=bool(np.all(np.isfinite(v)))))
    return out


def format_convergence_table(table):
    """The convergence table as text."""
    names = ["N %d, gap basis %s, M %d" % (r["n_modes"], r["p_max"], r["orders"]) for r in CONVERGENCE_LADDER]
    names.append("the default with 2M")
    default, finest = CONVERGENCE_LADDER.index(DEFAULT_RESOLUTION), len(CONVERGENCE_LADDER) - 1
    lines = ["Convergence (F over rho g A**2, moments in m; B11 over rho omega). Resolutions: "
             + "; ".join("[%d] %s" % (i, n) for i, n in enumerate(names)) + " ([%d] is DEFAULT_RESOLUTION)." % default,
             "Spreads: |[%d] - [%d]| (N, gap basis) and |[%d] - [%d]| (orders), over max(|[%d]|, floor), floor = %g of "
             "the largest |[%d]| of the quantity and motion; allowed: %g, or the stated uncertainty where larger."
             % (finest, default, len(names) - 1, default, default, CONVERGENCE_FLOOR, default, CONVERGENCE_SPREAD),
             "%-6s %5s %-6s %s %10s %10s %10s" % ("motion", "omega", "qty", " ".join("%17s" % ("[%d]" % i)
                                                                                 for i in range(len(names))),
                                                "spread NP", "spread M", "allowed")]
    for row in table:
        note = " known exceedance" if row["known"] else ""
        note += " (floor)" if row["floored"] else ""
        lines.append("%-6s %5.2f %-6s %s %10.1e %10.1e %10.1e%s" % (
            row["motion"], row["omega"], row["quantity"], " ".join("%17.10g" % v for v in row["values"]),
            row["spread_n"], row["spread_m"], row["allowed"], note))
    lines.append(FREE_BODY_UNCERTAINTY)
    return "\n".join(lines)


def check_convergence():
    """Every value of the convergence table is finite, and both spreads of every row are below CONVERGENCE_SPREAD, or
    below the row's stated uncertainty where KNOWN_UNCERTAINTY names a larger one (a known exceedance: it still fails if
    the spread grows past it). The message carries the table. See the module docstring, "Reference values"."""
    table = convergence_table()
    finite = all(row["finite"] for row in table)
    ok = finite and all(max(row["spread_n"], row["spread_m"]) < row["allowed"] for row in table)
    plain = [row for row in table if not row["known"]]
    known = ["%s %.2f %s %.1e (stated %.0e)" % (row["motion"], row["omega"], row["quantity"],
                                               max(row["spread_n"], row["spread_m"]), row["allowed"])
             for row in table if row["known"]]
    floored = ["%s %.2f %s" % (row["motion"], row["omega"], row["quantity"]) for row in table if row["floored"]]
    return ok, ("all finite %s; max spread %.1e (gate %.0e) outside the known exceedances; known exceedances: %s; "
                "near-zero floor applied at: %s\n%s"
                % (finite, max(max(r["spread_n"], r["spread_m"]) for r in plain), CONVERGENCE_SPREAD,
                   ", ".join(known) or "none", ", ".join(floored) or "none", format_convergence_table(table)))


UNCERTAINTY_FREQUENCIES = (0.60, 0.70, 0.80)   # check_uncertainty: the free body where the F_x residual is systematic
UNCERTAINTY_BASES = ((10, 2, 2), (20, 4, 4), (30, 6, 6))
UNCERTAINTY_MODES = 128000                       # N of the fit, where each basis is converged in N for F_x
UNCERTAINTY_LINEARITY = 0.1   # the fit's largest residual must be below this times the estimated error


def drift_uncertainty(cyl=None):
    """The estimates behind KNOWN_UNCERTAINTY. For each omega of UNCERTAINTY_FREQUENCIES, free body, at N =
    UNCERTAINTY_MODES and each gap basis of UNCERTAINTY_BASES: the F_x-identity residual r = (near - far)/far, and F_z
    and M_y. Each of F_z, M_y is fitted as c + s r over the bases; the estimated error of the reference value q (at
    DEFAULT_RESOLUTION) is abs(q - c)/abs(q), the distance to the value at zero residual. The F_x near estimate is abs(r)
    at the default. Assumption: the error vanishes where the residual does (not testable here; the residual itself
    does not fall with the basis). Returns a list of dicts per (omega, quantity) with the residuals, the slope, the
    intercept, the fit's largest residual and the relative estimate."""
    cyl = Cylinder() if cyl is None else cyl
    res = DEFAULT_RESOLUTION
    out = []
    for omega in UNCERTAINTY_FREQUENCIES:
        rows = [_reference_row(cyl, omega, "free", UNCERTAINTY_MODES, _spec_key(b), _orders_key(res["orders"]))
                for b in UNCERTAINTY_BASES]
        ref = _reference_row(cyl, omega, "free", res["n_modes"], _spec_key(res["p_max"]), _orders_key(res["orders"]))
        r = np.array([(row["Fx_near"] - row["Fx_far"]) / row["Fx_far"] for row in rows])
        out.append(dict(omega=omega, quantity="Fx_near", residuals=r, slope=float("nan"), intercept=ref["Fx_far"],
                        fit=0.0, estimate=abs(ref["Fx_near"] - ref["Fx_far"]) / abs(ref["Fx_far"]), value=ref["Fx_near"]))
        for quantity in ("Fz", "My"):
            q = np.array([row[quantity] for row in rows])
            slope, intercept = np.polyfit(r, q, 1)
            fit = float(np.max(np.abs(q - (intercept + slope * r))))
            out.append(dict(omega=omega, quantity=quantity, residuals=r, slope=float(slope), intercept=float(intercept),
                            fit=fit, estimate=abs(ref[quantity] - intercept) / abs(ref[quantity]),
                            value=ref[quantity]))
    return out


def check_uncertainty():
    """drift_uncertainty: every estimate is within the stated uncertainty (uncertainty(), from KNOWN_UNCERTAINTY or
    DEFAULT_UNCERTAINTY), and every fit is linear to UNCERTAINTY_LINEARITY of its estimate (or to 1e-6 relative where
    the estimate is below 1e-5). It reproduces KNOWN_UNCERTAINTY; it does not test the zero-intercept assumption."""
    ok, lines = True, []
    for e in drift_uncertainty():
        res = DEFAULT_RESOLUTION
        stated = uncertainty("free", e["quantity"], e["omega"], Cylinder(), res["n_modes"], res["p_max"], res["orders"])
        linear = e["fit"] <= max(UNCERTAINTY_LINEARITY * e["estimate"], 1e-6) * abs(e["value"])
        ok &= e["estimate"] <= stated and linear
        lines.append("%.2f %-7s residuals %s slope %s intercept %.7g fit %.1e estimate %.1e (stated %.0e)%s"
                     % (e["omega"], e["quantity"], " ".join("%+.3e" % x for x in e["residuals"]),
                        "%.4g" % e["slope"] if e["quantity"] != "Fx_near" else "-", e["intercept"], e["fit"],
                        e["estimate"], stated, "" if linear else " NOT LINEAR"))
    return ok, ("free body, N = %d, gap bases %s: F_x residuals per basis, F_z and M_y fitted against them, estimated "
                "relative errors against the stated ones\n%s" % (UNCERTAINTY_MODES, ", ".join(str(b) for b in
                                                                                          UNCERTAINTY_BASES),
                                                                  "\n".join(lines)))


# Comparison with HAMS (compare; module docstring, "Comparison with HAMS")
COMPARE_PREFIX = "a_sharp_"   # the case folders compared: the sharp hull of group A of drift_convergence.study_cases
COMPARE_HEADING = 0.0
# The study's results, resolved from this script's folder (not the working directory): outside the repository
COMPARE_RESULTS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                                                "validation-runs", "snellius-drift", "results"))
IRREGULAR_MARGIN = 0.05   # rows within 5 % of an irregular frequency of the hull are marked "irregular"
MISSING = "missing"
COMPARE_PAIRS = (11, 15, 33, 55)   # A and B compared; 13 and 35 are zero by symmetry in the reference
COMPARE_FIRST = (["A%d" % ij for ij in COMPARE_PAIRS] + ["B%d" % ij for ij in COMPARE_PAIRS]
                 + ["X%d" % j for j in MOTION_MODES] + ["RAO%d" % j for j in MOTION_MODES])
COMPARE_DRIFT = ["Fx_near", "Fx_far", "Fz", "My", "Fx_cs"]
COMPARE_TERMS = ["%s_%s" % (name, key) for name in ("Fx", "Fz", "My") for key in NEAR_FIELD_TERMS]
COMPARE_QUANTITIES = COMPARE_FIRST + COMPARE_DRIFT + COMPARE_TERMS
_UNCERTAINTY_OF = {"Fx_near": "Fx_near", "Fx_far": "Fx_far", "Fz": "Fz", "My": "My", "Fx_cs": "Fx_far"}
# compare: two frequencies are the same if within this relative distance (HAMS writes omega with 7 digits, %14.6E)
FREQUENCY_MATCH = 1e-6


def _match_frequency(omega, candidates):
    """The candidate nearest omega if within FREQUENCY_MATCH of it (relative), else None."""
    best = min(candidates, key=lambda w: abs(w - omega), default=None)
    return best if best is not None and abs(best - omega) <= FREQUENCY_MATCH * abs(omega) else None
_TOTAL_OF = {"Fx": "Fx_near", "Fz": "Fz", "My": "My"}


def irregular_frequencies(cyl, upto):
    """The irregular frequencies of a constant-panel BEM (no lid) on the cylinder below upto, rad/s: the eigenvalues of the
    interior problem with phi = 0 on the wetted surface and the free-surface condition on the waterplane, phi =
    J_m(k r) sinh(k (z + d)) cos(m theta) with J_m(k a) = 0 and nu = k coth(k d). For the default cylinder 1.5483
    (J_0) and 1.9394 (J_1) below 2 rad/s, the 1.548 and 1.939 of drift_convergence.IRREGULAR_FREQUENCIES (that tuple
    also holds 1.10, 1.50 and 1.60, the frequencies of the study's irregular-frequency runs, which are not irregular)."""
    from scipy.special import jn_zeros

    found = []
    for m in range(8):
        for j in jn_zeros(m, 8):
            k = j / cyl.radius
            omega = math.sqrt(make_mesh.G * k / math.tanh(k * cyl.draft))
            if omega <= upto:
                found.append(omega)
    return tuple(sorted(found))


def polygon_area_error(ntheta):
    """Area of the HAMS waterline over pi a**2, less 1: make_mesh.cylinder puts the ntheta vertices on the circle (an
    inscribed regular polygon), area (ntheta/2) a**2 sin(2 pi/ntheta)."""
    return ntheta * math.sin(2.0 * math.pi / ntheta) / (2.0 * math.pi) - 1.0


def _read_drift_terms(path):
    """MeanDriftTerms.txt (written by write_terms_row in src/SecondOrderOutput.f90; header in src/SecondOrderDriver.f90):
    rows PER BETA1 BETA2 body I term Re Im, term 1 waterline, 2 flare, 3 velocity, 4 displacement, 5 rotation,
    6 hydrostatic (NEAR_FIELD_TERMS I, Ia, II, III, IV, V), I = 1 ... 6, scaled and conjugated as MeanDriftPressure.9.
    Returns (omega, beta1, beta2, I) -> {term: Re}, body 1, or None if the file is absent. drift_convergence.read_terms
    reads I = 1 only."""
    if not os.path.exists(path):
        return None
    values = {}
    for line in open(path):
        c = line.split()
        if not c or line.startswith("#") or int(c[3]) != 1:
            continue
        key = (round(float(c[0]), 4), round(float(c[1]), 3), round(float(c[2]), 3), int(c[4]))
        values.setdefault(key, {})[int(c[5])] = float(c[6])
    return values


def _read_hams_case(folder):
    """omega -> {quantity: HAMS value or MISSING} of output/Wamit_format/ under a case folder, heading COMPARE_HEADING, body 1.
    AmssDamp.1, ExcForce.3, Motion.4 by drift_convergence.read_hams_first_order (A/rho, B/(rho omega), MOD columns);
    MeanDriftMomentum.8 (F_x far), MeanDriftPressure.9 (index 6(b - 1) + I: F_x 1, F_z 3, M_y 5) and
    MeanDriftControlSurface.7 (F_x, same format) by drift_convergence.read_drift (Re column); MeanDriftTerms.txt by
    _read_drift_terms. A file that is absent makes its quantities MISSING."""
    out = os.path.join(folder, "output", "Wamit_format")
    exists = lambda name: os.path.exists(os.path.join(out, name))
    first = all(exists(n) for n in ("AmssDamp.1", "ExcForce.3", "Motion.4"))
    a, b, x, rao = drift_convergence.read_hams_first_order(out) if first else ({}, {}, {}, {})
    drift = {name: drift_convergence.read_drift(os.path.join(out, name)) if exists(name) else None
             for name in ("MeanDriftMomentum.8", "MeanDriftPressure.9", "MeanDriftControlSurface.7")}
    terms = _read_drift_terms(os.path.join(out, "MeanDriftTerms.txt"))
    beta, beta1 = round(COMPARE_HEADING, 3), round(COMPARE_HEADING, 1)
    # The first-order reader keys omega to 2 decimals, the drift readers to 4: each file's value is taken at its
    # frequency within FREQUENCY_MATCH of the row's (_match_frequency), MISSING if it has none
    x_omegas = {key[0] for key in x if key[1] == beta1}
    first_omegas = x_omegas | {key[0] for key in a} | {key[0] for key in rao if key[1] == beta1}
    drift_omegas = {name: {key[0] for key in (table or {}) if key[1] == beta and key[2] == beta}
                    for name, table in drift.items()}
    terms_omegas = {key[0] for key in (terms or {}) if key[1] == beta and key[2] == beta}
    omegas = []                       # the frequencies of ExcForce.3 and the drift files, merged within FREQUENCY_MATCH
    for omega in sorted(x_omegas.union(*drift_omegas.values())):
        if _match_frequency(omega, omegas) is None:
            omegas.append(omega)
    values = {}
    for omega in omegas:
        w1 = _match_frequency(omega, first_omegas)
        v = {}
        for ij in COMPARE_PAIRS:
            i, j = divmod(ij, 10)
            v["A%d" % ij] = a.get((w1, i, j), MISSING) if first else MISSING
            v["B%d" % ij] = b.get((w1, i, j), MISSING) if first else MISSING
        for j in MOTION_MODES:
            v["X%d" % j] = x.get((w1, beta1, j), MISSING)
            v["RAO%d" % j] = rao.get((w1, beta1, j), MISSING)
        for q, name, index in (("Fx_far", "MeanDriftMomentum.8", 1), ("Fx_near", "MeanDriftPressure.9", 1),
                               ("Fz", "MeanDriftPressure.9", 3), ("My", "MeanDriftPressure.9", 5),
                               ("Fx_cs", "MeanDriftControlSurface.7", 1)):
            table, wd = drift[name], _match_frequency(omega, drift_omegas[name])
            v[q] = MISSING if table is None else table.get((wd, beta, beta, index), MISSING)
        wt = _match_frequency(omega, terms_omegas)
        for name, index in (("Fx", 1), ("Fz", 3), ("My", 5)):
            row = None if terms is None else terms.get((wt, beta, beta, index))
            for number, key in enumerate(NEAR_FIELD_TERMS, start=1):
                v["%s_%s" % (name, key)] = MISSING if row is None or number not in row else row[number]
        values[omega] = v
    return values


def _compare_relative(value, ref, scale):
    """(value - ref)/abs(ref); None where drift_convergence.percent prints an empty field (abs(ref) below 1e-4 of scale),
    MISSING where HAMS has no value."""
    if value == MISSING:
        return MISSING
    if scale <= 0.0 or drift_convergence.percent(value, ref, scale) == "":
        return None
    return (value - ref) / abs(ref)


def _reference_of(q):
    """The reference column a compared quantity is measured against: Fx_cs against Fx_far (the control surface is a
    far-field method), the others against their own."""
    return "Fx_far" if q == "Fx_cs" else q


def _floor_scales(by):
    """(motion, quantity) -> the largest abs(reference) of the quantity and motion over the frequencies of the reference
    rows, the scale of drift_convergence.percent's floor; for a near-field term, the largest of the term and its total
    (Fx_near, Fz, My), so that a term is empty only where it is negligible against the force it belongs to. RAO uses the
    free rows for both motions."""
    scales = {}
    for motion in ("fixed", "free"):
        for q in COMPARE_QUANTITIES:
            source = "free" if q.startswith("RAO") else motion
            cols = [_reference_of(q)]
            if q in COMPARE_TERMS:
                cols.append(_TOTAL_OF[q.split("_")[0]])
            vals = [abs(row[c]) for (m, _), row in by.items() if m == source for c in cols]
            scales[motion, q] = max(vals) if vals else 0.0
    return scales


def _sharp_cases():
    return {case.name: case for case in drift_convergence.study_cases() if case.name.startswith(COMPARE_PREFIX)}


def compare(results, reference_rows, names=None):
    """One dict per (case, omega) at heading 0 for every folder of results whose name starts with COMPARE_PREFIX (or the
    folders names), against the reference rows (reference's dicts, fixed and free, any frequencies; the RAOs of both
    motions are compared with the free rows, since HAMS writes Motion.4 for the fixed body too). Keys: case, motion,
    velocity, multiple, panels (drift_convergence.panels_of), ntheta, area_error (polygon_area_error), omega, flag
    ("irregular" within IRREGULAR_MARGIN of irregular_frequencies, "no reference" if the reference has no row within
    FREQUENCY_MATCH of omega, else ""), and for each q of COMPARE_QUANTITIES: q (HAMS value or MISSING), q_ref, q_rel
    ((HAMS - reference)/abs(reference), None below the floor, MISSING), and for COMPARE_DRIFT q_unc (the stated relative
    uncertainty that the reference row carries, its *_unc column; nan, not measured, unless the rows are the default
    hull at DEFAULT_RESOLUTION). Cases are sorted by motion, velocity, multiple. A COMPARE_PREFIX folder that is not a
    case of drift_convergence.study_cases, or has no output/Wamit_format, is skipped with a printed note."""
    cyl = Cylinder()
    by = {(r["motion"], r["omega"]): r for r in reference_rows}
    ref_omegas = {motion: [w for (m, w) in by if m == motion] for motion in ("fixed", "free")}
    scales = _floor_scales(by)
    cases = _sharp_cases()
    folders = []
    for n in sorted(names if names is not None else os.listdir(results)):
        if not n.startswith(COMPARE_PREFIX) or not os.path.isdir(os.path.join(results, n)):
            continue
        if n not in cases:
            print("compare: skipped %s: not a case of drift_convergence.study_cases" % n)
        elif not os.path.isdir(os.path.join(results, n, "output", "Wamit_format")):
            print("compare: skipped %s: no output/Wamit_format" % n)
        else:
            folders.append(n)
    folders.sort(key=lambda n: (cases[n].motion, cases[n].velocity, cases[n].multiple))
    top = max(drift_convergence.WIDE_FREQUENCIES) * (1.0 + IRREGULAR_MARGIN) / (1.0 - IRREGULAR_MARGIN)
    irregular = irregular_frequencies(cyl, top)
    if len(irregular) != 2:   # 1.5483 and 1.9394 rad/s on the default cylinder; IRREGULAR_MARGIN was set for these two
        raise ValueError("compare: expected 2 irregular frequencies below %.3f rad/s, found %s" % (top, irregular))
    table = []
    for name in folders:
        case = cases[name]
        ntheta = drift_convergence.scaled(drift_convergence.fixtures.DRIFT_CYL_MESH, ("ntheta",), case.multiple)["ntheta"]
        for omega, values in _read_hams_case(os.path.join(results, name)).items():
            row = dict(case=name, motion=case.motion, velocity=case.velocity, multiple=case.multiple,
                       panels=drift_convergence.panels_of(case), ntheta=ntheta, area_error=polygon_area_error(ntheta),
                       omega=omega, flag="")
            if any(abs(omega - w) <= IRREGULAR_MARGIN * w for w in irregular):
                row["flag"] = "irregular"
            ref = by.get((case.motion, _match_frequency(omega, ref_omegas[case.motion])))
            free = by.get(("free", _match_frequency(omega, ref_omegas["free"])))
            if ref is None or free is None:
                row["flag"] = "no reference"
            for q in COMPARE_QUANTITIES:
                source = free if q.startswith("RAO") else ref
                row[q] = values[q]
                if source is None:
                    row[q + "_ref"], row[q + "_rel"] = float("nan"), MISSING
                else:
                    row[q + "_ref"] = source[_reference_of(q)]
                    row[q + "_rel"] = _compare_relative(values[q], row[q + "_ref"], scales[case.motion, q])
                if q in COMPARE_DRIFT:   # the stated uncertainty of the reference row itself (nan: not measured)
                    row[q + "_unc"] = (float("nan") if ref is None
                                       else ref.get(_UNCERTAINTY_OF[q] + "_unc", float("nan")))
            table.append(row)
    return table


def _pct_text(rel):
    if rel == MISSING:
        return MISSING
    return "" if rel is None else "%+.3f" % (100.0 * rel)


def format_compare(table):
    """The comparison as text: one block per case, headed by its panels and polygon area error, with three tables, each a
    row per omega: first order, drift (each percentage followed by the reference's stated uncertainty in percent, [u],
    or [n/m] where the reference row states none),
    and the near-field terms (Ia, zero on the vertical wall in both, is left out here; it is in the CSV)."""
    lines = []
    shown_terms = [q for q in COMPARE_TERMS if not q.endswith("_Ia")]
    for name in sorted(set(r["case"] for r in table), key=[r["case"] for r in table].index):
        rows = [r for r in table if r["case"] == name]
        r0 = rows[0]
        lines.append("")
        lines.append("%s: %s, %s, multiple %g, %d panels, waterline a %d-sided inscribed polygon, area/(pi a**2) - 1 = "
                     "%+.3e" % (name, r0["motion"], r0["velocity"], r0["multiple"], r0["panels"], r0["ntheta"],
                                r0["area_error"]))
        lines.append("  first order, % difference (HAMS - reference)/|reference| (RAO against the free body):")
        lines.append("  %5s %-9s %s" % ("omega", "flag", " ".join("%8s" % q for q in COMPARE_FIRST)))
        for r in rows:
            lines.append("  %5.2f %-9s %s" % (r["omega"], r["flag"],
                                               " ".join("%8s" % _pct_text(r[q + "_rel"]) for q in COMPARE_FIRST)))
        lines.append("  drift, % difference [reference uncertainty, %]; Fx_cs against far-field F_x of the reference:")
        lines.append("  %5s %-9s %s" % ("omega", "flag", " ".join("%17s" % q for q in COMPARE_DRIFT)))
        for r in rows:
            cells = []
            for q in COMPARE_DRIFT:
                text = _pct_text(r[q + "_rel"])
                unc = r[q + "_unc"]
                cells.append("%17s" % (text if text in ("", MISSING) else "%s [%s]" % (
                    text, "n/m" if math.isnan(unc) else "%.2f" % (100.0 * unc))))
            lines.append("  %5.2f %-9s %s" % (r["omega"], r["flag"], " ".join(cells)))
        lines.append("  near-field terms (MeanDriftTerms.txt), % difference:")
        lines.append("  %5s %-9s %s" % ("omega", "flag", " ".join("%8s" % q for q in shown_terms)))
        for r in rows:
            lines.append("  %5.2f %-9s %s" % (r["omega"], r["flag"],
                                               " ".join("%8s" % _pct_text(r[q + "_rel"]) for q in shown_terms)))
    return "\n".join(lines)


def compare_columns():
    """The columns of the comparison CSV, in order."""
    cols = ["case", "motion", "velocity", "multiple", "panels", "ntheta", "area_error", "omega", "flag"]
    for q in COMPARE_QUANTITIES:
        cols += [q, q + "_ref", q + "_pct"] + ([q + "_unc_pct"] if q in COMPARE_DRIFT else [])
    return cols


def write_compare_csv(path, table):
    """The comparison as CSV (compare_columns): values with 10 digits, percentages with 3 decimals as
    drift_convergence.percent, empty below the floor, MISSING where HAMS has no value."""
    cols = compare_columns()
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(cols)
        for r in table:
            out = []
            for c in cols:
                if c.endswith("_unc_pct"):   # empty where the reference row states none (nan, not measured)
                    unc = r[c[:-8] + "_unc"]
                    out.append("" if math.isnan(unc) else "%.3g" % (100.0 * unc))
                elif c.endswith("_pct"):
                    out.append(_pct_text(r[c[:-4] + "_rel"]))
                elif c == "area_error":
                    out.append("%.6e" % r[c])
                else:
                    out.append(_csv_value(r[c]) if not isinstance(r[c], str) else r[c])
            writer.writerow(out)


def default_reference_rows():
    """The fixed and free rows of reference at DEFAULT_RESOLUTION and every WIDE_FREQUENCIES (cached in _reference_row)."""
    cyl, res = Cylinder(), DEFAULT_RESOLUTION
    rows = []
    for omega in drift_convergence.WIDE_FREQUENCIES:          # omega outermost: fixed and free share the solves
        for motion in ("fixed", "free"):
            rows += reference(cyl, [omega], motion, res["n_modes"], res["p_max"], res["orders"])
    return rows


def run_compare(results, out=None, reference_csv=None):
    """The compare subcommand: the reference rows (read from reference_csv, a CSV of the reference subcommand, or computed
    by default_reference_rows), compare, the text table, and the CSV if out is given. Returns 1 if results holds no
    COMPARE_PREFIX folder, else 0."""
    start = time.time()
    if not os.path.isdir(results):
        print("compare: %s is not a folder" % results)
        return 1
    if reference_csv:
        rows = read_reference_csv(reference_csv)
        source = "read from %s" % reference_csv
    else:
        rows = default_reference_rows()
        source = "computed"
    resolutions = sorted(set((r["n_modes"], str(r["p_max"]), r["orders"]) for r in rows))
    source += ", resolution %s" % "; ".join("N %d, gap basis %s, orders 0 ... %d" % res for res in resolutions)
    if not all(_is_default_resolution(r["n_modes"], r["p_max"], r["orders"]) for r in rows):
        source += " (not DEFAULT_RESOLUTION)"
    unstated = sum(math.isnan(r["%s_unc" % q]) for r in rows for q in UNCERTAIN_QUANTITIES)
    if unstated:
        source += "; %d stated uncertainties are nan (not measured: only the default hull at DEFAULT_RESOLUTION " \
                  "has them), shown [n/m]" % unstated
    table = compare(results, rows)
    if not table:
        print("compare: no %s* case folder with results in %s" % (COMPARE_PREFIX, results))
        return 1
    irregular = irregular_frequencies(Cylinder(), 2.0)
    print("HAMS against the reference (%s), heading %g, %d cases, %d rows. Irregular frequencies of the hull below "
          "2 rad/s: %s; rows within %g %% are marked irregular. Empty: |reference| below 1e-4 of the largest of its "
          "quantity and motion (drift_convergence.percent)."
          % (source, COMPARE_HEADING, len(set(r["case"] for r in table)), len(table),
             ", ".join("%.4f" % w for w in irregular), 100 * IRREGULAR_MARGIN))
    print(format_compare(table))
    if out:
        write_compare_csv(out, table)
        print("wrote %d rows to %s" % (len(table), out))
    print("compare: %.0f s" % (time.time() - start))
    return 0


COMPARE_FIXTURE_FREQUENCIES = (0.40, 0.70, 1.60)   # free F_x below the floor, an ordinary row, a row near 1.548 rad/s


def _fixture_value(value):
    return "%23.15E" % value


# The fixture writes every complex value with this phase, so that a reader that takes the wrong column (Re or Im for MOD,
# MOD or Im for Re) sees a relative difference of at least 1 - cos(phase) = 0.17
FIXTURE_PHASE = 0.6   # rad
FIXTURE_UNCERTAINTY = 0.0123   # check_compare_fixture: Fz_unc set on the free 0.70 rad/s row; compare must carry it


def _fixture_complex(mod=None, re=None):
    """MOD, PHS (degrees), Re, Im of a WAMIT-format row with Re = MOD cos(phase), Im = MOD sin(phase): from its modulus
    (first-order files, phase FIXTURE_PHASE) or from its real part (drift files, phase FIXTURE_PHASE or FIXTURE_PHASE +
    pi, MOD = abs(re)/cos(FIXTURE_PHASE))."""
    if mod is not None:
        phase = FIXTURE_PHASE
    else:
        mod = abs(re) / math.cos(FIXTURE_PHASE)
        phase = FIXTURE_PHASE + (math.pi if re < 0.0 else 0.0)
    return mod, math.degrees(phase), mod * math.cos(phase), mod * math.sin(phase)


def _write_hams_fixture(folder, rows, free_rows, terms=True, control_surface=True):
    """Write output/Wamit_format/ of a HAMS case under folder from reference rows (omega -> row), in the formats HAMS writes
    (see _read_hams_case), with 15 digits. free_rows give Motion.4. Rows at heading 30 hold other values, which compare
    must not read. Every complex value has the phase FIXTURE_PHASE (_fixture_complex), so the readers' columns are tested.
    terms and control_surface: whether MeanDriftTerms.txt and MeanDriftControlSurface.7 are written."""
    out = os.path.join(folder, "output", "Wamit_format")
    os.makedirs(out)
    v = _fixture_value
    with open(os.path.join(out, "AmssDamp.1"), "w") as f:
        for omega, row in sorted(rows.items()):
            for ij in REFERENCE_PAIRS:
                i, j = divmod(ij, 10)
                f.write("%14.6E %5d %5d %s %s\n" % (omega, i, j, v(row["A%d" % ij]), v(row["B%d" % ij])))
    for name, source, key in (("ExcForce.3", rows, "X"), ("Motion.4", free_rows, "RAO")):
        with open(os.path.join(out, name), "w") as f:
            for omega in sorted(rows):
                for beta, factor in ((0.0, 1.0), (30.0, 2.0)):
                    for j in MOTION_MODES:
                        cols = _fixture_complex(mod=factor * source[omega]["%s%d" % (key, j)])
                        f.write("%14.6E %14.6E %5d %s\n" % (omega, beta, j, " ".join(v(c) for c in cols)))

    def drift(name, entries):
        with open(os.path.join(out, name), "w") as f:
            for omega, row in sorted(rows.items()):
                for beta, factor in ((0.0, 1.0), (30.0, 2.0)):
                    for index, key in entries:
                        cols = _fixture_complex(re=factor * row[key])
                        f.write("%14.6E %14.6E %14.6E %5d %s\n" % (omega, beta, beta, index, " ".join(v(c) for c in cols)))

    drift("MeanDriftMomentum.8", ((1, "Fx_far"),))
    drift("MeanDriftPressure.9", ((1, "Fx_near"), (3, "Fz"), (5, "My")))
    if control_surface:
        drift("MeanDriftControlSurface.7", ((1, "Fx_far"),))
    if terms:
        with open(os.path.join(out, "MeanDriftTerms.txt"), "w") as f:
            f.write("# PER BETA1 BETA2 body I term Re Im, scaled and conjugated as MeanDriftPressure.9\n")
            f.write("# term: 1 waterline, 2 flare, 3 velocity, 4 displacement, 5 rotation, 6 hydrostatic\n")
            for omega, row in sorted(rows.items()):
                for beta, factor in ((0.0, 1.0), (30.0, 2.0)):
                    for index, name in ((1, "Fx"), (3, "Fz"), (5, "My")):
                        for number, key in enumerate(NEAR_FIELD_TERMS, start=1):
                            _, _, re, im = _fixture_complex(re=factor * row["%s_%s" % (name, key)])
                            f.write("%14.6E %14.6E %14.6E %5d %5d %5d %s %s\n"
                                    % (omega, beta, beta, 1, index, number, v(re), v(im)))


def check_compare_fixture():
    """compare on a temporary RESULTS written from the reference itself (LIGHT_RESOLUTION, COMPARE_FIXTURE_FREQUENCIES), the
    rows passed through the reference CSV (write_reference_csv, read_reference_csv): a free case with every file, where
    every difference is 0 to 1e-9 and the control surface is read; a fixed case without MeanDriftTerms.txt and
    MeanDriftControlSurface.7, whose columns are "missing" while the rest is 0; every complex value written with the phase
    FIXTURE_PHASE, so a reader that takes the wrong column fails; an a_sharp_ folder that is not a study case and one
    without output/Wamit_format, skipped; the 1.60 rad/s rows marked "irregular" and
    the others not; free F_x near and far at 0.40 rad/s (3e-6 of 0.11) empty, below the floor of
    drift_convergence.percent; the reference rows (LIGHT_RESOLUTION, not the default) all nan in *_unc, and compare's
    q_unc taken from the row (nan everywhere but FIXTURE_UNCERTAINTY, set on the free 0.70 rad/s row's Fz_unc); and the
    irregular frequencies of irregular_frequencies equal to the 1.548 and 1.939 rad/s of
    drift_convergence.IRREGULAR_FREQUENCIES to 1e-3."""
    import tempfile

    cyl, res = Cylinder(), LIGHT_RESOLUTION
    computed = []
    for omega in COMPARE_FIXTURE_FREQUENCIES:
        for motion in ("fixed", "free"):
            computed += reference(cyl, [omega], motion, res["n_modes"], res["p_max"], res["orders"])
    problems = []
    with tempfile.TemporaryDirectory() as root:
        path = os.path.join(root, "reference.csv")
        write_reference_csv(path, computed)
        rows = read_reference_csv(path)
        by = {(r["motion"], round(r["omega"], 2)): r for r in rows}
        free = {w: by["free", w] for w in COMPARE_FIXTURE_FREQUENCIES}
        fixed = {w: by["fixed", w] for w in COMPARE_FIXTURE_FREQUENCIES}
        results = os.path.join(root, "results")
        _write_hams_fixture(os.path.join(results, "a_sharp_free_src_m2"), free, free)
        _write_hams_fixture(os.path.join(results, "a_sharp_fixed_diff_m1"), fixed, free, terms=False,
                            control_surface=False)
        os.makedirs(os.path.join(results, "b_bottom_src_m1"))   # not a_sharp_: not read
        os.makedirs(os.path.join(results, "a_sharp_free_src_m9", "output", "Wamit_format"))   # not a study case: skipped
        os.makedirs(os.path.join(results, "a_sharp_free_src_m1"))   # no output/Wamit_format: skipped
        stated = sum(not math.isnan(r["%s_unc" % q]) for r in rows for q in UNCERTAIN_QUANTITIES)
        if stated:   # LIGHT_RESOLUTION is not the default: uncertainty states none
            problems.append("%d stated uncertainties in the LIGHT_RESOLUTION rows" % stated)
        by["free", 0.70]["Fz_unc"] = FIXTURE_UNCERTAINTY   # compare must carry the row's own *_unc
        table = compare(results, rows)
    cases = sorted(set(r["case"] for r in table))
    if cases != ["a_sharp_fixed_diff_m1", "a_sharp_free_src_m2"]:
        problems.append("cases %s" % cases)
    worst, empty, missing_ok = 0.0, [], True
    for r in table:
        for q in COMPARE_QUANTITIES:
            rel = r[q + "_rel"]
            no_files = r["case"] == "a_sharp_fixed_diff_m1" and (q == "Fx_cs" or q in COMPARE_TERMS)
            if no_files:
                missing_ok &= rel == MISSING and r[q] == MISSING
            elif rel == MISSING:
                problems.append("%s %.2f %s missing" % (r["case"], r["omega"], q))
            elif rel is None:
                empty.append((r["case"], r["omega"], q))
            else:
                worst = max(worst, abs(rel))
        for q in COMPARE_DRIFT:
            marked = r["case"] == "a_sharp_free_src_m2" and r["omega"] == 0.70 and q == "Fz"
            if marked != (r[q + "_unc"] == FIXTURE_UNCERTAINTY) or not (marked or math.isnan(r[q + "_unc"])):
                problems.append("%s %.2f %s_unc %r" % (r["case"], r["omega"], q, r[q + "_unc"]))
        flagged = r["flag"] == "irregular"
        if flagged != (r["omega"] == 1.60):
            problems.append("%s %.2f flag %r" % (r["case"], r["omega"], r["flag"]))
    if worst > 1e-9:
        problems.append("largest difference %.1e" % worst)
    if not missing_ok:
        problems.append("the fixed case's terms and control surface are not all 'missing'")
    for q in ("Fx_near", "Fx_far", "Fx_cs"):
        if ("a_sharp_free_src_m2", 0.40, q) not in empty:
            problems.append("free %s at 0.40 rad/s not empty" % q)
    irregular = irregular_frequencies(cyl, max(drift_convergence.WIDE_FREQUENCIES) * 1.2)
    for listed in (1.548, 1.939):
        if listed not in drift_convergence.IRREGULAR_FREQUENCIES or min(abs(w - listed) for w in irregular) > 1e-3:
            problems.append("irregular frequency %.3f: computed %s" % (listed, irregular))
    return not problems, ("%d rows, largest difference %.1e, %d empty fields (below the floor: %s), irregular "
                          "frequencies %s; %s" % (len(table), worst, len(empty),
                                                  ", ".join(sorted(set("%s %.2f" % (e[0][8:], e[1]) for e in empty))),
                                                  ", ".join("%.4f" % w for w in irregular),
                                                  "; ".join(problems) or "no problems"))


SIGN_CASES = ("a_sharp_free_src_m4", "a_sharp_fixed_src_m4")
SIGN_QUANTITIES = ("Fx_near", "Fx_far", "Fz", "My")
# compare_signs: a near-field total is sign-tested only where abs(reference) is at least this fraction of its largest
# reference term. HAMS's terms at 7168 panels are off by up to 3.4 % (free F_z II at 0.70 rad/s) and more (F_x II, 12 %),
# so a total that cancels its terms to below 5 % can take either sign without any convention error. Its terms are tested.
SIGN_CANCELLATION = 0.05


def check_compare_signs():
    """At the finest source mesh of the study (SIGN_CASES, if COMPARE_RESULTS holds a_sharp_free_src_m4; else SKIP:
    ok is None, which _run_checks prints as SKIP and counts as neither passed nor failed), outside the irregular rows
    and wherever abs(reference) is above the floor of drift_convergence.percent: every near-field term (MeanDriftTerms.txt) of F_x, F_z and M_y, and far-field F_x, has the reference's sign; and the
    near-field totals F_x, F_z, M_y do too where abs(reference) is at least SIGN_CANCELLATION of the largest abs(term) of
    the reference (the totals below that are listed, not tested). A flipped sign or a conjugate in a conversion would
    flip whole terms. COMPARE_RESULTS is resolved from this script's folder, not the working directory."""
    if not os.path.isdir(os.path.join(COMPARE_RESULTS, SIGN_CASES[0])):
        return None, "%s not found" % os.path.join(COMPARE_RESULTS, SIGN_CASES[0])
    table = compare(COMPARE_RESULTS, default_reference_rows(),
                    names=[n for n in SIGN_CASES if os.path.isdir(os.path.join(COMPARE_RESULTS, n))])
    tested, wrong, cancelling = 0, [], []
    for r in table:
        if r["flag"] == "irregular":
            continue
        quantities = list(SIGN_QUANTITIES) + COMPARE_TERMS
        for q in quantities:
            if r[q + "_rel"] in (None, MISSING):
                continue
            prefix = {"Fx_near": "Fx", "Fz": "Fz", "My": "My"}.get(q)
            if prefix:
                largest = max(abs(r["%s_%s_ref" % (prefix, key)]) for key in NEAR_FIELD_TERMS)
                if abs(r[q + "_ref"]) < SIGN_CANCELLATION * largest:
                    cancelling.append("%s %.2f %s %.4g (%.1f %% of its largest term; HAMS %.4g)"
                                      % (r["case"][8:], r["omega"], q, r[q + "_ref"],
                                         100.0 * abs(r[q + "_ref"]) / largest, r[q]))
                    continue
            tested += 1
            if r[q] * r[q + "_ref"] <= 0.0:
                wrong.append("%s %.2f %s HAMS %.4g reference %.4g" % (r["case"], r["omega"], q, r[q], r[q + "_ref"]))
    return not wrong, ("%d values tested (totals and terms) in %s; opposite sign: %s. Totals not tested (cancelling): %s"
                       % (tested, ", ".join(sorted(set(r["case"] for r in table))), "; ".join(wrong) or "none",
                          "; ".join(cancelling) or "none"))


CHECKS = {"roots": check_roots, "inner_products": check_inner_products, "polylog": check_polylog,
          "matching_tail": check_matching_tail, "bottom_limit_x1": check_bottom_limit_x1, "matching": check_matching,
          "scattering_energy": check_scattering_energy, "reciprocity": check_reciprocity,
          "radiation_energy": check_radiation_energy, "particular_bottom": check_particular_bottom,
          "bottom_limit_drift": check_bottom_limit_drift, "body_matrices": check_body_matrices,
          "free_energy": check_free_energy, "long_wave_motions": check_long_wave_motions,
          "corner_limit": check_corner_limit, "reference_point": check_reference_point,
          "pair_symmetry": check_pair_symmetry, "fx_identity": check_fx_identity, "convergence": check_convergence,
          "compare_fixture": check_compare_fixture, "qtf_diagonal": check_qtf_diagonal}
# Run by name (check uncertainty compare_signs) and by the reference subcommand, not by a bare check: uncertainty computes
# nine rows at N = 128000 that no other check needs; compare_signs reads the study's results outside the repository and
# computes the 22 reference rows (about 92 s)
SLOW_CHECKS = {"uncertainty": check_uncertainty, "compare_signs": check_compare_signs}
ALL_CHECKS = dict(CHECKS, **SLOW_CHECKS)


def reference_columns():
    """The columns of the reference CSV, in order (reference's keys). A tuple (p_max) is written space-separated,
    "20 4 4"; nan (the fixed body's RAOs) as an empty field."""
    cols = ["omega", "motion", "n_modes", "p_max", "orders"]
    cols += ["A%d" % ij for ij in REFERENCE_PAIRS] + ["B%d" % ij for ij in REFERENCE_PAIRS]
    cols += ["X%d" % j for j in MOTION_MODES] + ["RAO%d" % j for j in MOTION_MODES]
    cols += ["Fx_far", "Fx_near", "Fz", "My"]
    cols += ["%s_unc" % q for q in UNCERTAIN_QUANTITIES]
    cols += ["%s_%s" % (name, key) for name in ("Fx", "Fz", "My") for key in NEAR_FIELD_TERMS]
    return cols


def _csv_value(value):
    if isinstance(value, tuple):
        return " ".join(str(v) for v in value)
    if isinstance(value, float):
        return "" if math.isnan(value) else "%.10g" % value
    return str(value)


def write_reference_csv(path, rows):
    """The reference rows as CSV, columns reference_columns (_csv_value)."""
    cols = reference_columns()
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(cols)
        for row in rows:
            writer.writerow([_csv_value(row[c]) for c in cols])


def read_reference_csv(path):
    """The rows of a CSV of write_reference_csv, as reference's dicts: motion a string, n_modes and orders ints, p_max an
    int or a tuple ("20 4 4"), the rest floats (an empty field nan). The values carry the CSV's 10 digits."""
    rows = []
    with open(path, newline="") as handle:
        for record in csv.DictReader(handle):
            row = {}
            for key, text in record.items():
                if key == "motion":
                    row[key] = text
                elif key == "p_max":
                    parts = [int(p) for p in text.split()]
                    row[key] = parts[0] if len(parts) == 1 else tuple(parts)
                elif key in ("n_modes", "orders"):
                    row[key] = int(text)
                else:
                    row[key] = float(text) if text != "" else float("nan")
            rows.append(row)
    return rows


def format_reference(rows):
    """The main reference values as a text table (every column is in the CSV)."""
    cols = ["A11", "A15", "A33", "A55", "B11", "B15", "B33", "B55", "X1", "X3", "X5", "RAO1", "RAO3", "RAO5",
            "Fx_far", "Fx_near", "Fz", "My"]
    lines = ["%-6s %5s %s" % ("motion", "omega", " ".join("%12s" % c for c in cols))]
    for row in rows:
        lines.append("%-6s %5.2f %s" % (row["motion"], row["omega"],
                                         " ".join("%12s" % ("" if math.isnan(row[c]) else "%.6g" % row[c])
                                                  for c in cols)))
    return "\n".join(lines)


def _run_checks(names):
    """Run the checks; (failed, skipped). A check returns ok True (PASS), False (FAIL) or None (SKIP: it could not
    run for want of data, as compare_signs without the study's results)."""
    failed = skipped = 0
    for name in names:
        start = time.time()
        try:
            ok, message = ALL_CHECKS[name]()
        except Exception as exc:  # a check that cannot run has failed
            ok, message = False, "%s: %s" % (type(exc).__name__, exc)
        failed += ok is not None and not ok  # a check may return numpy.bool_, which "ok is False" misses
        skipped += ok is None
        status = "SKIP" if ok is None else "PASS" if ok else "FAIL"
        print("%s %s (%.0f s): %s" % (status, name, time.time() - start, message), flush=True)
    return failed, skipped


def run_reference(out=None):
    """The reference subcommand: the fixed and free rows at DEFAULT_RESOLUTION for every WIDE_FREQUENCIES (printed, and
    written to out as CSV if given), then the convergence table (the convergence check), then every other check.
    Returns the number of failed checks."""
    start = time.time()
    cyl, res = Cylinder(), DEFAULT_RESOLUTION
    rows = {}
    for omega in drift_convergence.WIDE_FREQUENCIES:          # omega outermost: fixed and free share the solves
        for motion in ("fixed", "free"):
            rows[motion, omega] = reference(cyl, [omega], motion, res["n_modes"], res["p_max"], res["orders"])[0]
    ordered = [rows[motion, omega] for motion in ("fixed", "free") for omega in drift_convergence.WIDE_FREQUENCIES]
    print("Reference values at DEFAULT_RESOLUTION (N %d, gap basis %s, orders 0 ... %d); A/rho, B/(rho omega), "
          "abs(X)/(rho g A), abs(xi) (rad/m for 5), F/(rho g A**2), M_y in m about the origin:"
          % (res["n_modes"], res["p_max"], res["orders"]))
    print(format_reference(ordered))
    print(FREE_BODY_UNCERTAINTY)
    if out:
        write_reference_csv(out, ordered)
        print("wrote %d rows to %s" % (len(ordered), out))
    print("reference values: %.0f s" % (time.time() - start), flush=True)
    failed, skipped = _run_checks(["convergence"] + [name for name in ALL_CHECKS if name != "convergence"])
    print("reference: %d check(s) failed, %d skipped; run time %.0f s" % (failed, skipped, time.time() - start))
    return failed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p_check = sub.add_parser("check", help="run the named checks, or all of them")
    p_check.add_argument("names", nargs="*", metavar="NAME")
    p_ref = sub.add_parser("reference", help="reference values, the convergence table and all checks")
    p_ref.add_argument("--out", metavar="CSV", help="write the reference rows to this CSV file")
    p_qtf = sub.add_parser("qtf", help="quadratic QTF of the near-field terms on the grid of the Nemoh runs")
    p_qtf.add_argument("--out", metavar="CSV", required=True)
    p_qtf.add_argument("--motion", choices=("fixed", "free", "both"), default="both")
    p_cmp = sub.add_parser("compare", help="compare the HAMS a_sharp_* case folders under RESULTS with the reference")
    p_cmp.add_argument("results", metavar="RESULTS", help="folder holding the case folders (a_sharp_*/output/Wamit_format)")
    p_cmp.add_argument("--out", metavar="CSV", help="write the comparison to this CSV file")
    p_cmp.add_argument("--reference", metavar="CSV", help="read the reference rows from a CSV of the reference subcommand "
                                                          "instead of computing them (93 s); under 1 s")
    args = parser.parse_args(argv)
    if args.command == "reference":
        return 1 if run_reference(args.out) else 0
    if args.command == "compare":
        return run_compare(args.results, args.out, args.reference)
    if args.command == "qtf":
        return run_qtf(args.out, ("fixed", "free") if args.motion == "both" else (args.motion,))
    names = args.names or list(CHECKS)
    unknown = [n for n in names if n not in ALL_CHECKS]
    if unknown:
        parser.error("unknown check(s) %s; known: %s" % (", ".join(unknown), ", ".join(ALL_CHECKS)))
    start = time.time()
    failed, skipped = _run_checks(names)
    print("check: %d of %d failed, %d skipped; run time %.0f s" % (failed, len(names), skipped, time.time() - start))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
