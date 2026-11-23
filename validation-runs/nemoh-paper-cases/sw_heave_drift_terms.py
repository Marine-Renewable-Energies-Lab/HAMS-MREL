"""Desk check of the SOFTWIND heave mean drift (QTF-_3 diagonal, heading 0) near the pitch resonance, where HAMS-MREL and
HydroStar have opposite signs at 0.314 rad/s. Splits HAMS's diagonal into its near-field terms (MeanDriftTerms.txt),
recomputes the rotation and hydrostatic terms by hand from the first-order motions, and fits the difference
HydroStar - HAMS to those two terms over all frequencies.

Finding: the difference is HAMS's term V. The Nemoh 3.0.3 documentation (input-output.rst, QTF flag "hydrostatic terms of
the quadratic first order motion") states that HydroStar does not include -[K] xi2 with xi2_3 = z_G (theta_1^2 + theta_2^2)/2,
which for heave is HAMS's term V (Lee's F_S with Z_o = z_G); the test case runs Nemoh with that flag at 0 as well.

Usage: python -I sw_heave_drift_terms.py HAMS_WAMIT_FORMAT_DIR HYDROSTAR_DIR OUT_DIR
Writes sw_heave_drift_terms.png and .txt. Values are QTF/(rho g), per unit amplitude squared (HAMS's scaling).
Hand formulas (reference point at the CoG, heading 0, so alpha_1 = 0 and the surge-pitch coupling of the mass matrix is
zero): rotation IV_z = -1/2 Re(alpha_2 conj(G_1)) with G_1 = -omega^2 m xi_1 from the equation of motion (C_1j = 0), so
IV_z = 1/2 (omega^2 / g) V Re(xi_5 conj(xi_1)) with V = m / rho; hydrostatic V_z = -1/4 |xi_5|^2 A_wp z_o.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import readers as R  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

G = 9.80665
NAMES = {1: "I waterline", 2: "Ia flare", 3: "II velocity", 4: "III displacement", 5: "IV rotation", 6: "V hydrostatic"}


def motions(d):
    a = np.loadtxt(os.path.join(d, "Motion.4"))
    out = {}
    for i in (1, 3, 5):
        s = a[(a[:, 2] == i) & (a[:, 1] == 0.0)]
        out[i] = (s[:, 0], s[:, 5] + 1j * s[:, 6])
    return out


def terms(d):
    a = np.loadtxt(os.path.join(d, "MeanDriftTerms.txt"))
    a = a[(a[:, 1] == 0) & (a[:, 2] == 0) & (a[:, 3] == 1) & (a[:, 4] == 3)]
    w = np.unique(a[:, 0])
    return w, {t: np.array([a[(a[:, 0] == x) & (a[:, 5] == t), 6][0] for x in w]) for t in range(1, 7)}


def main(hdir, rdir, out):
    w, T = terms(hdir)
    grid = np.round(w, 4)
    ham = np.real(np.diag(R.hams_qtf(os.path.join(hdir, "QTF_Indirect.11d"), 3, 0.0, 0.0, grid)))
    q = R.hs_qtf(os.path.join(rdir, "MQTFFz.rao"))
    ref = np.real(np.diag(R.hs_pair_matrix(q, (0.0, 0.0), grid))) / (R.RHO_HS * R.G_HS)
    mo = motions(hdir)
    xi1 = np.interp(w, mo[1][0], mo[1][1].real) + 1j * np.interp(w, mo[1][0], mo[1][1].imag)
    xi5 = np.interp(w, mo[5][0], mo[5][1].real) + 1j * np.interp(w, mo[5][0], mo[5][1].imag)
    # V from the hand formula's fit to HAMS's own rotation term (checks the formula, gives V); A_wp z_o likewise
    k_rot = 0.5 * w**2 / G * np.real(xi5 * np.conj(xi1))
    k_hs = -0.25 * abs(xi5)**2
    V = np.sum(T[5] * k_rot) / np.sum(k_rot**2)
    AZ = np.sum(T[6] * k_hs) / np.sum(k_hs**2)
    rot_hand, hs_hand = V * k_rot, AZ * k_hs
    tsum = sum(T.values())
    diff = ref - ham
    m = np.isfinite(diff)
    X = np.c_[T[5][m], T[6][m]]
    coef, *_ = np.linalg.lstsq(X, diff[m], rcond=None)
    res = diff[m] - X @ coef
    i0 = np.argmin(abs(w - 0.314))
    lines = [
        "SOFTWIND heave mean drift, heading 0, /(rho g) per unit amplitude^2",
        "sum of the six terms vs QTF_Indirect.11d diagonal: max |difference| %.3g (the potential QTF part is 0 on the diagonal)"
        % np.nanmax(abs(tsum - ham)),
        "hand formulas: rotation fits with V = %.0f m3 (residual max %.3g), hydrostatic with A_wp z_o = %.1f m3 (residual max %.3g)"
        % (V, np.max(abs(rot_hand - T[5])), AZ, np.max(abs(hs_hand - T[6]))),
        "fit HydroStar - HAMS = a IV + b V over %d frequencies: a = %.3f, b = %.3f, residual max %.3g, rms %.3g (|diff| max %.3g)"
        % (m.sum(), coef[0], coef[1], np.max(abs(res)), np.sqrt(np.mean(res**2)), np.max(abs(diff[m]))),
        "at %.3f rad/s: HAMS %.2f, HydroStar %.2f, terms %s; |xi1| %.3f, |xi5| %.4f, phase(xi5) - phase(xi1) %.0f deg"
        % (w[i0], ham[i0], ref[i0], ", ".join("%s %.2f" % (NAMES[t], T[t][i0]) for t in range(1, 7)), abs(xi1[i0]),
           abs(xi5[i0]), np.degrees(np.angle(xi5[i0] / xi1[i0]))),
        "HAMS without term V (HydroStar's convention, Nemoh doc: HydroStar omits -[K] xi2, xi2_3 = z_G (theta_1^2 + theta_2^2)/2):"
        " residual against HydroStar max %.3g at %.3f rad/s, median %.3g"
        % (np.nanmax(abs(ref - (ham - T[6]))), w[np.nanargmax(abs(ref - (ham - T[6])))], np.nanmedian(abs(ref - (ham - T[6])))),
    ]
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    a = ax[0]
    a.plot(w, ref, "k--o", ms=3, label="HydroStar (Nemoh 3.0.3 test case 10)")
    a.plot(w, ham, color="tab:red", lw=2, label="HAMS-MREL, total")
    for t, c in ((3, "tab:blue"), (4, "tab:cyan"), (5, "tab:purple"), (6, "tab:orange")):
        a.plot(w, T[t], color=c, lw=1, label="HAMS-MREL term %s" % NAMES[t])
    a.plot(w, ham - T[6], color="tab:green", ls=":", lw=2.2,
           label="HAMS-MREL total without term V (HydroStar omits -[K] ξ̃⁽²⁾, Nemoh 3.0.3 doc)")
    a.set_xlim(0, 1.0)
    a.set_title("QTF⁻_3 diagonal (heave mean drift) / ρg, β = 0°, and HAMS-MREL's terms", fontsize=9)
    a.set_xlabel("ω [rad/s]")
    a.grid(alpha=0.3)
    a.legend(fontsize=7)
    a = ax[1]
    a.plot(w, diff, "k-o", ms=3, label="HydroStar − HAMS-MREL")
    a.plot(w[m], X @ coef, color="tab:green", lw=1.5, label="fit a·IV + b·V (a = %.2f, b = %.2f)" % tuple(coef))
    a.plot(w, T[5], color="tab:purple", lw=1, ls="--", label="HAMS-MREL term IV rotation")
    a.plot(w, -T[6], color="tab:orange", lw=1.5, ls="--", label="minus HAMS-MREL term V hydrostatic (the term HydroStar omits)")
    a.set_xlim(0, 1.0)
    a.set_title("HydroStar − HAMS-MREL against the rotation and hydrostatic terms", fontsize=9)
    a.set_xlabel("ω [rad/s]")
    a.grid(alpha=0.3)
    a.legend(fontsize=7)
    fig.suptitle("SOFTWIND, heave mean drift near the pitch resonance (0.314 rad/s): HAMS-MREL against HydroStar", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "sw_heave_drift_terms.png"), dpi=100)
    with open(os.path.join(out, "sw_heave_drift_terms.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
