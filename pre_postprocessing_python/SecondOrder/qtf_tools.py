"""
Read the quadratic QTFs of HAMS-MREL, apply Newman's approximation, and compute the spectrum of the slowly varying force.

Usage:
    python qtf_tools.py check
    python qtf_tools.py newman <QTF_Quadratic.10d> --index I [--beta DEG] [--out CSV]

Files. QTF_Quadratic.10d and .10s (written when qtf = 1, 2 or 3 in SecondOrder.in) hold one row per pair:

    OUFR1 OUFR2 BETA1 BETA2 I Mod Pha Re Im

OUFR is the frequency in the form of Output_frequency_type in ControlFile.in, I = 6(b - 1) + 1 ... 6 for body b, and
the values are scaled by rho g L**k (k = 1 for forces, 2 for moments), in the e^{+i omega t} convention of the WAMIT-format
files: the conjugates of the internal T. With the internal time factor e^{-i omega t} and components a_i of the wave
amplitude, the quadratic second-order load is (docs/output_files.md)

    F(t) = sum_i sum_j Re[ T-_ij a_i conj(a_j) exp(-i (w_i - w_j) t) + T+_ij a_i a_j exp(-i (w_i + w_j) t) ],

over all i and j, so T-_ij = conj(T-_ji), T+_ij = T+_ji, and T-_ii is the mean drift per unit amplitude squared. The rows
of the file give conj(T): read_qtf returns the file values; internal() conjugates them back.

Newman's approximation (Newman 1974) replaces the off-diagonal difference QTF by the mean of the diagonal values,
T-_ij ~ (T-_ii + T-_jj)/2. It needs only the mean drift, and holds when the QTF varies slowly away from the diagonal.

Slow-drift spectrum. For a random sea with one-sided wave spectrum S(w), a_i = sqrt(2 S(w_i) dw) exp(i eps_i) with
independent uniform phases. The slowly varying part (i /= j) has the one-sided spectrum, at mu = k dw,

    S_F(mu) = 8 sum_i S(w_i) S(w_i + mu) |T-(w_i, w_i + mu)|**2 dw,

so that its variance is sum_k S_F(mu_k) dw [derived; the continuous form 8 int S(w) S(w + mu) |T-|**2 dw is the one of
Pinkster 1980]. check_spectrum_variance tests the constant 8 against a Monte Carlo time series of F(t).

Sources:
    J. N. Newman, Second-order, slowly-varying forces on vessels in irregular waves, Proc. Int. Symp. Dynamics of Marine
        Vehicles and Structures in Waves, London, 1974, pp. 182-186.
    J. A. Pinkster, Low frequency second order wave exciting forces on floating structures, PhD thesis, Delft University
        of Technology, 1980.
"""
import argparse
import csv
import math
import os
import sys
import tempfile

import numpy as np


def read_qtf(path):
    """{(oufr1, oufr2, beta1, beta2, I): complex} of the file values (e^{+i omega t}, scaled), keys rounded to 1e-6."""
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) != 9:
                continue
            key = (round(float(p[0]), 6), round(float(p[1]), 6), round(float(p[2]), 6), round(float(p[3]), 6), int(p[4]))
            out[key] = complex(float(p[7]), float(p[8]))
    return out


def internal(qtf):
    """The same dict in the internal e^{-i omega t} convention: the conjugates."""
    return {key: value.conjugate() for key, value in qtf.items()}


def frequencies(qtf):
    """Sorted first-column frequencies of a QTF dict."""
    return sorted({key[0] for key in qtf})


def matrix(qtf, index, beta1=0.0, beta2=None):
    """(omegas, M) with M[i, j] the value at (omega_i, omega_j, beta1, beta2, index), NaN where the file has no row."""
    beta2 = beta1 if beta2 is None else beta2
    omegas = frequencies(qtf)
    m = np.full((len(omegas), len(omegas)), np.nan + 1j * np.nan, dtype=complex)
    for i, wi in enumerate(omegas):
        for j, wj in enumerate(omegas):
            value = qtf.get((wi, wj, round(beta1, 6), round(beta2, 6), index))
            if value is not None:
                m[i, j] = value
    return np.array(omegas), m


def newman(diagonal):
    """Newman's approximation from the diagonal values d (n,): M[i, j] = (d_i + d_j)/2, real."""
    d = np.real(np.asarray(diagonal))
    return 0.5 * (d[:, None] + d[None, :])


def slow_drift_spectrum(spectrum, qtf_minus, dw):
    """(mu, S_F): mu_k = k dw for k = 1 ... n - 1 and S_F(mu_k) = 8 sum_i S_i S_{i+k} |T-_{i,i+k}|**2 dw, on a uniform grid
    with spacing dw; spectrum (n,) is S(w_i), qtf_minus (n, n) the difference QTF (either convention: only |T| is used)."""
    s = np.asarray(spectrum, dtype=float)
    t = np.abs(np.asarray(qtf_minus)) ** 2
    n = s.size
    mu = dw * np.arange(1, n)
    sf = np.array([8.0 * np.sum(s[:n - k] * s[k:] * np.diagonal(t, offset=k)) * dw for k in range(1, n)])
    return mu, sf


# ------------------------------------------------------------------------------------------------------------- checks

def check_newman_diagonal():
    """Newman's approximation is exact on the diagonal and symmetric."""
    d = np.array([0.1, 0.4, 0.9, 1.3])
    m = newman(d)
    assert np.allclose(np.diagonal(m), d) and np.allclose(m, m.T)
    return "max diagonal error %.1e" % np.max(np.abs(np.diagonal(m) - d))


def check_spectrum_variance(realisations=200, seed=7):
    """The variance of the slowly varying part of F(t), from random phases and the double sum with a complex QTF that
    varies over the grid, against sum_k S_F(mu_k) dw. Agreement within 3% over the ensemble."""
    rng = np.random.default_rng(seed)
    n, w0, dw = 24, 0.4, 0.04
    w = w0 + dw * np.arange(n)
    s = np.exp(-((w - 0.8) / 0.2) ** 2)                               # a smooth one-sided spectrum, m^2 s
    ti, tj = np.meshgrid(w, w, indexing="ij")
    t = (1.0 + 0.3 * ti + 0.2 * tj) + 1j * 0.4 * (ti - tj)              # T-_ij, Hermitian
    t = 0.5 * (t + t.conj().T)
    _, sf = slow_drift_spectrum(s, t, dw)
    expected = np.sum(sf) * dw
    period = 2.0 * math.pi / dw
    times = np.linspace(0.0, period, 256, endpoint=False)
    phase = np.exp(-1j * np.outer(times, w))                           # (nt, n): e^{-i w_i t}
    off = ~np.eye(n, dtype=bool)
    total = 0.0
    for _ in range(realisations):
        a = np.sqrt(2.0 * s * dw) * np.exp(1j * rng.uniform(0.0, 2.0 * math.pi, n))
        b = a[None, :] * phase                                          # a_i e^{-i w_i t}
        f = np.real(np.einsum("ti,ij,tj->t", b, np.where(off, t, 0.0), b.conj()))
        total += np.mean(f ** 2)
    measured = total / realisations
    rel = abs(measured - expected) / expected
    assert rel < 0.03, "variance %.4g against %.4g" % (measured, expected)
    return "variance %.5g, expected %.5g, relative %.2e" % (measured, expected, rel)


def check_read_round_trip():
    """read_qtf reads the format of write_qtf_row back, and internal() conjugates."""
    value = complex(0.25, -1.5)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "QTF_Quadratic.10d")
        with open(path, "w") as f:
            mod, pha = abs(value), math.degrees(math.atan2(value.imag, value.real))
            f.write("%14.6E%14.6E%14.6E%14.6E%6d%14.6E%14.6E%14.6E%14.6E\n" % (0.7, 0.95, 0.0, 30.0, 5, mod, pha,
                                                                              value.real, value.imag))
        q = read_qtf(path)
    assert q[(0.7, 0.95, 0.0, 30.0, 5)] == value and internal(q)[(0.7, 0.95, 0.0, 30.0, 5)] == value.conjugate()
    return "ok"


CHECKS = {"newman_diagonal": check_newman_diagonal, "spectrum_variance": check_spectrum_variance,
          "read_round_trip": check_read_round_trip}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check")
    p_new = sub.add_parser("newman")
    p_new.add_argument("path")
    p_new.add_argument("--index", type=int, required=True)
    p_new.add_argument("--beta", type=float, default=0.0)
    p_new.add_argument("--out")
    args = parser.parse_args(argv)
    if args.command == "check":
        failed = 0
        for name, func in CHECKS.items():
            try:
                print("%-20s %s" % (name, func()))
            except AssertionError as e:
                failed += 1
                print("%-20s FAILED %s" % (name, e))
        return 1 if failed else 0
    omegas, m = matrix(read_qtf(args.path), args.index, args.beta)
    approx = newman(np.diagonal(m))
    rows = [(wi, wj, m[i, j].real, m[i, j].imag, approx[i, j]) for i, wi in enumerate(omegas) for j, wj in enumerate(omegas)]
    out = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.writer(out)
    writer.writerow(["omega_i", "omega_j", "re", "im", "newman"])
    writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
