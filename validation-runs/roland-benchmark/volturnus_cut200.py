"""VolturnUS-S, difference frequency, heading 0: HAMS-MREL with the free-surface integral cut at 200 m against WAMIT.

HAMS cut at 200 m = QTF_Indirect.11d of q_volturnus_free_full_ext (inner free surface to 72 m, partial waves beyond) minus
FS_outer of q_volturnus_r200 (the integral beyond R = 200 m; with qtf_fs_circle = panels its partial waves come from the hull
integral, so they do not depend on the inner mesh). The body terms of the two runs must agree; the script checks it.
Pairs hit by the FS_outer defect (non-finite, or a tail above 10 times the larger of |.12d| and |WAMIT|) are left out and
counted. Along lines of constant omega_1 - omega_2 it fits X of exp(-i (k1 - k2) X), finite depth h = 200 m.

WAMIT: v6.107S, IEA-15-240-RWT-UMaineSemi.12d (IEA Wind Task 37, IEA-15-240-RWT v1.1.18, Apache-2.0), RINNER = 200 m,
NAL = 0. Not in the repository.
Usage: python -I volturnus_cut200.py WAMIT_12D FULL_WAMIT_FORMAT_DIR R200_WAMIT_FORMAT_DIR OUT_PNG OUT_TXT
"""
import cmath
import math
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RHO_G = 1025.0 * 9.80665


def fnum(x):
    """Fortran drops the E of a three-digit exponent (-1.18+103): an overflowed value, read as inf."""
    try:
        return float(x)
    except ValueError:
        return float("inf")


def qtf(path, periods=False, modes=(1, 3, 5)):
    """{(mode, omega_1, omega_2): QTF in N/m^2 (N m/m^2 for moments)}, heading 0, both orders of the pair."""
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) != 9:
            continue
        v = [fnum(x) for x in p]
        m = int(v[4])
        if v[2] or v[3] or m not in modes:
            continue
        a, b = (2 * math.pi / v[0], 2 * math.pi / v[1]) if periods else (v[0], v[1])
        a, b = round(a, 2), round(b, 2)
        c = complex(v[7], v[8]) * RHO_G
        out[(m, a, b)] = c
        out.setdefault((m, b, a), c.conjugate())
    return out


def terms(path, n):
    """The n complex terms per row of QTF_PotentialTerms.txt (or QTF_DirectTerms.txt), heading 0, in N/m^2."""
    out = {}
    for line in open(path):
        if line.lstrip().startswith("#"):
            continue
        p = line.split()
        if len(p) != 5 + 2 * n or fnum(p[2]) or fnum(p[3]):
            continue
        m, a, b = int(p[4]), round(float(p[0]), 2), round(float(p[1]), 2)
        t = [complex(fnum(p[5 + 2 * k]), fnum(p[6 + 2 * k])) * RHO_G for k in range(n)]
        out[(m, a, b)] = t
        out.setdefault((m, b, a), [x.conjugate() for x in t])
    return out


def kk(w, h=200.0, g=9.80665):
    k = w * w / g
    for _ in range(60):
        k -= (k * math.tanh(k * h) - w * w / g) / (math.tanh(k * h) + k * h / math.cosh(k * h) ** 2)
    return k


def fit_x(get, m, dw, top=2.25):
    """Length X of the phase law exp(-i (k1 - k2) X) along the line omega_1 - omega_2 = dw, and its coherence."""
    num, den, dk = 0j, 0.0, []
    for i in range(6, 51):
        w2 = round(i * 0.05, 2)
        w1 = round(w2 + dw, 2)
        a, b = (m, w1, w2), (m, round(w1 + 0.05, 2), round(w2 + 0.05, 2))
        if b[1] > top + 1e-9:
            continue
        try:
            qa, qb = get(a), get(b)
        except KeyError:
            continue
        z = qb * qa.conjugate()
        num += z
        den += abs(z)
        dk.append(kk(b[1]) - kk(b[2]) - kk(a[1]) + kk(a[2]))
    if not dk:
        return float("nan"), 0.0
    return -cmath.phase(num) / (sum(dk) / len(dk)), abs(num) / den


INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
COL = {1: "#2a78d6", 3: "#eb6834", 5: "#1baf7a"}
PINK = "#e87ba4"
med = statistics.median


def main(wamit, full, r200, out_png, out_txt):
    W = qtf(wamit, periods=True)
    H12, H11 = qtf(os.path.join(full, "QTF_Direct.12d")), qtf(os.path.join(full, "QTF_Indirect.11d"))
    P = terms(os.path.join(full, "QTF_PotentialTerms.txt"), 7)
    N = terms(os.path.join(r200, "QTF_PotentialTerms.txt"), 7)
    lines = []
    # 1. the body terms (FK_I2, incident diffraction, B1, B2, B3) of the two runs must agree
    rel = []
    for k, t in N.items():
        if k in P and k[0] in (1, 3, 5):
            a, b = sum(t[:5]), sum(P[k][:5])
            if abs(b) > 0:
                rel.append(abs(a - b) / abs(b))
    lines.append("Body terms (1-5) of q_volturnus_r200 against q_volturnus_free_full_ext: median %.2e, max %.2e over %d rows"
                 % (med(rel), max(rel), len(rel)))
    # 2. the cut and the defect pairs
    T, bad = {}, {1: [], 3: [], 5: []}
    for k in H11:
        if k[0] not in (1, 3, 5) or k not in N or k not in W or k not in P or not k[1] > k[2]:
            continue
        o200, o72 = N[k][6], P[k][6]
        ok = all(math.isfinite(v) for v in (o200.real, o200.imag, o72.real, o72.imag))
        # a defect pair: non-finite, or a tail larger than 10 times the direct HAMS QTF (.12d, no partial-wave sum) or WAMIT of that pair
        if not ok or abs(o200) > 10 * max(abs(H12[k]), abs(W[k])) or abs(o72) > 10 * max(abs(H12[k]), abs(W[k])):
            bad[k[0]].append(k)
            continue
        T[k] = H11[k] - o200
    for m in (1, 3, 5):
        inband = [k for k in bad[m] if 0.3 - 1e-9 < k[1] - k[2] <= 1.2 + 1e-9 and k[2] >= 0.30 - 1e-9]
        lines.append("QTF-_%d: %d defect pairs left out (%d of them with 0.3 < dw <= 1.2 and both omega >= 0.30); lowest/highest "
                     "omega of them: %s" % (m, len(bad[m]), len(inband),
                                            "%.2f / %.2f" % (min(k[2] for k in bad[m]), max(k[1] for k in bad[m]))
                                            if bad[m] else "-"))

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.facecolor": SURFACE, "figure.facecolor": SURFACE, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.6), layout="constrained")
    edges = [0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.5]

    def sel(m, lo, hi):
        return [k for k in T if k[0] == m and k[2] >= 0.30 - 1e-9 and k[1] <= 2.50 + 1e-9 and lo + 1e-9 < k[1] - k[2] <= hi + 1e-9]

    # (a) median |WAMIT| / |HAMS| per band, full and cut at 200 m
    a = ax[0, 0]
    lines.append("")
    lines.append("Per band of omega_1 - omega_2 (both omega 0.30-2.50): median |W|/|H full|, median |W|/|H cut 200|, "
                 "median |W - H cut| / |W - H full|, median |FS beyond 200| / |W - H full|, pairs")
    for m in (1, 3, 5):
        xs, rf, rc = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            ks = sel(m, lo, hi)
            if not ks:
                continue
            xs.append(0.5 * (lo + hi))
            rf.append(med(abs(W[k]) / abs(H11[k]) for k in ks))
            rc.append(med(abs(W[k]) / abs(T[k]) for k in ks))
            lines.append("  QTF-_%d %4.2f-%4.2f  %6.2f  %6.2f  %6.2f  %6.3f  %d" % (
                m, lo, hi, rf[-1], rc[-1], med(abs(W[k] - T[k]) / abs(W[k] - H11[k]) for k in ks),
                med(abs(N[k][6]) / abs(W[k] - H11[k]) for k in ks), len(ks)))
        a.plot(xs, rf, color=COL[m], lw=1.4, ls=(0, (5, 2)), marker="s", ms=3,
               label=r"$QTF^{-}_{%d}$: WAMIT / HAMS .11d (free surface to infinity)" % m)
        a.plot(xs, rc, color=COL[m], lw=2.0, marker="o", ms=3.5,
               label=r"$QTF^{-}_{%d}$: WAMIT / HAMS .11d cut at 200 m" % m)
    a.set_yscale("log")
    a.axhline(1, color=INK2, lw=0.8)
    a.set_xlabel(r"$\omega_1 - \omega_2$ (rad/s), band centre")
    a.set_ylabel("median |WAMIT| / |HAMS-MREL| over the band")
    a.set_title("(a) Does cutting HAMS at 200 m close the gap?", loc="left", color=INK)
    a.legend(fontsize=6.6, loc="upper left")

    # (b) QTF-_3 along dw = 0.70
    b = ax[0, 1]
    m, dw = 3, 0.70
    ks = [(m, round(i * 0.05 + dw, 2), round(i * 0.05, 2)) for i in range(6, 37)]
    ks = [k for k in ks if k[1] <= 2.5 + 1e-9 and k in T]
    w1 = [k[1] for k in ks]
    b.plot(w1, [abs(W[k]) for k in ks], color=INK, lw=2.0, marker="o", ms=3.5,
           label="WAMIT v6.107S, IEA-15-240-RWT-UMaineSemi.12d (IEA Wind Task 37)")
    b.plot(w1, [abs(H11[k]) for k in ks], color=COL[1], lw=1.4, ls=(0, (5, 2)), marker="s", ms=3,
           label="HAMS-MREL QTF_Indirect.11d, free surface to infinity")
    b.plot(w1, [abs(T[k]) for k in ks], color=COL[3], lw=2.0, marker="o", ms=3.5, label="HAMS-MREL .11d cut at 200 m")
    b.plot(w1, [abs(N[k][6]) for k in ks], color=PINK, lw=1.2, ls=":", label="HAMS-MREL free surface beyond 200 m alone")
    b.set_yscale("log")
    b.set_xlabel(r"$\omega_1$ (rad/s), with $\omega_2 = \omega_1 - 0.70$ rad/s")
    b.set_ylabel(r"$|QTF^{-}_{3}|$ (N/m$^2$)")
    b.set_title(r"(b) $|QTF^{-}_{3}|$ along $\omega_1-\omega_2$ = 0.70 rad/s", loc="left", color=INK)
    b.legend(fontsize=7, loc="lower left")

    # (c) phase-law X of the tail beyond 200 m against that of WAMIT - HAMS
    c = ax[1, 0]
    dws = [round(0.1 * i, 2) for i in range(1, 13)]
    lines.append("")
    lines.append("X (m) of exp(-i (k1-k2) X) along constant dw, dw:X(coherence)")

    def safe(f):
        def g(k):
            if k not in T:
                raise KeyError(k)
            return f(k)
        return g
    for m in (1, 3, 5):
        xo = [fit_x(safe(lambda k: N[k][6]), m, d, 2.5) for d in dws]
        xr = [fit_x(safe(lambda k: W[k] - H11[k]), m, d, 2.5) for d in dws]
        c.plot(dws, [v[0] for v in xo], color=COL[m], lw=1.4, ls=(0, (5, 2)), marker="s", ms=3,
               label=r"HAMS free surface beyond 200 m, $QTF^{-}_{%d}$" % m)
        c.plot(dws, [v[0] for v in xr], color=COL[m], lw=2.0, marker="o", ms=3.5,
               label=r"WAMIT $-$ HAMS .11d, $QTF^{-}_{%d}$" % m)
        for name, xv in (("beyond 200 m", xo), ("WAMIT - HAMS .11d", xr)):
            lines.append("  QTF-_%d %-18s %s" % (m, name, " ".join("%.1f:%.0f(%.2f)" % (d, v[0], v[1]) for d, v in zip(dws, xv))))
    c.axhline(200, color=INK2, lw=0.9, ls="--")
    c.text(1.2, 232, "200 m (WAMIT RINNER; HAMS circle about (-12.9, 0))", fontsize=7.5, color=INK2, ha="right", va="bottom")
    c.set_ylim(-60, 300)
    c.set_xlabel(r"$\omega_1 - \omega_2$ (rad/s)")
    c.set_ylabel("X (m)")
    c.set_title(r"(c) X fitted to $e^{-i(k_1-k_2)X}$ along lines of constant $\omega_1-\omega_2$", loc="left", color=INK)
    c.legend(fontsize=6.4, ncol=2, loc="lower center")

    # (d) sizes per band
    d = ax[1, 1]
    for m in (1, 3):
        xs, y200, y72, yex = [], [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            ks = sel(m, lo, hi)
            if not ks:
                continue
            xs.append(0.5 * (lo + hi))
            y200.append(med(abs(N[k][6]) for k in ks))
            y72.append(med(abs(P[k][6]) for k in ks))
            yex.append(med(abs(W[k] - H11[k]) for k in ks))
        d.plot(xs, yex, color=COL[m], lw=2.0, marker="o", ms=3.5, label=r"$|$WAMIT $-$ HAMS .11d$|$, $QTF^{-}_{%d}$" % m)
        d.plot(xs, y72, color=COL[m], lw=1.2, ls=":", marker="^", ms=3,
               label=r"HAMS free surface beyond 72 m, $QTF^{-}_{%d}$" % m)
        d.plot(xs, y200, color=COL[m], lw=1.4, ls=(0, (5, 2)), marker="s", ms=3,
               label=r"HAMS free surface beyond 200 m, $QTF^{-}_{%d}$" % m)
    d.set_yscale("log")
    d.set_xlabel(r"$\omega_1 - \omega_2$ (rad/s), band centre")
    d.set_ylabel(r"median modulus over the band (N/m$^2$)")
    d.set_title("(d) The tail beyond 200 m against the WAMIT excess", loc="left", color=INK)
    d.legend(fontsize=6.6, loc="lower right")

    fig.suptitle("VolturnUS-S, difference-frequency QTF, heading 0: HAMS-MREL with the free-surface integral cut at 200 m "
                 "against WAMIT (defect pairs left out)", x=0.01, ha="left", color=INK)
    fig.savefig(out_png, dpi=150)
    open(out_txt, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
