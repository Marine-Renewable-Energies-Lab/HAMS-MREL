"""OC4: the change of QTF_Indirect.11s and of its potential terms B2 and FS_inner against the old run, per sum frequency.

Usage: python -I oc4_by_omega.py NEW_RUN OLD_RUN OUT_DIR

Measure as qtf_compare.py: |new - old| / max |old| over the force group (modes 1-3) or the moment group (modes 4-6);
for the terms, over all modes of the term. Writes oc4_by_omega.txt and oc4_by_omega.png. The old run warns from
Omega = 1.60 rad/s (nu h > 50, Limes in the finite-depth series), so the table splits at Omega = 1.60.
"""
import sys
from pathlib import Path

import numpy as np

new, old, out = (Path(p) for p in sys.argv[1:4])
wf = Path("output", "Wamit_format")
a = np.loadtxt(old / wf / "QTF_Indirect.11s")
b = np.loadtxt(new / wf / "QTF_Indirect.11s")
ta = np.loadtxt(old / wf / "QTF_PotentialTermsSum.txt")
tb = np.loadtxt(new / wf / "QTF_PotentialTermsSum.txt")
if a.shape != b.shape or ta.shape != tb.shape or not (a[:, :5] == b[:, :5]).all():
    sys.exit("the row keys of the two runs differ")
W = np.round(a[:, 0] + a[:, 1], 3)
Ws = np.unique(W)
za, zb, m = a[:, 7] + 1j * a[:, 8], b[:, 7] + 1j * b[:, 8], a[:, 4]
curves, lines = {}, ["|new - old| / max |old| in %, max per sum frequency Omega (rad/s)"]
for mm in range(1, 7):
    g = (m <= 3) if mm <= 3 else (m >= 4)
    d = 100 * np.abs(za - zb) / np.abs(za[g]).max()
    curves[f".11s mode {mm}"] = [d[(m == mm) & (W == x)].max() for x in Ws]
TW = np.round(ta[:, 0] + ta[:, 1], 3)
for t, name in ((4, "B2"), (6, "FS_inner"), (7, "FS_outer")):
    c = 5 + 2 * (t - 1)
    xa, xb = ta[:, c] + 1j * ta[:, c + 1], tb[:, c] + 1j * tb[:, c + 1]
    d = 100 * np.abs(xa - xb) / np.abs(xa).max()
    curves[f"term {name}"] = [d[TW == x].max() for x in Ws]
lines.append("%-16s" % "Omega" + "".join("%7.2f" % x for x in Ws))
for k, v in curves.items():
    lines.append("%-16s" % k + "".join("%7.2f" % x for x in v))
lo = Ws < 1.575
lines.append("")
for k, v in curves.items():
    v = np.array(v)
    lines.append("%-16s max %6.2f %% for Omega <= 1.55, %6.2f %% for Omega >= 1.60" % (k, v[lo].max(), v[~lo].max()))
txt = "\n".join(lines)
(out / "oc4_by_omega.txt").write_text(txt + "\n")
print(txt)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
cols = ["#2a6fdb", "#e0892b", "#1b9e77", "#7a5bc7", "#c43c5c", "#7f7f7f"]
for i, mm in enumerate(range(1, 7)):
    axs[0].plot(Ws, curves[f".11s mode {mm}"], "-o", ms=3, lw=1.5, color=cols[i], label=f"QTF⁺_{mm}")
for i, n in enumerate(("B2", "FS_inner", "FS_outer")):
    axs[1].plot(Ws, curves[f"term {n}"], "-o", ms=3, lw=1.5, color=cols[i], label=n)
for ax, t in zip(axs, ("QTF_Indirect.11s, per mode (force or moment group)", "potential terms (all modes)")):
    ax.axvline(1.575, color="#999", lw=1, ls="--")
    ax.text(1.59, ax.get_ylim()[1] * 0.92, "old run warns:\nnu h > 50", fontsize=8, color="#555", va="top")
    ax.set_title(t, fontsize=10)
    ax.set_xlabel("sum frequency Omega (rad/s)")
    ax.grid(color="#ddd", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
axs[0].set_ylabel("max |new - old| / max |old|  (%)")
fig.suptitle("OC4 WAMIT model, QTF_Indirect.11s: HAMS with the speed-ups (deep-water rows on FS, C and F) "
             "against HAMS at 22086da (reference)", fontsize=10)
fig.tight_layout()
fig.savefig(out / "oc4_by_omega.png", dpi=150)
