"""Cost figure of the runs that would extend HAMS to the paper's full frequency range (snellius/COST.md, "Extension to
the full range"). Usage: python -I plot_cost_extension.py OUT_PNG"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
parts = ["First order", "Field passes / body part", "Free-surface part", "Direct method"]
colours = ["#9ec5f4", "#2a78d6", "#1c5cab", "#eb6834"]
# (expected hours per part, upper-bound hours per part, limit in hours)
options = {
    "OC4, 9208 panels\n0.25-2.55 rad/s (rec.)": ([3.5, 4.9, 0, 0], [3.5, 13.5, 0, 0], 22),
    "OC4, 5022 panels\n0.25-2.55 rad/s": ([1.0, 2.1, 0, 0], [1.0, 5.9, 0, 0], 9),
    "VolturnUS-S A, 0.25-2.50\ninner FS 1.25 m (rec.)": ([0.9, 1.6, 5.2, 2.2], [0.9, 5.8, 13.0, 4.0], 30),
    "VolturnUS-S A, 0.25-2.50\ninner FS 2 m": ([0.9, 1.6, 2.0, 2.2], [0.9, 5.8, 5.0, 4.0], 18),
    "VolturnUS-S A as approved\n0.25-1.90": ([0.8, 1.2, 1.5, 1.6], [0.8, 4.3, 7.7, 3.0], 18),
}
fig, ax = plt.subplots(figsize=(10.5, 6.6), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ticks = []
for k, (name, (exp, upper, limit)) in enumerate(options.items()):
    y = -1.0 * k
    for row, values in ((0, exp), (1, upper)):
        left = 0.0
        for v, c in zip(values, colours):
            if v:
                ax.barh(y - 0.24 * row, v * 24, left=left, height=0.22, color=c, edgecolor=SURFACE, linewidth=1.5)
            left += v * 24
        ax.text(left + 6, y - 0.24 * row, "%s %.0f" % ("expected" if row == 0 else "upper", left), va="center",
                fontsize=8, color=INK)
    ax.plot([limit * 24] * 2, [y - 0.38, y + 0.14], color=INK2, lw=1.2, ls=(0, (4, 2)))
    ax.text(limit * 24 + 6, y + 0.14, "limit %d h = %d (max billed)" % (limit, limit * 24), fontsize=7.5, color=INK2,
            va="center")
    ticks.append((y - 0.12, name))
ax.set_yticks([t for t, _ in ticks])
ax.set_yticklabels([n for _, n in ticks], fontsize=8.5)
ax.set_xlabel("core-hours (24 CPUs, fat_genoa, tdsei3290); one job each, none submitted")
ax.set_xlim(0, 900)
ax.grid(axis="x", color=GRID)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in colours], parts, loc="upper center",
          bbox_to_anchor=(0.5, -0.1), frameon=False, fontsize=8, ncol=4)
ax.set_title("Extending HAMS to the frequency range of Roland et al. (2026): full reruns with irregular-frequency removal",
             loc="left", fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=150)
