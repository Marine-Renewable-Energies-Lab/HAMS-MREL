"""Cost figure of the VolturnUS-S run (snellius/COST.md). Usage: python -I plot_cost.py OUT_PNG"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
parts = ["First order", "Body part", "Free-surface part", "Direct method"]
colours = ["#9ec5f4", "#2a78d6", "#1c5cab", "#eb6834"]
options = {"A: indirect + direct\n(recommended)": ([0.8, 1.2, 1.5, 1.6], [0.8, 4.3, 7.7, 3.0], 18),
           "B: indirect only": ([0.8, 1.2, 1.5, 0.0], [0.8, 4.3, 7.7, 0.0], 14)}
fig, ax = plt.subplots(figsize=(9.5, 4.2), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
y = 0
labels = []
for name, (exp, upper, limit) in options.items():
    for row, values in ((0, exp), (1, upper)):
        left = 0.0
        for v, c in zip(values, colours):
            if v:
                ax.barh(y - 0.2 * row, v * 24, left=left, height=0.18, color=c, edgecolor=SURFACE, linewidth=1.5)
            left += v * 24
        ax.text(left + 6, y - 0.2 * row, "%s %.0f core-h" % ("expected" if row == 0 else "upper bound", left),
                va="center", fontsize=8.5, color=INK)
    ax.plot([limit * 24] * 2, [y - 0.33, y + 0.13], color=INK2, lw=1.2, ls=(0, (4, 2)))
    ax.text(limit * 24 + 6, y + 0.12, "limit %d h = %d core-h (max billed)" % (limit, limit * 24), fontsize=8,
            color=INK2, va="center")
    labels.append((y - 0.1, name))
    y -= 1.0
ax.set_yticks([p for p, _ in labels])
ax.set_yticklabels([n for _, n in labels])
ax.set_xlabel("core-hours (24 CPUs, fat_genoa)")
ax.set_xlim(0, 620)
ax.grid(axis="x", color=GRID)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colours]
ax.legend(handles, parts, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=8, ncol=4)
ax.set_title("VolturnUS-S HAMS run with the free-surface integral: one job, no build (existing binary)", loc="left",
             fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=150)
