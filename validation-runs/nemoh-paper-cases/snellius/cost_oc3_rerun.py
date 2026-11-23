"""Cost of an OC3 difference-QTF rerun on 62e46e2 (direct field sum above nu*h = 50), from the stage times of J5B (the same
case on 1515e89, 24 threads on fat_genoa: 45 527 s wall, of which 39 210 s for the B2 control-surface fields at the 50
frequencies, 5 422 s of them at the four frequencies up to 0.251 rad/s). Writes COST_OC3_rerun.png and prints the table.
fat_genoa bills 24 CPUs at least, 1.5 SBU per core-hour. The direct sum costs 8.4 us per near pair against 26 us for
PIDCOCK, so the field times above nu*h = 50 are taken as in J5B (an upper estimate).

Usage: python -I cost_oc3_rerun.py OUT_DIR
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CORES, SBU_PER_CH, BUILD_S = 24, 1.5, 300
FIELDS, REST, LOW4, NF = 39210.0, 45527.0 - 39210.0, 5422.0, 50
per_other = (FIELDS - LOW4) / (NF - 4)


def job(name, wall_s, limit_h):
    h = (wall_s + BUILD_S) / 3600
    return dict(name=name, h=h, ch=h * CORES, sbu=h * CORES * SBU_PER_CH, lim=limit_h, chl=limit_h * CORES,
                sbul=limit_h * CORES * SBU_PER_CH)


# (ii) full rerun, same grid as J5B
full = job("(ii) full OC3 rerun, 50 frequencies (as J5B)", FIELDS + REST, 16)
# (i-a) rows omega_1 <= 0.25 on the same grid: every omega_2 still needs its fields, so the same cost
rows_same = job("(i-a) rows ω₁ ≤ 0.25, same 50-frequency grid", FIELDS + REST, 16)
# (i-b) rows omega_1 <= 0.25 with omega_2 on every third grid frequency from 0.314: 4 + 15 = 19 frequencies; the rest of
# the run (first-order and psi solves, pair sums) scaled by the frequency count, plus up to 15 extra psi solves (~60 s)
nf = 4 + 15
rows_coarse = job("(i-b) rows ω₁ ≤ 0.25, ω₂ every 3rd frequency (19 frequencies)",
                  LOW4 + 15 * per_other + REST * nf / NF + 15 * 60, 8)
jobs = [full, rows_same, rows_coarse]
for j in jobs:
    sys.stdout.buffer.write(("%-62s %5.1f h %6.0f core-h %6.0f SBU | limit %2d h: %4.0f core-h %4.0f SBU\n"
                             % (j["name"], j["h"], j["ch"], j["sbu"], j["lim"], j["chl"], j["sbul"])).encode("utf-8"))

fig, ax = plt.subplots(figsize=(10, 4.2))
y = range(len(jobs))
ax.barh([v + 0.18 for v in y], [j["sbu"] for j in jobs], height=0.35, color="tab:blue", label="expected (J5B stage times)")
ax.barh([v - 0.18 for v in y], [j["sbul"] for j in jobs], height=0.35, color="tab:gray", alpha=0.6,
        label="at the time limit")
for v, j in zip(y, jobs):
    ax.text(j["sbu"] + 5, v + 0.18, "%.0f SBU, %.1f h" % (j["sbu"], j["h"]), va="center", fontsize=8)
    ax.text(j["sbul"] + 5, v - 0.18, "%.0f SBU, %d h limit" % (j["sbul"], j["lim"]), va="center", fontsize=8)
ax.set_yticks(list(y))
ax.set_yticklabels([j["name"] for j in jobs], fontsize=8)
ax.set_xlabel("SBU (fat_genoa, 24 CPUs, 1.5 SBU per core-hour)")
ax.set_xlim(0, max(j["sbul"] for j in jobs) * 1.3)
ax.legend(fontsize=8, loc="lower right")
ax.set_title("OC3 difference-QTF rerun on 62e46e2: cost estimated from J5B (not submitted)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(sys.argv[1], "COST_OC3_rerun.png"), dpi=100)
