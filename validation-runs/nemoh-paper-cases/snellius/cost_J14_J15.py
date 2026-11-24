"""Cost of J14 (Nemoh v3.0.3 on OC3, 20 frequencies) and J15 (HAMS-MREL OC3 bottom-edge meshes e0-e2, 20 frequencies),
fat_genoa, 24 CPUs billed per job (the minimum), 1.5 SBU per core-hour. Writes COST_J14_J15.png and prints the table.

HAMS-MREL (J15), from the stage times of J5C (same binary, 50 frequencies, 4710 hull + 96 lid panels, 24 threads, 32 127 s):
the B2 control-surface field passes of the first 20 frequencies (0.0628-1.256 rad/s) took 16 640 s, the rest of the run
6 328 s for 50 frequencies (127 s per frequency). The field pass costs one hull influence row per point of C and F, so it
is scaled by the number of hull panels (C is held at J5C's); the rest (solves, QTF sums) by (panels with lid)^2.5.

Nemoh (J14): one thread. Measured on Snellius (J8, J9; deep water, LU, 11 frequencies): 3495 panels 6.3 min, 4713 panels
13.3 min, about 73 s per frequency at 4700 panels. OC3 is in finite depth (h = 320 m): a factor 3 to 6 on the influence
matrices is assumed (not measured), and QTFpreProc is taken as one more pass of the same cost; QTFsolver for 20 x 20 pairs
as minutes. Range 2.6-5.2 h, expected 4 h.

Usage: python -I cost_J14_J15.py OUT_DIR
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CORES, SBU_PER_CH = 24, 1.5
FIELD20, REST_PER_F, NF = 16640.0, 6328.0 / 50, 20
J5C_HULL, J5C_ALL = 4710, 4806


def hams(name, hull, lid, limit):
    s = FIELD20 * hull / J5C_HULL + REST_PER_F * NF * ((hull + lid) / J5C_ALL) ** 2.5
    return dict(name=name, h=s / 3600, lim=limit)


nemoh_per_f = 73.0 * (4842 / 4713) ** 2
lo, hi = 2 * NF * nemoh_per_f * 3 / 3600, 2 * NF * nemoh_per_f * 6 / 3600
jobs = [dict(name="J14 Nemoh 3.0.3, OC3 test case (4842 panels), 1 thread", h=4.0, lim=12, rng=(lo, hi)),
        hams("J15 oc3_e0, uniform (5138 hull panels)", 5138, 196, 9),
        hams("J15 oc3_e1, edge row 0.375 m (5306)", 5306, 196, 10),
        hams("J15 oc3_e2, edge row 0.094 m (5866)", 5866, 196, 11)]
for j in jobs:
    j["sbu"] = j["h"] * CORES * SBU_PER_CH
    j["sbul"] = j["lim"] * CORES * SBU_PER_CH
tot = dict(h=sum(j["h"] for j in jobs), sbu=sum(j["sbu"] for j in jobs), sbul=sum(j["sbul"] for j in jobs))
lines = ["| Job | Expected wall (h) | Core-h | SBU | Limit (h) | SBU at limit |", "|---|---:|---:|---:|---:|---:|"]
for j in jobs:
    extra = " (range %.1f-%.1f)" % j["rng"] if "rng" in j else ""
    lines.append("| %s | %.1f%s | %.0f | %.0f | %d | %.0f |" % (j["name"], j["h"], extra, j["h"] * CORES, j["sbu"], j["lim"], j["sbul"]))
lines.append("| **All four** | | %.0f | **%.0f** | | **%.0f** |" % (tot["h"] * CORES, tot["sbu"], tot["sbul"]))
sys.stdout.buffer.write(("\n".join(lines) + "\n").encode("utf-8"))

fig, ax = plt.subplots(figsize=(10.5, 4.4))
y = list(range(len(jobs)))
ax.barh([v + 0.18 for v in y], [j["sbu"] for j in jobs], height=0.35, color="tab:blue", label="expected")
ax.barh([v - 0.18 for v in y], [j["sbul"] for j in jobs], height=0.35, color="tab:gray", alpha=0.6, label="at the time limit")
for v, j in zip(y, jobs):
    ax.text(j["sbu"] + 5, v + 0.18, "%.0f SBU, %.1f h" % (j["sbu"], j["h"]), va="center", fontsize=8)
    ax.text(j["sbul"] + 5, v - 0.18, "%.0f SBU, %d h limit" % (j["sbul"], j["lim"]), va="center", fontsize=8)
ax.set_yticks(y)
ax.set_yticklabels([j["name"] for j in jobs], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("SBU (fat_genoa, 24 CPUs billed, 1.5 SBU per core-hour)")
ax.set_xlim(0, max(j["sbul"] for j in jobs) * 1.35)
ax.legend(fontsize=8, loc="lower right")
ax.set_title("OC3 low-row checks, not submitted: all four %.0f SBU expected, %.0f SBU at the limits" % (tot["sbu"], tot["sbul"]),
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(sys.argv[1], "COST_J14_J15.png"), dpi=100)
