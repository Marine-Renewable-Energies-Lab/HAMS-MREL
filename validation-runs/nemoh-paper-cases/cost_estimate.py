"""Cost estimate of the Snellius jobs of the Nemoh-paper cases (fat_genoa, 24 CPUs billed, SBU = 1.5 x core-h).

Basis: the measured log of the OC4 run of 9208 hull panels and a 1674-panel lid (47 frequencies, difference QTF with
qtf_terms = body, 24 threads, 22378 s; validation-runs/roland-benchmark): one B2 field pass over C and F (4454 points)
took 536 s where the finite-depth Green function was used and 8.2 s with the deep-water rows; the first-order solve and
the rest took 295 s per frequency. Scaling [estimate]:
- B2 pass: 536 s x (N_hull / 9208) x (N_points / 4454), N_hull the whole hull (a half mesh loops over both images);
  one pass per frequency of the list and one per radiation-only Omega (none on these uniform grids, except the sums
  above the list). Deep-water rows (1/65 of that) only where field_deep_draft allows them (depth >= 10 drafts) and
  k0 h >= 10.
- First order: 295 s x (N_eff / 9208)^2.5 per frequency, N_eff = hull + lid, x 0.6 for a half mesh (two half-size
  solves).
- Free-surface integral (sw_full): the free-surface part of the OC4 sum run of the old code (finite depth everywhere,
  5022 panels, 4044 inner points, 21 frequencies: 3.83 h), scaled by N_hull, the inner points and the frequency count.
- Nemoh (one thread): its OC4 first-order stage, 2.33 h for 21 frequencies on 5022 panels, scaled by N^2.5.
The limit is 2.5 x the expected wall time, rounded up to the hour, at least 3 h; a failed start stops in minutes.

Usage: python -I cost_estimate.py CASES_JSON OUT_PNG OUT_MD
"""
import json
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PASS_REF = 536.0 / (9208 * 4454)
DEEP_FACTOR = 8.2 / 536.0
FO_REF = 295.0
CPUS = 24


def fo(n_eff, nfreq, half=False):
    return FO_REF * (n_eff * (0.6 if half else 1.0) / 9208.0) ** 2.5 * nfreq


def passes(n_hull, pts, freqs, depth, draft, extra_omegas=()):
    """B2 passes at every frequency of the list and at the extra (radiation-only) Omega."""
    t = 0.0
    deep_ok = depth < 0 or depth >= 10.0 * draft
    for w in list(freqs) + list(extra_omegas):
        k0h = w * w / 9.80665 * depth if depth > 0 else 1e9
        f = DEEP_FACTOR if (deep_ok and k0h >= 10.0) else 1.0
        t += PASS_REF * n_hull * pts * f
    return t


def grid(n, step=0.0628):
    return [step * (k + 1) for k in range(n)]


def estimate(c):
    rows = []
    # J1: small cases in one job (HAMS on 23 threads, Nemoh on 1 thread beside it)
    hp1 = fo(c["hp_m1"]["hull"] + c["hp_m1"]["lid"]["panels"], 20)
    hp2 = fo(c["hp_m2"]["hull"] + c["hp_m2"]["lid"]["panels"], 20)
    swfo = fo(2 * c["sw_fo"]["hull"] + 2 * c["sw_fo"]["lid"]["panels"], 100, half=True)
    cy = c["cyl8b_diff"]
    cyl = fo(cy["hull"] + cy["lid"]["panels"], 50) + passes(cy["hull"], cy["cs"]["points"], grid(50), -1.0, 14.0)
    ham = hp1 + hp2 + swfo + cyl
    nem = 2.33 * 3600 / 21 * (c["nemoh_hp_m1"]["panels"] / 5022.0) ** 2.5 * 20
    rows.append(("J1 small", "hp_m1, hp_m2 (heave plate), cyl8b_diff, sw_fo; Nemoh nemoh_hp_m1 beside",
                 max(ham, nem) + 600))
    # J2: SOFTWIND difference, no FS
    s = c["sw_diff"]
    nfull = 2 * s["hull"]
    t = fo(nfull + 2 * s["lid"]["panels"], 50, half=True) + passes(nfull, s["cs"]["points"], grid(50), 200.0, 91.4)
    rows.append(("J2 sw_diff", "SOFTWIND QTF-, 50 freq. to 3.14 rad/s, no FS", t + 600))
    # J3: SOFTWIND difference + sum with FS (32 frequencies to 2.01 rad/s)
    s = c["sw_full"]
    sums_off = [0.0628 * k for k in range(33, 65)]
    t_b = fo(nfull + 2 * s["lid"]["panels"], 32 + len(sums_off), half=True)
    t_c = passes(nfull, s["cs"]["points"], grid(32), 200.0, 91.4, sums_off)
    n_inner = 3300
    t_fs = 3.83 * 3600 * (nfull / 5022.0) * (n_inner / 4044.0) * (32 / 21.0)
    rows.append(("J3 sw_full", "SOFTWIND QTF- and QTF+, FS integral, 32 freq. to 2.01 rad/s", t_b + t_c + t_fs + 600))
    # J4: OC4 bi-directional
    o = c["oc4_bi"]
    t = fo(o["hull"] + o["lid"]["panels"], 50) * 1.3 + passes(o["hull"], o["cs"]["points"], grid(50), 200.0, 20.0) * 1.3
    rows.append(("J4 oc4_bi", "OC4 QTF-, headings 0 and 30, all pairs, 50 freq. to 3.14", t + 600))
    # J5: OC3
    o = c["oc3_diff"]
    t = fo(o["hull"] + o["lid"]["panels"], 50) + passes(o["hull"], o["cs"]["points"], grid(50), 320.0, 120.0)
    rows.append(("J5 oc3_diff", "OC3 QTF-, 50 freq. to 3.14 rad/s (hydrosCal first)", t + 600))
    # cheaper options
    o = c["oc3_diff"]
    t = fo(o["hull"] + o["lid"]["panels"], 32) + passes(o["hull"], o["cs"]["points"], grid(32), 320.0, 120.0)
    rows.append(("J5b oc3 (option)", "OC3 as J5 but 32 freq. to 2.01 rad/s", t + 600))
    s = c["sw_diff"]
    t = fo(nfull + 2 * s["lid"]["panels"], 32, half=True) + passes(nfull, s["cs"]["points"], grid(32), 200.0, 91.4)
    rows.append(("J2b sw_diff (option)", "SOFTWIND QTF- as J2 but 32 freq. to 2.01", t + 600))
    out = []
    for name, what, sec in rows:
        h = sec / 3600.0
        lim = max(3, int(math.ceil(2.5 * h)))
        out.append(dict(job=name, what=what, wall_h=round(h, 2), core_h=round(CPUS * h, 1), sbu=round(1.5 * CPUS * h),
                        limit_h=lim, limit_core_h=CPUS * lim, limit_sbu=round(1.5 * CPUS * lim)))
    return out


def main(cases_json, png, md):
    c = json.load(open(cases_json))
    rows = estimate(c)
    base = [r for r in rows if "option" not in r["job"]]
    tot = dict(core_h=sum(r["core_h"] for r in base), sbu=sum(r["sbu"] for r in base),
               limit_core_h=sum(r["limit_core_h"] for r in base), limit_sbu=sum(r["limit_sbu"] for r in base))
    alt = [r for r in rows if r["job"] not in ("J2 sw_diff", "J5 oc3_diff")]
    tot_alt = dict(core_h=sum(r["core_h"] for r in alt), sbu=sum(r["sbu"] for r in alt),
                   limit_core_h=sum(r["limit_core_h"] for r in alt), limit_sbu=sum(r["limit_sbu"] for r in alt))
    with open(md, "w") as f:
        f.write("| Job | Content | Expected wall (h) | Core-h | SBU | Limit (h) | Core-h at limit | SBU at limit |\n")
        f.write("|---|---|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            f.write("| %s | %s | %.1f | %.0f | %d | %d | %d | %d |\n" % (r["job"], r["what"], r["wall_h"], r["core_h"],
                                                                     r["sbu"], r["limit_h"], r["limit_core_h"],
                                                                     r["limit_sbu"]))
        f.write("| **All (J1-J5)** | | | **%.0f** | **%d** | | **%d** | **%d** |\n" % (tot["core_h"], tot["sbu"],
                                                                                      tot["limit_core_h"], tot["limit_sbu"]))
        f.write("| All with the options J2b, J5b | | | %.0f | %d | | %d | %d |\n" % (tot_alt["core_h"], tot_alt["sbu"],
                                                                                    tot_alt["limit_core_h"],
                                                                                    tot_alt["limit_sbu"]))
    fig, ax = plt.subplots(figsize=(10, 4.8))
    names = [r["job"] for r in rows]
    y = range(len(rows))
    ax.barh(y, [r["limit_core_h"] for r in rows], color="#d9d9d9", label="limit (2.5 x expected, rounded up, at least 3 h)")
    ax.barh(y, [r["core_h"] for r in rows], color="#2f6db5", label="expected")
    for i, r in enumerate(rows):
        ax.text(r["limit_core_h"] + 5, i, "%.0f / %d core-h (%d / %d SBU)" % (r["core_h"], r["limit_core_h"], r["sbu"],
                                                                         r["limit_sbu"]), va="center", fontsize=8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("core-hours billed (fat_genoa, 24 CPUs per job); SBU = 1.5 x core-h")
    ax.set_title("Nemoh-paper cases on Snellius: J1-J5 %.0f core-h expected (%d SBU), %d at the limits (%d SBU)\n"
                 "with options J2b and J5b instead of J2 and J5: %.0f expected, %d at the limits"
                 % (tot["core_h"], tot["sbu"], tot["limit_core_h"], tot["limit_sbu"], tot_alt["core_h"],
                    tot_alt["limit_core_h"]), fontsize=9)
    ax.set_xlim(0, max(r["limit_core_h"] for r in rows) * 1.45)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    for r in rows:
        print(r)
    print("total", tot, "alt", tot_alt)


if __name__ == "__main__":
    main(*sys.argv[1:4])
