"""Accuracy gate G1 of the free-surface speed-ups, laptop part: the table G1_table.md and the figure G1_errors.png.

Usage: python make_g1_table.py   (writes both files next to this script)

Only numbers already measured at unit-test scale are used (no HAMS case). Sources:
- fs_inner_source_exact.txt (this folder; unit test fs_inner_source_matches_exact, Task 6, exact potential on the sharp m1
  hull, fixed body): the inner fields and the inner integral of the stencil and source routes, circle samples on the edge.
- The same test with the outer-edge partial waves from qtf_fs_circle = panels (fs_circle_modes graf = .true., a temporary
  test edit, Task 7, 2026-11-20): sum 3.581e-2, difference 1.782e-2, the same digits as with the samples.
- pw_from_panels_matches_samples (Task 5, 448-panel cylinder): the field on the circle against the exact direct-sum
  centroid field, samples 1.68e-6, panels 3.59e-9 (MacCamy-Fuchs state; re-run in Task 7 with the same values).
- fs_circle_panels_close (integration, qtf_fs_main): FS_inner + FS_outer of panels against samples, 4.715e-6 of max |FS|.
- Task 6 review, hull refinement of qtf_fs_main (220 against 960 hull panels, fixed body).
The stencil route's inner fields and inner integral do not use the circle (the partial waves enter only the outer
integral and, on the source route, the line integral on the outer edge), so their panels values are those of samples by
construction; the circle field does not depend on the inner route.
G1 rule (the accuracy gate of G1_table.md): a route passes a quantity when its error is at most today's + max(0.5 percentage points, 10 % of today's).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROUTES = ("stencil + samples (today)", "stencil + panels (item 3)", "source + samples (item 2)", "source + panels (2 + 3)")
SHORT = ("today", "item 3", "item 2", "items 2 + 3")

# quantity, unit note, errors in percent per route (ROUTES order), how each was obtained per route
ROWS = [
    ("phi, inner points 0.1-3.5 m from the waterline (worst ring)", [2.501, 2.501, 1.552, 1.552],
     ["measured", "= today, no circle", "measured", "= item 2, no circle"]),
    ("grad phi, same rings (worst: 0.1 m)", [4.918, 4.918, 5.838, 5.838],
     ["measured", "= today, no circle", "measured", "= item 2, no circle"]),
    ("inner integral, sum pair (1.4, 1.0) rad/s, max over psi_1,3,5", [2.049, 2.049, 3.581, 3.581],
     ["measured", "= today, no circle", "measured", "measured"]),
    ("inner integral, difference pair, max over psi_1,3,5", [2.284, 2.284, 1.782, 1.782],
     ["measured", "= today, no circle", "measured", "measured"]),
    ("field on the circle against the exact direct sum", [1.68e-4, 3.59e-7, 1.68e-4, 3.59e-7],
     ["measured", "measured", "= today, no inner route", "= item 3, no inner route"]),
]


def margin(today):
    return max(0.5, 0.1 * today)


def table():
    out = ["# Accuracy gate G1, laptop part (unit-test scale, exact references)", "",
           "Errors in % of the largest exact value (fields: of the ring; integral: of the three psi_k; circle: of the",
           "largest sample). Exact data: the sharp m1 hull (448 panels), fixed body, the exact potential on the panels",
           "(fs_inner_source_matches_exact, fs_inner_source_exact.txt in this folder) and the 448-panel cylinder for the",
           "circle (pw_from_panels_matches_samples). G1 margin = max(0.5 percentage points, 10 % of today's error); a route",
           "passes a quantity when its error is at most today's + margin. Made by make_g1_table.py; figure G1_errors.png.", "",
           "| Quantity | " + " | ".join(ROUTES) + " | G1 margin | Item 3 (panels) | Item 2 (source) | Items 2 + 3 |",
           "|---|" + "---|" * (len(ROUTES) + 4)]
    verdict = {1: True, 2: True, 3: True}
    for name, errs, how in ROWS:
        today = errs[0]
        m = margin(today)
        cells = []
        for e, h in zip(errs, how):
            note = "" if h == "measured" else " (%s)" % h
            cells.append(("%.3g" % e) + note)
        res = []
        for j in (1, 2, 3):
            ok = errs[j] <= today + m
            verdict[j] = verdict[j] and ok
            res.append("pass (%+.3g pp)" % (errs[j] - today) if ok else "**FAIL** (%+.3g pp)" % (errs[j] - today))
        out.append("| %s | %s | %.3g pp | %s |" % (name, " | ".join(cells), m, " | ".join(res)))
    out += ["", "![G1 errors](G1_errors.png)", "",
            "G1_errors.png: error of each quantity against its exact value, per route (log axis); the dashed line is the "
            "G1 limit of each quantity. The panels circle (hatched) changes only the circle field, which it puts on the "
            "exact series; the source route (orange) passes phi and the difference integral but fails grad phi at 0.1 m "
            "and the sum integral. Bars marked \"=\" in the table are equal by construction, not separate runs."]
    out += ["", "Verdict on the laptop numbers: item 3 (qtf_fs_circle = panels) %s; item 2 (qtf_fs_inner = source) %s;"
            " both together %s." % tuple("passes" if verdict[j] else "fails" for j in (1, 2, 3)), ""]
    out += [
        "## Without an exact laptop reference", "",
        "| Quantity | Measure | Value | Reading for G1 |",
        "|---|---|---|---|",
        "| FS_inner + FS_outer, qtf_fs_main (integration) | panels against samples, of max abs FS | 4.7e-6 | the panels error "
        "is within 4.7e-4 pp of the samples error (triangle inequality): far inside the 0.5 pp margin |",
        "| non-FS files, qtf_fs_main | panels against samples | line for line | unchanged |",
        "| FS_inner + FS_outer, qtf_fs_main (integration) | source against stencil, of max abs FS | 0.559 (free), 0.067-0.268 "
        "per mode (fixed) | not an error against an exact value, but 13-112 times the 0.5 pp G1 margin |",
        "| hull refinement 220 -> 960 panels, fixed body (Task 6 review) | change of FS_inner + FS_outer, of max abs FS | "
        "stencil 0.2-1.1 % (heave 4.9 %), source 1.5-24.9 % | the source route is not hull-converged; the stencil route is |",
        "| .11s/.12s FS part and totals, sharp m2 fixed and free, round m2 | against the 7c-2 exact references | Snellius | "
        "prepared in ../snellius (G1 item 3 job set), not run |", "",
        "## Decision", "",
        "- Item 2 (qtf_fs_inner = source) stays off (ruling 1, decided on these numbers): grad phi next to the waterline "
        "(5.84 % against 4.92 %, margin 0.5 pp) and the sum inner integral (3.58 % against 2.05 %) fail, and the route moves "
        "16-25 % with the hull mesh where the stencil route moves 0.3-1 %. It is also blocked with qtf_direct = 1, "
        "IRSP = 1 and ISOLV = 2. No Snellius case for it.",
        "- Item 3 (qtf_fs_circle = panels) passes every laptop quantity: it changes only the partial waves, which it puts "
        "on the exact series (3.6e-9 against 1.7e-6 for the samples), and the end-to-end FS terms move by 4.7e-6. It "
        "becomes the default only after the Snellius part (sharp m2 fixed and free, round m2, sum and difference, against "
        "the 7c-2 exact references; ../snellius/README.md).",
        "- Items 2 + 3 together: the panels circle does not change the source route's inner integral (3.581e-2 and "
        "1.782e-2 with either circle, measured), so the combination fails where item 2 fails.",
        "- The defaults are not changed here (stencil and samples stay); the defaults commit waits for the Snellius results.",
    ]
    return "\n".join(out) + "\n"


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Okabe-Ito colours (colour-blind safe), one per route, fixed order; hatching as a second cue
    colours = ["#0072B2", "#56B4E9", "#E69F00", "#D55E00"]
    hatches = ["", "//", "", "//"]
    labels = ["phi\n(near waterline)", "grad phi\n(near waterline)", "inner integral\nsum pair", "inner integral\ndifference pair",
              "field on\nthe circle"]
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    w = 0.2
    for j, route in enumerate(ROUTES):
        xs = [i + (j - 1.5) * w for i in range(len(ROWS))]
        vals = [row[1][j] for row in ROWS]
        ax.bar(xs, vals, w * 0.92, color=colours[j], hatch=hatches[j], edgecolor="white", linewidth=0.8, label=route)
    for i, row in enumerate(ROWS):
        today = row[1][0]
        lim = today + margin(today)
        ax.plot([i - 2 * w, i + 2 * w], [lim, lim], color="#333333", linewidth=1.2, linestyle="--")
    ax.plot([], [], color="#333333", linewidth=1.2, linestyle="--", label="gate limit: today + max(0.5 pp, 10 % of today)")
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 30)
    ax.set_xticks(range(len(ROWS)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("error against the reference (% of its largest value)")
    ax.grid(alpha=0.3, axis="y")
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)
    ax.set_title("Accuracy gate, laptop, sharp m1 hull, fixed body. Reference for phi, grad phi and the integrals: the\n"
                 "semi-analytical potential (make_fs_fields_exact_fixture.py, truncated_cylinder_reference.py; Garrett 1971;\n"
                 "Yeung 1981); for the circle: the direct Green-function sum over the hull (pw_from_panels_matches_samples)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "G1_errors.png"), dpi=150)


def main():
    with open(os.path.join(HERE, "G1_table.md"), "w", newline="\n") as f:
        f.write(table())
    figure()
    print("wrote G1_table.md and G1_errors.png")


if __name__ == "__main__":
    main()
