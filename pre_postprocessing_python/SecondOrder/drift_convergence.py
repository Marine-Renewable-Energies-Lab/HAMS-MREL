"""
Convergence study of the mean wave drift: generate the cases, write Slurm scripts, run one case, and collect a table.

Usage:
    python drift_convergence.py list [--tier full|small] [--markdown] [--control [full|cheap]]
    python drift_convergence.py generate <root> [--case NAME | --group GROUP] [--control [full|cheap]]
    python drift_convergence.py run <root> --case NAME --hams <path> [--hams-grn8 <path>] [--nemoh-build <path>] [--threads N]
                                [--control [full|cheap]]
    python drift_convergence.py calibrate <root> --hams <path> [--threads N]
    python drift_convergence.py slurm <root on the cluster> --repo <path on the cluster> [--out <local folder>]
                                [--tier full|small] [--account NAME] [--control [full|cheap]]
    python drift_convergence.py collect <root> [--out table.csv] [--control [full|cheap]] [--reference <root>]

Each case is one run of HAMS-MREL or Nemoh, under <root>/<name>. Its input is written when it is first run. A run writes
result.json with its wall time and peak memory.

The table gives F_x, F_y and M_z of the near-field method against the far-field method of the same run. On the bottom-mounted
cylinder it also gives F_x of both methods against the MacCamy-Fuchs value, where the exact F_y and M_z are zero.

The frequencies are set here before any result is seen, and stay as they are: a uniform grid with the three frequencies of the
fixtures added. 1.6 rad/s is 3% from the first irregular frequency of the cylinder, 1.548 rad/s, without a lid.

--control selects the study of the control-surface method (drift_methods = far near control). It runs the cases of groups A and B,
the boxes of group C, and the two m3 cases of the rounded cylinder that have the green route. A case with velocity_method
differencing takes cs_velocity = green, and one with velocity_method = source takes cs_velocity = source, so the near-field
column of every case is the one of the first study and the source solve is shared. `cheap` gives the m4 cases of the green route
only the three frequencies of the fixtures. `list --control` and `slurm --control` only print and write: nothing is submitted.
The table then has the columns Fx_ctrl_green or Fx_ctrl_src (and Fy, Mz) next to the near-field ones, and Fz, Mx, My of the
control-surface method against near-field of the same run and against near-field at m4 of the same hull and motion, the converged
reference. That m4 case is read from <root>, or from --reference (the results of the first study) when <root> does not hold it.
`collect --control` also prints the errors per case over the tolerance of the gates, |control - reference| over
max(1% |reference|, 0.001 rho g A^2 a^k), a = 10 m, k = 1 for forces and 2 for moments, and the Phase B gate: F_z, M_x, M_y
within that tolerance of near-field on the rounded cylinder at m3.

Nemoh runs as a separate program. Only its documented input and output files are used. Its source and binaries stay out of this
repository. --nemoh-build is the folder holding its built programs.

The estimates use the constants below, from `calibrate` on a laptop, and are checked against the wall time of a real run by
`collect`. The collector needs scipy for the MacCamy-Fuchs value: on Snellius, load a SciPy-bundle module or collect elsewhere.
"""

import argparse
import csv
import json
import math
import os
import re
import socket
import subprocess
import sys
import time

import make_drift_fixtures as fixtures
import make_mesh

# Frequencies in rad/s
FIXTURE_FREQUENCIES = fixtures.DRIFT_CYL_FREQUENCIES
WIDE_FREQUENCIES = (0.40, 0.60, 0.70, 0.80, 0.95, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80)
# The first two irregular frequencies of the cylinder of radius 10 and draft 10, from J0(ka) = 0 and J1(ka) = 0
IRREGULAR_FREQUENCIES = (1.10, 1.50, 1.548, 1.60, 1.939)
BOX_FINE_FREQUENCIES = (0.80, 1.00, 1.40)
# Nemoh needs a uniform grid. All of these lie below the first irregular frequency, and contain WIDE_FREQUENCIES up to 1.40
NEMOH_FREQUENCY_RANGE = (0.40, 1.40, 21)
# The QTF study (--qtf) runs HAMS and Nemoh on a finer grid, so that the QTF curves along omega_j = omega_i + dw are smooth
QTF_FREQUENCY_RANGE = (0.40, 1.40, 101)
# The QTF potential study (--qtf-potential): HAMS on QTF_FREQUENCY_RANGE with qtf_terms = body and qtf_max_dw, and Nemoh from
# 0.01 rad/s, since Nemoh takes the first-order solutions at the difference frequencies from its own frequency list
QTF_POTENTIAL_MAX_DW = 0.3
QTF_POTENTIAL_NEMOH_RANGE = (0.01, 1.40, 140)
# Nemoh.cal, QTF section: the contributing terms (1 DUOK, 2 DUOK + HASBO, 3 full) and the flags of the terms the
# post-processor adds to the total, (DUOK, HASBO, HASFS + ASYMP)
NEMOH_QTF_TERMS = {1: (1, 0, 0), 2: (1, 1, 0)}
NEMOH_STAGES = (("preProcessor", "preProc"), ("Solver", "solver"), ("postProcessor", "postProc"),
                (os.path.join("QTF", "PreProcessor"), "QTFpreProc"), (os.path.join("QTF", "Solver"), "QTFsolver"),
                (os.path.join("QTF", "PostProcessor"), "QTFpostProc"))
# Holds a Nemoh body still, against the wave forces of these cases: under 1e6 N
NEMOH_FIXED_STIFFNESS = 1.0e13

# Laptop, 8 threads, ifort 2022, seconds per frequency for two headings: A + B N^2 + C N^3. From `calibrate` at 1008 and 1792
# panels, with the logged 153 s per frequency at 4032 panels. The N^3 part is what the small sizes cannot fix on their own.
TIME_N0 = 0.5
TIME_N2 = 3.0e-6
TIME_N3 = 1.6e-9
# The source solve adds this. Its assembly and factorization are logged, 0.82, 3.77 and 16.5 s at 448, 1008 and 1792 panels
SOURCE_N2 = 3.5e-6
SOURCE_N3 = 0.9e-9
# ISOLV=2 against ISOLV=1, measured at 480 panels only
GMRES_FACTOR = 1.0
# Laptop time over Snellius time at 32 threads. See the job table for the evidence
SNELLIUS_SPEEDUP = 4.0
# Nemoh, one thread on the laptop, seconds per frequency, from runs of two frequencies at 448, 1008 and 1792 panels. They do
# not fix the N^3 part, which is the factorization of N panels at a rate similar to the HAMS one on one thread.
NEMOH_N2 = 6.0e-6
NEMOH_N3 = 0.5e-9
NEMOH_QTF_N2 = 6.2e-6
NEMOH_BYTES = 102.0   # per N^2, measured at 1792 panels
# Snellius time of a Nemoh run over nemoh_seconds: nemoh_sharp_*_m2 (1792 panels, 140 frequencies, DUOK + HASBO) took
# 2 h 34 min on one thread (job 29017289) against 1.635 h from nemoh_seconds. One point only: an N^3 share above the modelled
# one would make larger runs slower, which the time limit of twice the estimate covers
NEMOH_SNELLIUS_FACTOR = 2.567 / 1.635

# Bytes per N^2 per symmetry image, from the ALLOCATE statements of HAMS_Prog.f90 and SurfaceVelocity.f90
BYTES_AMAT = 16.0
BYTES_GREEN = {4: 48.0, 8: 96.0}   # four parts, each of a complex and a real array
BYTES_SOURCE = 64.0   # lu 16, gradient 3 x 16
BYTES_GMRES = 320.0   # measured at 480 panels, where the H-matrix compresses least
BYTES_BASE = 0.02e9   # all the rest, from `calibrate`

# Slurm settings. 48 CPUs are charged for this request, with a weight of 1.5 on this partition (sacct: billing=72)
PARTITION = "fat_genoa"
CPUS = 32
MEMORY_GB = 64
CHARGED_CPUS = 48
SBU_PER_JOB_HOUR = 72.0
# Nemoh runs on one thread per case, side by side in one allocation, which asks for the 24 CPUs that fat_genoa bills at least,
# and little memory (a few GB per run)
NEMOH_CPUS = 24
NEMOH_MEMORY_GB = 24
# ifx 2024.2 with -heap-arrays 0 stops with SIGSEGV when it frees the temporary of an array constructor of a type with allocatable
# parts, as in far_field_pair. A small program with no HAMS code does the same. So the jobs leave that flag out, and give the
# threads a large stack instead. A stack overflow stops the run, so it cannot change a number silently.
IFX_FLAGS = "-qmkl -qopenmp"


class Case:
    """One run. tier is "both", "full" or "small": whether it belongs to the full study, the smaller option, or both."""

    def __init__(self, name, group, question, hull, multiple, motion="fixed", velocity="differencing", waterline="extrapolate",
                 frequencies=WIDE_FREQUENCIES, headings=(0.0, 30.0), grn=4, isolv=1, irsp=0, mean_drift=1, tier="both",
                 solver="hams", where="snellius", control=None, qtf=0, nemoh_range=None, qtf_terms=None, qtf_max_dw=None,
                 nemoh_contrib=1, edge_power=None, qtf_edge_weights=None, extra_options=None, damping=None):
        self.name = name
        self.group = group
        self.question = question
        self.hull = hull
        self.multiple = multiple
        self.motion = motion
        self.velocity = velocity
        self.waterline = waterline
        self.frequencies = frequencies
        self.headings = headings
        self.grn = grn
        self.isolv = isolv
        self.irsp = irsp
        self.mean_drift = mean_drift
        self.tier = tier
        self.solver = solver
        self.where = where
        self.control = control   # None, "green" or "source": the cs_velocity of drift_methods = control
        self.qtf = qtf           # qtf of SecondOrder.in: 0 off, 3 both quadratic QTFs
        self.nemoh_range = nemoh_range or NEMOH_FREQUENCY_RANGE   # (low, high, count) of a Nemoh run
        self.qtf_terms = qtf_terms     # qtf_terms of SecondOrder.in: None (not written, quadratic) or "body"
        self.qtf_max_dw = qtf_max_dw   # qtf_max_dw of SecondOrder.in, or None (not written)
        self.nemoh_contrib = nemoh_contrib   # Nemoh's QTF terms: 1 DUOK, 2 DUOK + HASBO (NEMOH_QTF_TERMS)
        self.edge_power = edge_power   # sharp hull only: power-law rows toward the bilge (make_mesh.cylinder), None for cosine
        self.qtf_edge_weights = qtf_edge_weights   # qtf_edge_weights of SecondOrder.in: None (not written), "on" or "off"
        self.extra_options = extra_options or {}   # further SecondOrder.in keys, written as given
        self.damping = damping or {}   # sharp and oc4 hulls: external linear damping {mode 1..6: value} (Hydrostatic.in)
        if self.damping and hull not in ("sharp", "oc4"):
            raise ValueError("%s: damping is written for the sharp and oc4 hulls only, not %s" % (name, hull))


def in_tier(case, tier):
    return case.tier == "both" or case.tier == tier


# Hulls: water depth, panels at multiple 1, and lid panels at multiple 1 where a lid exists
HULLS = {
    "sharp": dict(depth=100.0, panels=448, lid=192),
    "round": dict(depth=100.0, panels=480, lid=192),
    "bottom": dict(depth=10.0, panels=256),
    "box": dict(depth=50.0, panels=440),
    "round_box": dict(depth=50.0, panels=500),
    "corner_box": dict(depth=50.0, panels=440),
    # OC4 semisubmersible (make_mesh.OC4_MESH): the panels per multiple, since the transition rows do not scale as m^2. Its
    # C and F (CSMesh_1.pnl, FSMesh_1.pnl, cs_shape = file) are the same at every multiple: C 1400 panels, F per multiple.
    # cf_points counts the panels of C and F only: B2 evaluates no field at the segments of Gamma_c (b2_control_setup in
    # SecondOrderDriver.f90), so the 1416 points of B2_PASS_SECONDS are likewise C and F points only.
    # Its Nemoh runs have allocations per multiple (nemoh_batches).
    # A fixed Nemoh body needs a stiffer hold than NEMOH_FIXED_STIFFNESS: wave moments of order 1e9 N m per metre of amplitude
    "oc4": dict(depth=make_mesh.OC4_DEPTH, panels={1: 1224, 2: 4608}, cf_points={1: 1400 + 750, 2: 1400 + 1384},
                fixed_stiffness=1.0e15, nemoh_batches=True),
    # The same meshes with the matrices of the OpenFAST WAMIT run (make_mesh.OC4W_*, write_oc4 model "wamit"): HAMS only
    "oc4w": dict(depth=make_mesh.OC4_DEPTH, panels={1: 1224, 2: 4608}, cf_points={1: 1400 + 750, 2: 1400 + 1384}),
}
# Both radii are for a box 20 x 10 x 5. 1.5 is the bilge radius of the rounded cylinder over its draft, times the draft here
ROUND_BOX_MESH = dict(lx=20.0, ly=10.0, draft=5.0, corner_radius=2.5, bilge_radius=1.5, nx=16, ny=6, ncorner=4, nz=3, narc=2)
CORNER_BOX_MESH = dict(lx=20.0, ly=10.0, draft=5.0, corner_radius=2.5, bilge_radius=0.0, nx=16, ny=6, ncorner=4, nz=4, narc=0)
BOTTOM_RADIUS = fixtures.BOTTOM_MOUNTED_MESH["radius"]
BOTTOM_DEPTH = fixtures.BOTTOM_MOUNTED_MESH["draft"]


def scaled(mesh, keys, multiple):
    out = dict(mesh)
    for key in keys:
        out[key] = int(round(mesh[key] * multiple))
    return out


def panels_of(case):
    panels = HULLS[case.hull]["panels"]
    if isinstance(panels, dict):
        return panels[case.multiple]
    return int(round(panels * case.multiple ** 2))


def lid_panels_of(case):
    return int(round(HULLS[case.hull]["lid"] * case.multiple ** 2)) if case.irsp else 0


def tag(multiple):
    return "m%g" % multiple


def study_cases():
    """All the cases, by the question they answer. Combinations that answer none of them are not here."""
    cases = []
    short = {"differencing": "diff", "source": "src"}

    # A: truncated cylinder, sharp and rounded bilge. Which velocity method, by hull and by frequency, whether the rounded hull's
    # term II error goes to zero, and the sharp free hull at 0.70 rad/s. Multiples 1, 2, 4 double the panels per side on the wide
    # range. Multiple 3 with the source solve completes the table of 4032 and 4320 panels at the frequencies of the fixtures.
    for hull in ("sharp", "round"):
        for motion in ("fixed", "free"):
            for velocity in ("differencing", "source"):
                for m in (1, 2, 4):
                    cases.append(Case("a_%s_%s_%s_%s" % (hull, motion, short[velocity], tag(m)), "A", "1,2,7", hull, m, motion,
                                      velocity, tier="full" if m == 4 else "both"))
            cases.append(Case("a_%s_%s_src_m3" % (hull, motion), "A", "1,2", hull, 3, motion, "source",
                              frequencies=FIXTURE_FREQUENCIES))

    # B: bottom-mounted cylinder, the only exact reference. Motion plays no part, and heading 0 is enough for a circle
    for velocity in ("differencing", "source"):
        for m in (1, 2, 4):
            cases.append(Case("b_bottom_%s_%s" % (short[velocity], tag(m)), "B", "1,7", "bottom", m, "fixed", velocity,
                              headings=(0.0,)))

    # C: box at heading 30, the residual of M_z and F_y. The sharp box is the second sharp hull for the velocity method. The
    # rounded box asks whether having no edges removes the residual, and corners alone whether the vertical edges carry it.
    for hull, motions in (("box", ("fixed", "free")), ("round_box", ("fixed", "free"))):
        for velocity in ("differencing", "source"):
            for motion in motions:
                for m in (1, 2):
                    cases.append(Case("c_%s_%s_%s_%s" % (hull, motion, short[velocity], tag(m)), "C", "3,1,7", hull, m, motion,
                                      velocity, headings=(30.0,)))
            cases.append(Case("c_%s_fixed_%s_m4" % (hull, short[velocity]), "C", "3,7", hull, 4, "fixed", velocity,
                              frequencies=BOX_FINE_FREQUENCIES, headings=(30.0,), tier="full"))
            cases.append(Case("c_%s_fixed_%s_m3" % (hull, short[velocity]), "C", "3,7", hull, 3, "fixed", velocity,
                              frequencies=BOX_FINE_FREQUENCIES, headings=(30.0,), tier="small"))
    for m in (1, 2):
        cases.append(Case("c_corner_box_fixed_diff_%s" % tag(m), "C", "3", "corner_box", m, "fixed", headings=(30.0,)))

    # D: settings. Each is compared with the same mesh and frequency in group A
    for hull, velocity, m in (("round", "differencing", 2), ("round", "source", 2), ("sharp", "differencing", 2),
                              ("sharp", "source", 2), ("round", "source", 3)):
        cases.append(Case("d_grn8_%s_free_%s_%s" % (hull, short[velocity], tag(m)), "D", "6", hull, m, "free", velocity,
                          frequencies=FIXTURE_FREQUENCIES, grn=8))
    # The source solve refuses ISOLV=2, so this is differencing only
    cases.append(Case("d_base_round_free_diff_m3", "D", "6", "round", 3, "free", frequencies=FIXTURE_FREQUENCIES))
    for hull, m in (("round", 1), ("round", 2), ("round", 3), ("sharp", 2)):
        cases.append(Case("d_gmres_%s_free_%s" % (hull, tag(m)), "D", "6", hull, m, "free", frequencies=FIXTURE_FREQUENCIES,
                          isolv=2))
    # Near the irregular frequencies: the source solve refuses IRSP=1, so it is run with IRSP=0 to see what it does there
    for hull, m in (("round", 1), ("round", 2), ("sharp", 2)):
        for label, velocity, irsp in (("diff_irsp0", "differencing", 0), ("diff_irsp1", "differencing", 1),
                                      ("src_irsp0", "source", 0)):
            cases.append(Case("d_%s_%s_free_%s" % (label, hull, tag(m)), "D", "6", hull, m, "free", velocity,
                              frequencies=IRREGULAR_FREQUENCIES, irsp=irsp))

    # E: a small check of the field waterline. Each is compared with its extrapolate run in group A
    for hull in ("round", "sharp"):
        cases.append(Case("e_field_%s_free_m2" % hull, "E", "5", hull, 2, "free", waterline="field"))

    # N: code comparison with Nemoh. The HAMS runs of group A serve on the frequencies shared with the Nemoh grid. These add
    # what group A does not have: differencing on the fixed rounded hull at multiple 3, and the heading pairs.
    cases.append(Case("n_round_fixed_diff_m3", "N", "8", "round", 3, "fixed", frequencies=FIXTURE_FREQUENCIES))
    for m in (1, 2):
        cases.append(Case("n_pair_round_free_%s" % tag(m), "N", "8", "round", m, "free", mean_drift=2,
                          frequencies=(0.70, 0.95, 1.10), tier="full"))

    # Nemoh itself. Multiple 1 is short enough for the laptop. Finer runs need a Nemoh build on Snellius
    nemoh = dict(group="N", question="8", solver="nemoh", headings=(0.0,))
    for hull, multiples in (("round", (1, 2, 3)), ("sharp", (1, 2))):
        for motion in ("fixed", "free"):
            for m in multiples:
                tier = "both" if m == 1 or (hull == "round" and m == 2) else "full"
                cases.append(Case("nemoh_%s_%s_%s" % (hull, motion, tag(m)), hull=hull, multiple=m, motion=motion, tier=tier,
                                  where="laptop" if m == 1 else "snellius", **nemoh))
    for m in (1, 2):
        pair = dict(nemoh, headings=(0.0, 30.0))
        cases.append(Case("nemoh_pair_round_free_%s" % tag(m), hull="round", multiple=m, motion="free", mean_drift=2, tier="full",
                          where="laptop" if m == 1 else "snellius", **pair))

    names = [case.name for case in cases]
    if len(names) != len(set(names)):
        raise ValueError("two cases have the same name")
    return cases


# The cases that carry the control-surface method. The m3 cases of the rounded cylinder with differencing are in groups D and N
CONTROL_EXTRA = ("d_base_round_free_diff_m3", "n_round_fixed_diff_m3")


def control_cases(subset="full"):
    """The cases of the control-surface study. Each is a case of the first study with `far near control` in its SecondOrder.in,
    so its near-field and far-field columns are those of the first study. The differencing cases take the green route and the
    source cases take the source route. subset "cheap" gives the m4 cases of the green route the three frequencies of the fixtures."""
    chosen = []
    for case in study_cases():
        if not in_tier(case, "full") or case.solver != "hams":
            continue
        if case.group in ("A", "B") or (case.group == "C" and case.hull in ("box", "round_box")) or case.name in CONTROL_EXTRA:
            chosen.append(case)
    # The sharp cylinder has no differencing case at m3: the first study ran the source solve there only
    for motion in ("fixed", "free"):
        chosen.append(Case("a_sharp_%s_diff_m3" % motion, "A", "1,2", "sharp", 3, motion, "differencing",
                           frequencies=FIXTURE_FREQUENCIES))
    for case in chosen:
        case.control = "source" if case.velocity == "source" else "green"
        if subset == "cheap" and case.control == "green" and case.multiple == 4 and len(case.frequencies) > len(FIXTURE_FREQUENCIES):
            case.frequencies = FIXTURE_FREQUENCIES
    return chosen


def qtf_cases():
    """The cases of the QTF study (--qtf): quadratic QTFs (qtf = 3) on the grid of the Nemoh runs, heading 0, on the
    meshes of the first study, so the diagonal is the mean drift of the a_* cases and Nemoh's nemoh_* runs have the same
    panels. Velocity by source solve (the method that agreed with Nemoh and the exact reference for F_z), and the
    differencing velocity once, at m2 of the free sharp hull. `collect --qtf` collects their mean drift only;
    qtf_compare.py reads the QTF files in place, under q_*/output/Wamit_format."""
    grid = tuple(round(w, 4) for w in nemoh_frequencies(QTF_FREQUENCY_RANGE))
    chosen = []
    for motion in ("fixed", "free"):
        for m in (1, 2, 3, 4):
            chosen.append(Case("q_sharp_%s_src_%s" % (motion, tag(m)), "Q", "4", "sharp", m, motion, "source",
                               frequencies=grid, headings=(0.0,), qtf=3))
    for m in (1, 2):
        chosen.append(Case("q_round_free_src_%s" % tag(m), "Q", "4", "round", m, "free", "source", frequencies=grid,
                           headings=(0.0,), qtf=3))
    chosen.append(Case("q_sharp_free_diff_m2", "Q", "4", "sharp", 2, "free", "differencing", frequencies=grid,
                       headings=(0.0,), qtf=3))
    # Nemoh's DUOK part on the same grid and meshes
    for motion in ("fixed", "free"):
        for m in (1, 2):
            chosen.append(Case("nemoh_sharp_%s_%s" % (motion, tag(m)), "Q", "4", "sharp", m, motion, headings=(0.0,),
                               solver="nemoh", nemoh_range=QTF_FREQUENCY_RANGE))
    return chosen


# The default control surface of the free sharp hull at m1 (CSMesh_1.pnl of q_sharp_free_src_m1: radius 17.84, bottom
# -17.84, so the margins are 7.84 m = 4 x the median hull panel size 1.96 m)
CS_FIXED_M1 = dict(cs_margin=7.84, cs_bottom_margin=7.84, cs_panel_size=1.96)

# External linear damping of the damped free sharp hull: B33 = 0.05 x 2 sqrt((M + A33) C33) and B55 = 0.05 x 2
# sqrt((I55 + A55) C55), moments about the origin as in Hydrostatic.in. M = 3.12145e6 kg, I55 = 1.81584e8 kg m^2,
# C33 = 3.06109e6 N/m, C55 = 7.60372e7 N m (Hydrostatic.in of the m1 case); A33 = 1.67200e6 kg at 0.80 rad/s (the heave
# peak) and A55 = 5.51253e7 kg m^2 at 0.77 rad/s (the coupled surge-pitch peak), from
# truncated_cylinder_reference.added_mass_damping. So B55 takes A55 at the coupled peak, not at the uncoupled pitch
# frequency sqrt(C55 / (I55 + A55)) (about 0.57 rad/s): it is "5 % of critical" by this formula only, not the critical
# damping of an uncoupled pitch mode. At 0.77 rad/s it is 1.4 times (heave) and 2.1 times (pitch) the radiation damping
# rho omega B (2.75e5 N s/m, 6.50e6 N m s; truncated_cylinder_reference check damped_motions, LIGHT_RESOLUTION). A choice
# that damps the resonances, not a model of viscous damping.
# truncated_cylinder_reference reads it for its motion "free_damped".
SHARP_DAMPING = {3: 3.8306e5, 5: 1.3416e7}

# Grid of the q_oc4w cases: the omega of the WAMIT .12d of OpenFAST (marin_semi.12d: 0.25 to 3.00 rad/s in steps of 0.05) up
# to 1.90, below the first irregular frequency of OC4, about 1.98 rad/s (upper columns, J0 zero under a deep lid [derived]):
# velocity_method = source rejects IRSP = 1, so the grid stops there. 34 frequencies, 33 distinct differences
OC4W_FREQUENCIES = tuple(round(0.25 + 0.05 * i, 4) for i in range(34))

# Powers of the edge-graded meshes of the QTF potential study (make_mesh.cylinder edge_power)
EDGE_POWERS = (2.5, 3)


def qtf_potential_cases():
    """The cases of the QTF potential study (--qtf-potential): the sharp hull of qtf_cases, fixed and free, source velocity,
    m1 to m4, with qtf_terms = body (QTF_Indirect.11d and QTF_PotentialTerms.txt) on QTF_FREQUENCY_RANGE, qtf_max_dw =
    QTF_POTENTIAL_MAX_DW, heading 0, difference QTF only (qtf = 1); and Nemoh with DUOK + HASBO (body forcing, no free
    surface integral) at m1 and m2, from 0.01 rad/s (QTF_POTENTIAL_NEMOH_RANGE). The exact reference is
    truncated_cylinder_reference.py qtf-potential; qtf_compare.py --potential compares them. The OC4 semisubmersible
    (q_oc4_*, nemoh_oc4_*) has no exact reference: qtf_compare.py --against-nemoh compares HAMS with Nemoh on it. q_oc4w_*
    is OC4 with the matrices of the OpenFAST WAMIT run: qtf_compare.py --against-wamit compares it with that run's .12d."""
    grid = tuple(round(w, 4) for w in nemoh_frequencies(QTF_FREQUENCY_RANGE))
    chosen = []
    for motion in ("fixed", "free"):
        for m in (1, 2, 3, 4):
            chosen.append(Case("q_sharp_%s_src_%s" % (motion, tag(m)), "P", "5", "sharp", m, motion, "source",
                               frequencies=grid, headings=(0.0,), qtf=1, qtf_terms="body",
                               qtf_max_dw=QTF_POTENTIAL_MAX_DW))
    for motion in ("fixed", "free"):
        for m in (1, 2):
            chosen.append(Case("nemoh_sharp_%s_%s" % (motion, tag(m)), "P", "5", "sharp", m, motion, headings=(0.0,),
                               solver="nemoh", nemoh_range=QTF_POTENTIAL_NEMOH_RANGE, nemoh_contrib=2))
    # Edge grading: the free hull with power-law rows toward the sharp bilge (EDGE_POWERS), the same panel count, with the
    # sharp-edge weights off and on (q_sharp_free_g<q>[w]_m<k>, 2.5 as g2p5). The cosine rows of the
    # other cases are close to q = 2, so the powers start above it
    for q in EDGE_POWERS:
        for weights in ("off", "on"):
            for m in (1, 2):
                name = "q_sharp_free_g%s%s_%s" % (("%g" % q).replace(".", "p"), "w" if weights == "on" else "", tag(m))
                chosen.append(Case(name, "P", "5", "sharp", m, "free", "source", frequencies=grid, headings=(0.0,), qtf=1,
                                   qtf_terms="body", qtf_max_dw=QTF_POTENTIAL_MAX_DW, edge_power=q,
                                   qtf_edge_weights=weights))
    # Fixed control surface: C and F of the default at m1 (CS_FIXED_M1) at the finer hulls, where the default C moves in
    # (margin 4 x the median hull panel size) and refines with the hull, so its field passes cost N_hull x N_C
    for m in (2, 3):
        chosen.append(Case("q_sharp_free_fixc_%s" % tag(m), "P", "5", "sharp", m, "free", "source", frequencies=grid,
                           headings=(0.0,), qtf=1, qtf_terms="body", qtf_max_dw=QTF_POTENTIAL_MAX_DW,
                           extra_options=CS_FIXED_M1))
    # Damped free hull (SHARP_DAMPING; the reference's motion "free_damped"): m1 with the default C, m2, m3 with C fixed
    for m in (1, 2, 3):
        chosen.append(Case("q_sharp_free_damp_%s" % tag(m), "P", "5", "sharp", m, "free", "source", frequencies=grid,
                           headings=(0.0,), qtf=1, qtf_terms="body", qtf_max_dw=QTF_POTENTIAL_MAX_DW,
                           extra_options=CS_FIXED_M1 if m > 1 else None, damping=SHARP_DAMPING))
    # OC4 semisubmersible (make_mesh.write_oc4), HAMS against Nemoh on the same mesh, mass, CM and mooring. Four separate
    # waterlines are not star-shaped, so C and F come from files (CSMesh_1.pnl, FSMesh_1.pnl: cs_shape = file), the same C at
    # both meshes. Surge (about 0.06 rad/s) and pitch (about 0.3) lie well below the grid; heave at about 0.37 to 0.43 rad/s
    # [derived from the restoring and a heave-plate added mass of 0.7e7 to 1.4e7 kg] may reach its lowest frequencies
    for motion in ("fixed", "free"):
        for m in (1, 2):
            chosen.append(Case("q_oc4_%s_%s" % (motion, tag(m)), "P", "5", "oc4", m, motion, "source", frequencies=grid,
                               headings=(0.0,), qtf=1, qtf_terms="body", qtf_max_dw=QTF_POTENTIAL_MAX_DW,
                               extra_options=dict(cs_shape="file")))
    for motion in ("fixed", "free"):
        for m in (1, 2):
            chosen.append(Case("nemoh_oc4_%s_%s" % (motion, tag(m)), "P", "5", "oc4", m, motion, headings=(0.0,),
                               solver="nemoh", nemoh_range=QTF_POTENTIAL_NEMOH_RANGE, nemoh_contrib=2))
    # OC4 matched to the WAMIT run of OpenFAST (marin_semi.12d): the same meshes, C and F, with the full-system mass, the
    # mooring and the hydrostatic restoring of that run (make_mesh.oc4w_matrices), free body only, on the .12d grid up to
    # 1.90 rad/s (OC4W_FREQUENCIES), every pair (qtf_max_dw = 0). qtf_compare.py --against-wamit compares QTF_Indirect.11d
    # and QTF_Quadratic.10d with the .12d
    for m in (1, 2):
        chosen.append(Case("q_oc4w_free_%s" % tag(m), "P", "5", "oc4w", m, "free", "source", frequencies=OC4W_FREQUENCIES,
                           headings=(0.0,), qtf=1, qtf_terms="body", extra_options=dict(cs_shape="file", qtf_max_dw=0)))
    return chosen


def qtf_flag(qtf):
    """The command-line flag of a study selected by cases_of's qtf: " --qtf", " --qtf-potential" or empty."""
    return " --qtf-potential" if qtf == "potential" else " --qtf" if qtf else ""


def study_of(args):
    """cases_of's qtf from the parsed arguments: "potential" (--qtf-potential), True (--qtf) or False."""
    return "potential" if getattr(args, "qtf_potential", False) else getattr(args, "qtf", False)


def cases_of(control, qtf=False):
    """The cases selected by --control (the control-surface study), --qtf (qtf True, the QTF study) or --qtf-potential
    (qtf "potential"); otherwise the first study."""
    if qtf == "potential":
        return qtf_potential_cases()
    if qtf:
        return qtf_cases()
    return control_cases(control) if control else study_cases()


def nemoh_frequencies(frequency_range=None):
    low, high, count = frequency_range or NEMOH_FREQUENCY_RANGE
    return [low + (high - low) * i / (count - 1) for i in range(count)]


def laptop_seconds(case):
    """Estimated laptop time of a HAMS run at 8 threads."""
    n = panels_of(case) + lid_panels_of(case)
    per_frequency = TIME_N0 + TIME_N2 * n ** 2 + TIME_N3 * n ** 3
    if case.isolv == 2:
        per_frequency *= GMRES_FACTOR
    if case.velocity == "source":
        per_frequency += SOURCE_N2 * n ** 2 + SOURCE_N3 * n ** 3
    seconds = per_frequency * len(case.frequencies)
    if case.qtf_terms == "body":
        # one radiation-only solve per difference frequency up to qtf_max_dw, taken at the cost of a frequency [inferred]
        seconds += (TIME_N0 + TIME_N2 * n ** 2 + TIME_N3 * n ** 3) * len(difference_frequencies(case))
    return seconds


def difference_frequencies(case):
    """The distinct omega_i - omega_j > 0 of a case's grid within its qtf_max_dw (all if None): the frequencies of the
    radiation-only solves of qtf_terms = body."""
    top = case.qtf_max_dw if case.qtf_max_dw else float("inf")
    return sorted({round(b - a, 6) for a in case.frequencies for b in case.frequencies if 0.0 < b - a <= top * (1 + 1e-9)})


def nemoh_seconds(case):
    """Estimated single thread time of a Nemoh run: the first-order solver, and the QTF stages."""
    n = panels_of(case)
    count = case.nemoh_range[2]
    problems = (len(case.headings) + 6) / 7.0
    return count * problems * (NEMOH_N2 * n ** 2 + NEMOH_N3 * n ** 3 + NEMOH_QTF_N2 * n ** 2)


# The field pass of qtf_b2 = control on Snellius at 32 threads: 25 s at 1792 hull panels and 1416 points on C and F
# (q_sharp_free_fixc_m2), taken proportional to the product [inferred from the Green calls]. One pass per frequency and per
# difference frequency. Counted for the hulls whose C and F are fixed by files (HULLS cf_points); for the others it is not
# modelled here
B2_PASS_SECONDS = 25.0 / (1792 * 1416)


def b2_snellius_hours(case):
    points = HULLS[case.hull].get("cf_points")
    if case.qtf_terms != "body" or not points:
        return 0.0
    passes = len(case.frequencies) + len(difference_frequencies(case))
    return B2_PASS_SECONDS * panels_of(case) * points[case.multiple] * passes / 3600.0


def wall_hours(case):
    """Estimated wall time on Snellius, or on the laptop for a case that runs there."""
    if case.solver == "nemoh":
        return nemoh_seconds(case) / 3600.0 * (NEMOH_SNELLIUS_FACTOR if case.where == "snellius" else 1.0)
    return laptop_seconds(case) / SNELLIUS_SPEEDUP / 3600.0 + control_snellius_hours(case) + b2_snellius_hours(case)


def memory_gb(case):
    """Estimated peak memory, from the arrays that grow with N^2."""
    n = panels_of(case)
    total = n + lid_panels_of(case)
    if case.solver == "nemoh":
        return (BYTES_BASE + NEMOH_BYTES * n ** 2) / 1e9
    size = BYTES_AMAT * total ** 2 + BYTES_GREEN[case.grn] * n ** 2 + BYTES_BASE
    if case.irsp:
        size += BYTES_AMAT * n ** 2 + BYTES_GREEN[case.grn] * n * lid_panels_of(case)
    if case.velocity == "source":
        size += BYTES_SOURCE * n ** 2
    if case.isolv == 2:
        size += BYTES_GMRES * total ** 2
    if case.control:
        n_c, n_g = control_surface_counts(case)
        size += BYTES_CONTROL_POINT * (n_c + n_g) * len(case.headings)
    return size / 1e9

# Control-surface fields. The surface C is a cylinder around each body (cs_shape = cylinder, the default), of radius the largest
# horizontal extent of the hull plus 4 median panel sizes, with panels of the median hull panel size. Its panel count is modelled
# from the hull area, and fixed to the counts logged by the fixtures at one mesh (the `control_*` fixture logs, test/data/drift_cyl/README.md). The field time is
# proportional to the Green calls [derived]: (6 N_C + 5 N_Gamma_c) N_hull for the green route, and (N_C + N_Gamma_c) N_hull
# for the source route. Its rate per call is fixed to the logged field time of the same hull and water depth.
CS_MARGIN_PANELS = 4.0
# Hull: wall area, largest horizontal extent, draft (metres), whether the wall reaches the sea bed (then C has no bottom), and the
# reference fixture whose log fixes the model. The rounded hulls take the factors of the sharp ones with their own area and extent
CS_HULLS = {
    "sharp": dict(area=942.5, radius=10.0, draft=10.0, bed=False, ref="cyl"),
    "round": dict(area=857.5, radius=10.0, draft=10.0, bed=False, ref="cyl"),
    "bottom": dict(area=314.2, radius=5.0, draft=10.0, bed=True, ref="bottom"),
    "box": dict(area=500.0, radius=11.18, draft=5.0, bed=False, ref="box"),
    "round_box": dict(area=500.0, radius=10.4, draft=5.0, bed=False, ref="box"),
}
# Logged by the fixtures: hull panels, panels on C, segments on Gamma_c, and the mean field time per frequency in seconds on the
# laptop with 4 threads. cyl: control_cyl_free (2 headings, 45.3, 39.3 and 33.7 s). box: control_box (30.8 and 39.0 s). bottom:
# control_bottom (33.1 and 24.8 s). The source route was logged on cyl only: 3.7, 4.7 and 3.5 s
CS_REFERENCE = {
    "cyl": dict(hull="sharp", panels=448, n_c=1160, n_g=58, green_s=39.4, source_s=3.95),
    "box": dict(hull="box", panels=440, n_c=1804, n_g=82, green_s=34.9, source_s=None),
    "bottom": dict(hull="bottom", panels=1024, n_c=1292, n_g=76, green_s=29.0, source_s=None),
}
# The source field time over the green one, from the cyl logs. Taken for the other hulls, which were logged on the green route only
CS_SOURCE_OVER_GREEN = CS_REFERENCE["cyl"]["source_s"] / CS_REFERENCE["cyl"]["green_s"]
# Bytes per point for the point list, and per point and heading for the fields (a complex potential and gradient)
BYTES_CONTROL_POINT = 164.0


def _surface_model(hull, panels):
    """Unscaled panels on C and segments on Gamma_c of a cylinder surface around the hull."""
    h = CS_HULLS[hull]
    size2 = h["area"] / panels
    size = math.sqrt(size2)
    radius = h["radius"] + CS_MARGIN_PANELS * size
    depth = h["draft"] + (0.0 if h["bed"] else CS_MARGIN_PANELS * size)
    area = 2.0 * math.pi * radius * depth + (0.0 if h["bed"] else math.pi * radius ** 2)
    return area / size2, 2.0 * math.pi * radius / size


def control_surface_counts(case):
    """Panels on C and segments on Gamma_c of a case. Modelled, and equal to the log of the fixture at its reference mesh."""
    ref = CS_REFERENCE[CS_HULLS[case.hull]["ref"]]
    model_c, model_g = _surface_model(case.hull, panels_of(case))
    ref_c, ref_g = _surface_model(ref["hull"], ref["panels"])
    return model_c * ref["n_c"] / ref_c, model_g * ref["n_g"] / ref_g


def control_laptop_seconds(case):
    """Field time of the control method on the laptop with 4 threads: the logged rate per Green call, times the calls."""
    if not case.control:
        return 0.0
    ref = CS_REFERENCE[CS_HULLS[case.hull]["ref"]]
    n_c, n_g = control_surface_counts(case)
    n = panels_of(case)
    if case.control == "green":
        calls = (6.0 * n_c + 5.0 * n_g) * n
        rate = ref["green_s"] / ((6.0 * ref["n_c"] + 5.0 * ref["n_g"]) * ref["panels"])
    else:
        calls = (n_c + n_g) * n
        rate = CS_SOURCE_OVER_GREEN * ref["green_s"] / ((ref["n_c"] + ref["n_g"]) * ref["panels"])
    return calls * rate * len(case.frequencies)


def control_snellius_hours(case):
    """The same on Snellius at 32 threads, at the speed-up of the first study. That speed-up is for 8 laptop threads, and the
    fields were logged at 4, so this is on the high side."""
    return control_laptop_seconds(case) / SNELLIUS_SPEEDUP / 3600.0



def write_hams_hull(case, path):
    """Write the hull of a case at the folder path, in the HAMS format, and any lid. Returns the mesh."""
    mesh_path = os.path.join(path, "HullMesh.pnl")
    hydro_path = os.path.join(path, "Hydrostatic.in")
    m = case.multiple
    if case.hull == "sharp":
        params = scaled(fixtures.DRIFT_CYL_MESH, ("ntheta", "nz", "nr"), m)
        if case.edge_power is not None:
            params.update(graded=False, edge_power=case.edge_power)
        if case.damping:
            params.update(external_damping={j - 1: v for j, v in case.damping.items()})
        mesh = make_mesh.write_cylinder(mesh_path, hydro_path, **params)
    elif case.hull == "round":
        mesh = make_mesh.write_rounded_cylinder(mesh_path, hydro_path,
                                                **scaled(fixtures.ROUND_MESH, ("ntheta", "nz", "narc", "nr"), m))
    elif case.hull == "bottom":
        mesh = make_mesh.write_cylinder(mesh_path, hydro_path,
                                        **scaled(dict(fixtures.BOTTOM_MOUNTED_MESH, ntheta=32, nz=8), ("ntheta", "nz"), m))
    elif case.hull == "box":
        mesh = make_mesh.write_box(mesh_path, hydro_path, **scaled(fixtures.DRIFT_BOX_MESH, ("nx", "ny", "nz"), m))
    elif case.hull == "round_box":
        mesh = make_mesh.write_rounded_box(mesh_path, hydro_path,
                                           **scaled(ROUND_BOX_MESH, ("nx", "ny", "ncorner", "nz", "narc"), m))
    elif case.hull == "corner_box":
        mesh = make_mesh.write_rounded_box(mesh_path, hydro_path,
                                           **scaled(CORNER_BOX_MESH, ("nx", "ny", "ncorner", "nz", "narc"), m))
    elif case.hull in ("oc4", "oc4w"):
        # HullMesh.pnl, Hydrostatic.in (mooring as external restoring), CSMesh_1.pnl and FSMesh_1.pnl, checked; oc4w has the
        # matrices of the WAMIT run
        damping = {j - 1: v for j, v in case.damping.items()} if case.damping else None
        written = make_mesh.write_oc4(path, m, external_damping=damping, model="wamit" if case.hull == "oc4w" else "platform")
        mesh = written["hull"]
        points = len(written["control"].panels) + len(written["free"].panels)
        if points != HULLS[case.hull]["cf_points"][m]:
            raise ValueError("%s has %d panels on C and F, not the %d of the job table" % (
                case.name, points, HULLS[case.hull]["cf_points"][m]))
    else:
        raise ValueError("unknown hull " + case.hull)
    if len(mesh.panels) != panels_of(case):
        raise ValueError("%s has %d panels, not the %d of the job table" % (case.name, len(mesh.panels), panels_of(case)))
    if case.irsp:
        lid = make_mesh.write_waterplane(os.path.join(path, "WaterPlaneMesh.pnl"), 10.0, int(round(32 * m)), int(round(6 * m)))
        if len(lid.panels) != lid_panels_of(case):
            raise ValueError("the lid of %s has %d panels" % (case.name, len(lid.panels)))
    return mesh


def write_hams_case(case, root, threads):
    path = os.path.join(root, case.name, "input")
    os.makedirs(path, exist_ok=True)
    write_hams_hull(case, path)
    fixtures.write_control(os.path.join(path, "ControlFile.in"), depth=HULLS[case.hull]["depth"], frequencies=case.frequencies,
                           headings=case.headings, solver=2 if case.isolv == 2 else None, irsp=case.irsp, threads=threads)
    options = dict(mean_drift=case.mean_drift, drift_methods="far near", motion_state=case.motion,
                   write_terms=1 if case.mean_drift == 1 else 0, velocity_method=case.velocity,
                   waterline_method=case.waterline)
    if case.control:
        options.update(drift_methods="far near control", cs_velocity=case.control)
    if case.qtf:
        options.update(qtf=case.qtf)
    if case.qtf_terms:
        options.update(qtf_terms=case.qtf_terms)
    if case.qtf_max_dw:
        options.update(qtf_max_dw=case.qtf_max_dw)
    if case.qtf_edge_weights:
        options.update(qtf_edge_weights=case.qtf_edge_weights)
    options.update(case.extra_options)
    fixtures.write_second_order(os.path.join(path, "SecondOrder.in"), **options)
    return path


def read_matrices(path):
    """The matrices of a Hydrostatic.in file, by the order they are written in."""
    matrices = []
    current = None
    for line in open(path):
        try:
            row = [float(text) for text in line.split()]
        except ValueError:
            current = None
            continue
        if len(row) == 6:
            if current is None:
                current = []
                matrices.append(current)
            current.append(row)
    return {"mass": matrices[0], "restoring": matrices[3], "external": matrices[4]}


def nemoh_matrix(path, matrix):
    with open(path, "w") as f:
        for row in matrix:
            f.write(" ".join("%.7E" % value for value in row) + "\n")


def write_nemoh_case(case, root):
    """Write the input of a Nemoh run. The nodes and panels are those of the HAMS case. No lid and no external loads, and the
    matrices are those given to HAMS, about the same point: the external restoring of Hydrostatic.in (the OC4 mooring) is
    Km.dat of a free body. A fixed body has a large external stiffness instead (HULLS fixed_stiffness, else
    NEMOH_FIXED_STIFFNESS)."""
    base = os.path.join(root, case.name)
    for folder in ("Mechanics", "Motion", os.path.join("results", "QTF"), os.path.join("results", "sources"), "QTFPreprocOut",
                   "mesh", "hams_hull"):
        os.makedirs(os.path.join(base, folder), exist_ok=True)
    mesh = write_hams_hull(case, os.path.join(base, "hams_hull"))
    make_mesh.write_nemoh_mesh(os.path.join(base, "hull.dat"), mesh)
    matrices = read_matrices(os.path.join(base, "hams_hull", "Hydrostatic.in"))
    zero = [[0.0] * 6 for _ in range(6)]
    if case.motion == "fixed":
        fixed = HULLS[case.hull].get("fixed_stiffness", NEMOH_FIXED_STIFFNESS)
        stiffness = [[fixed if i == j else 0.0 for j in range(6)] for i in range(6)]
    else:
        stiffness = matrices["external"]
    nemoh_matrix(os.path.join(base, "Mechanics", "Inertia.dat"), matrices["mass"])
    nemoh_matrix(os.path.join(base, "Mechanics", "Kh.dat"), matrices["restoring"])
    nemoh_matrix(os.path.join(base, "Mechanics", "Km.dat"), stiffness)
    nemoh_matrix(os.path.join(base, "Mechanics", "Badd.dat"), zero)
    # One range for the first-order and the QTF frequencies: the QTF part needs first-order solutions at every difference
    # frequency, down to the grid step (the Nemoh user manual asks for first-order coverage of [dw, w_max])
    low, high, count = case.nemoh_range
    modes = ("1 1. 0. 0. 0. 0. 0.", "1 0. 1. 0. 0. 0. 0.", "1 0. 0. 1. 0. 0. 0.",
             "2 1. 0. 0. 0. 0. 0.", "2 0. 1. 0. 0. 0. 0.", "2 0. 0. 1. 0. 0. 0.")
    lines = ["--- Environment ---", "%.1f" % make_mesh.RHO, "%.5f" % make_mesh.G, "%.1f" % HULLS[case.hull]["depth"], "0. 0.",
             "--- Description of floating bodies ---", "1", "--- Body 1 ---", "hull.dat",
             "%d %d" % (len(mesh.nodes), len(mesh.panels)), "6"]
    lines += list(modes) + ["6"] + list(modes) + ["0", "--- Load cases ---", "1 %d %.4f %.4f" % (count, low, high),
                                                  "%d %.4f %.4f" % (len(case.headings), case.headings[0], case.headings[-1]),
                                                  "--- Post processing ---", "0 0.1 10.", "0", "0 0. 180.", "0 50 400. 400.",
                                                  "1", "1", "--- QTF ---", "1", "%d %.4f %.4f" % (count, low, high),
                                                  "%d" % (1 if case.mean_drift == 2 else 0), "%d" % case.nemoh_contrib,
                                                  "NA", "0 0 0", "1", "1"]
    lines += ["%d" % flag for flag in NEMOH_QTF_TERMS[case.nemoh_contrib]] + ["---"]
    with open(os.path.join(base, "Nemoh.cal"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(base, "input_solver.txt"), "w") as f:
        f.write("2\n0.001\n1\n10 1e-5 1000\n")
    return base


def peak_memory_bytes(process):
    """Peak memory of a finished child process, or None where it cannot be read."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]

        counters = Counters()
        counters.cb = ctypes.sizeof(Counters)
        psapi = ctypes.WinDLL("psapi")
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if psapi.GetProcessMemoryInfo(int(process._handle), ctypes.byref(counters), counters.cb):
            return float(counters.PeakWorkingSetSize)
        return None
    import resource
    # Of all the children so far. Each Slurm task runs one solver, so this is that solver
    return 1024.0 * resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss


def run_program(command, log_path, cwd=None):
    """Run a program with its output in a log. Returns the return code, wall seconds and peak memory in bytes."""
    start = time.time()
    with open(log_path, "w") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
        code = process.wait()
    return code, time.time() - start, peak_memory_bytes(process)


def run_hams(case, root, hams, threads):
    path = write_hams_case(case, root, threads)
    base = os.path.join(root, case.name)
    code, wall, peak = run_program([hams, path + os.sep, os.path.join(base, "output") + os.sep], os.path.join(base, "log.txt"))
    # The solver stops on a bad input with return code 0, so its last message is the check
    finished = "completes successfully" in open(os.path.join(base, "log.txt"), errors="replace").read()
    return dict(code=code if finished else (code or 1), wall_s=wall, peak_bytes=peak)


def run_nemoh(case, root, build, threads):
    base = write_nemoh_case(case, root)
    stages = {}
    peak = 0.0
    code = 0
    for folder, program in NEMOH_STAGES:
        command = [os.path.abspath(os.path.join(build, folder, program + (".exe" if sys.platform == "win32" else ""))), "."]
        code, wall, memory = run_program(command, os.path.join(base, "log_%s.txt" % program), cwd=base)
        stages[program] = wall
        peak = max(peak, memory or 0.0)
        if code:
            break
    return dict(code=code, wall_s=sum(stages.values()), peak_bytes=peak, stages=stages)


def binary_for(case, args):
    return args.hams_grn8 if case.grn == 8 else args.hams


def run_case(case, args):
    if case.solver == "nemoh":
        if not args.nemoh_build:
            raise SystemExit("%s needs --nemoh-build" % case.name)
        result = run_nemoh(case, args.root, args.nemoh_build, args.threads)
    else:
        binary = binary_for(case, args)
        if not binary:
            raise SystemExit("%s needs --hams-grn8, a build with GRN_KIND = 8" % case.name)
        result = run_hams(case, args.root, os.path.abspath(binary), args.threads)
    result.update(name=case.name, panels=panels_of(case), lid=lid_panels_of(case), threads=args.threads,
                  host=socket.gethostname(), estimate_hours=wall_hours(case), estimate_gb=memory_gb(case))
    with open(os.path.join(args.root, case.name, "result.json"), "w") as f:
        json.dump(result, f, indent=1)
    print("%s: code %d, %.1f s, peak %.3f GB" % (case.name, result["code"], result["wall_s"], (result["peak_bytes"] or 0.0) / 1e9))
    return result


def calibrate(args):
    """Time small runs on this machine. Prints the fit of time per frequency to A N^2 + B N^3."""
    rows = []
    for multiple in (1, 1.5, 2):
        for velocity in ("differencing", "source"):
            case = Case("cal_%s_%g" % (velocity, multiple), "cal", "", "sharp", multiple, "free", velocity,
                        frequencies=(0.95,))
            rows.append((panels_of(case), velocity, run_case(case, args)))
            log = open(os.path.join(args.root, case.name, "log.txt"), errors="replace").read()
            found = re.search(r"source solve:\s+([\d.]+) s to assemble and factorize", log)
            print("   source assembly and factorization: %s s" % (found.group(1) if found else "not logged"))
    for velocity in ("differencing", "source"):
        points = [(n, result["wall_s"]) for n, v, result in rows if v == velocity]
        a, b = fit_n2_n3(points)
        print("%s: %.4g N^2 + %.4g N^3 s per frequency, residuals %s" % (
            velocity, a, b, ["%+.1f%%" % (100.0 * (a * n ** 2 + b * n ** 3 - t) / t) for n, t in points]))


def fit_n2_n3(points):
    """Least squares of t/n^2 = a + b n, which weighs the sizes evenly in relative terms."""
    xs = [float(n) for n, _ in points]
    ys = [t / n ** 2 for n, t in points]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return my - b * mx, b


def read_drift(path):
    """Rows of a mean drift file: real part by (omega, beta1, beta2, index)."""
    values = {}
    for line in open(path):
        columns = line.split()
        if columns and not line.startswith("#"):
            values[(round(float(columns[0]), 4), round(float(columns[1]), 3), round(float(columns[2]), 3),
                    int(columns[3]))] = float(columns[6])
    return values


def read_terms(path):
    """Terms of F_x, by (omega, beta1, beta2), in the order of the file."""
    values = {}
    if not os.path.exists(path):
        return values
    for line in open(path):
        if line.startswith("#"):
            continue
        columns = line.split()
        if int(columns[4]) == 1:
            key = (round(float(columns[0]), 4), round(float(columns[1]), 3), round(float(columns[2]), 3))
            values.setdefault(key, []).append(float(columns[6]))
    return values


def wave_number(omega, depth):
    """Solve omega^2 = g k tanh(k h) by Newton from the deep water value."""
    nu = omega ** 2 / make_mesh.G
    k = max(nu, math.sqrt(nu / depth))
    for _ in range(100):
        step = (k * math.tanh(k * depth) - nu) / (math.tanh(k * depth) + k * depth / math.cosh(k * depth) ** 2)
        k -= step
        if abs(step) < 1e-14 * k:
            break
    return k


def maccamy_fuchs_fx(omega, radius, depth, count=720):
    """Mean surge drift on a bottom-mounted cylinder, over rho g, per unit amplitude squared, from the far-field integral of
    the exact scattered amplitude. The same formula as the unit tests."""
    from scipy.special import h1vp, jvp

    k = wave_number(omega, depth)
    x = 2.0 * k * depth
    ratio = 0.5 * (1.0 + 2.0 * x * math.exp(-x) / (1.0 - math.exp(-2.0 * x)))
    orders = int(k * radius) + 20
    coefficient = [(1.0 if n == 0 else 2.0) * jvp(n, k * radius) / h1vp(n, k * radius) for n in range(orders + 1)]
    total = 0.0
    for j in range(count):
        theta = 2.0 * math.pi * j / count
        amplitude = sum(c * math.cos(n * theta) for n, c in enumerate(coefficient))
        total += abs(amplitude) ** 2 * 2.0 / (math.pi * k) * (1.0 - math.cos(theta))
    return 0.5 * ratio * total * 2.0 * math.pi / count


def percent(value, reference, scale):
    """Relative error in percent, or an empty field where the reference is too small to divide by."""
    if abs(reference) < 1e-4 * scale:
        return ""
    return "%+.3f" % (100.0 * (value - reference) / abs(reference))


# Tolerance of the control-surface gates: 1% of the reference, with a floor for loads near zero. The .7, .8 and .9 values are per
# rho g A^2 L^k with L = 1 m, so the floor is 0.001 a^k, a = 10 m the radius of the cylinders, k = 1 for forces and 2 for moments
GATE_RELATIVE = 0.01
GATE_FLOOR = 0.001
GATE_LENGTH = 10.0
VERTICAL = (("Fz", 3), ("Mx", 4), ("My", 5))


def gate_ratio(value, reference, index):
    """|value - reference| over the tolerance of the gates. Above 1 fails."""
    power = 1 if index <= 3 else 2
    return abs(value - reference) / max(GATE_RELATIVE * abs(reference), GATE_FLOOR * GATE_LENGTH ** power)


def m4_reference_name(case):
    """The near-field reference of a green-route case: the m4 differencing case of the same cylinder and motion in group A."""
    if case.hull not in ("round", "sharp") or case.velocity != "differencing":
        return None
    return "a_%s_%s_diff_m4" % (case.hull, case.motion)


def read_m4_near(case, roots):
    """Near-field rows of the m4 reference, from the first of the roots that holds them, or {}."""
    name = m4_reference_name(case)
    for root in roots:
        if name and root:
            path = os.path.join(root, name, "output", "Wamit_format", "MeanDriftPressure.9")
            if os.path.exists(path):
                return read_drift(path)
    return {}


def collect_hams(case, root, result, reference_roots=()):
    folder = os.path.join(root, case.name, "output", "Wamit_format")
    far = read_drift(os.path.join(folder, "MeanDriftMomentum.8"))
    near = read_drift(os.path.join(folder, "MeanDriftPressure.9"))
    terms = read_terms(os.path.join(folder, "MeanDriftTerms.txt"))
    # Same format as MeanDriftPressure.9: rows I = 1 to 6 of the control-surface method, or I = 1, 2, 6 alone where the
    # vertical loads are not computed (C reaches the sea bed, or F could not be built). Rows are keyed by I, so either works
    control_path = os.path.join(folder, "MeanDriftControlSurface.7")
    control = read_drift(control_path) if case.control and os.path.exists(control_path) else {}
    near_m4 = read_m4_near(case, (root,) + tuple(reference_roots)) if control else {}
    route = "green" if case.control == "green" else "src"
    rows = []
    for omega, beta1, beta2 in sorted(set(key[:3] for key in far)):
        if beta1 != beta2 and case.mean_drift == 1:
            continue
        get = lambda table, index: table.get((omega, beta1, beta2, index), float("nan"))
        scale = math.hypot(get(far, 1), get(far, 2))
        row = dict(case=case.name, group=case.group, hull=case.hull, motion=case.motion, velocity=case.velocity,
                   waterline=case.waterline, grn=case.grn, isolv=case.isolv, irsp=case.irsp, panels=panels_of(case), omega=omega,
                   beta1=beta1, beta2=beta2, wall_s="%.1f" % result["wall_s"],
                   peak_gb="%.3f" % ((result.get("peak_bytes") or 0.0) / 1e9))
        for label, index in (("Fx", 1), ("Fy", 2), ("Mz", 6)):
            row[label + "_near"] = "%.6g" % get(near, index)
            row[label + "_far"] = "%.6g" % get(far, index)
            row[label + "_err_pct"] = percent(get(near, index), get(far, index), scale if index < 3 else 10.0 * scale)
            if control:
                column = "%s_ctrl_%s" % (label, route)
                row[column] = "%.6g" % get(control, index)
                row[column + "_err_pct"] = percent(get(control, index), get(far, index), scale if index < 3 else 10.0 * scale)
        if case.hull == "bottom" and beta1 == 0.0:
            exact = maccamy_fuchs_fx(omega, BOTTOM_RADIUS, BOTTOM_DEPTH)
            row["Fx_exact"] = "%.6g" % exact
            row["Fx_near_exact_pct"] = percent(get(near, 1), exact, exact)
            row["Fx_far_exact_pct"] = percent(get(far, 1), exact, exact)
            if control:
                row["Fx_ctrl_%s_exact_pct" % route] = percent(get(control, 1), exact, exact)
        # Phase B: the vertical loads of the control-surface method, against near-field of this run and at m4
        for label, index in VERTICAL:
            key = (omega, beta1, beta2, index)
            if key not in control:
                continue
            column = "%s_ctrl_%s" % (label, route)
            scale_v = scale if index == 3 else 10.0 * scale
            row[label + "_near"] = "%.6g" % get(near, index)
            row[column] = "%.6g" % control[key]
            row[column + "_near_pct"] = percent(control[key], get(near, index), scale_v)
            row[column + "_near_gate"] = "%.4f" % gate_ratio(control[key], get(near, index), index)
            if key in near_m4:
                row[label + "_near_m4"] = "%.6g" % near_m4[key]
                row[label + "_near_m4_pct"] = percent(get(near, index), near_m4[key], scale_v)
                row[column + "_m4_pct"] = percent(control[key], near_m4[key], scale_v)
                row[column + "_m4_gate"] = "%.4f" % gate_ratio(control[key], near_m4[key], index)
        for number, value in enumerate(terms.get((omega, beta1, beta2), []), start=1):
            row["Fx_term%d" % number] = "%.5f" % value
        rows.append(row)
    return rows


def read_nemoh_drift(path):
    """The diagonal of OUT_QTFM_N.dat: real part over rho g, by (omega, beta1, beta2, dof)."""
    values = {}
    for line in open(path):
        columns = line.split()
        if len(columns) < 9 or not columns[0][0].isdigit():
            continue
        if abs(float(columns[0]) - float(columns[1])) < 1e-6:
            values[(round(float(columns[0]), 2), round(float(columns[2]), 1), round(float(columns[3]), 1),
                    int(columns[4]))] = float(columns[7])
    return values


def read_nemoh_matrix(path):
    """Added mass or damping of a Nemoh output file, by rounded frequency, as a 6 x 6 list."""
    numbers = open(path).read().split("\n")
    out = {}
    i = 1
    while i + 6 < len(numbers) and numbers[i].strip():
        out[round(float(numbers[i]), 2)] = [[float(text) for text in numbers[i + 1 + r].split()] for r in range(6)]
        i += 7
    return out


def read_nemoh_tec(path, zone=1):
    """Numeric rows of one zone of a Tecplot output file of Nemoh, by rounded frequency. zone counts from 1."""
    rows = {}
    current = 0
    for line in open(path, errors="replace"):
        if line.lstrip().lower().startswith("zone"):
            if "I=" in line:
                current += 1
            continue
        if current != zone:
            continue
        try:
            values = [float(text) for text in line.split()]
        except ValueError:
            continue
        if values:
            rows[round(values[0], 2)] = values[1:]
    return rows


def read_hams_first_order(folder):
    """A, B, |X| and |RAO| of the WAMIT format output."""
    a, b, x, rao = {}, {}, {}, {}
    for line in open(os.path.join(folder, "AmssDamp.1")):
        c = line.split()
        if len(c) >= 5:
            a[(round(float(c[0]), 2), int(c[1]), int(c[2]))] = float(c[3])
            b[(round(float(c[0]), 2), int(c[1]), int(c[2]))] = float(c[4])
    for name, table in (("ExcForce.3", x), ("Motion.4", rao)):
        for line in open(os.path.join(folder, name)):
            c = line.split()
            if len(c) >= 5:
                table[(round(float(c[0]), 2), round(float(c[1]), 1), int(c[2]))] = float(c[3])
    return a, b, x, rao


def hams_case_for_nemoh(case, cases, velocity):
    """The HAMS case on the same mesh and motion as a Nemoh case, or None."""
    short = {"differencing": "diff", "source": "src"}[velocity]
    if case.mean_drift == 2:
        name = "n_pair_%s_%s_%s" % (case.hull, case.motion, tag(case.multiple))
        return cases.get(name) if velocity == "differencing" else None
    for name in ("a_%s_%s_%s_%s" % (case.hull, case.motion, short, tag(case.multiple)),
                 "n_%s_%s_%s_%s" % (case.hull, case.motion, short, tag(case.multiple)),
                 "d_base_%s_%s_%s_%s" % (case.hull, case.motion, short, tag(case.multiple))):
        if name in cases:
            return cases[name]
    return None


def collect_nemoh(case, root, cases, result):
    """Rows of Nemoh against the HAMS case on the same mesh, for the shared frequencies. A second list has first-order values."""
    base = os.path.join(root, case.name)
    drift = read_nemoh_drift(os.path.join(base, "results", "QTF", "OUT_QTFM_N.dat"))
    rows = []
    first_order = []
    for velocity in ("differencing", "source"):
        partner = hams_case_for_nemoh(case, cases, velocity)
        if partner is None:
            continue
        folder = os.path.join(root, partner.name, "output", "Wamit_format")
        near_path = os.path.join(folder, "MeanDriftPressure.9")
        if not os.path.exists(near_path):
            continue
        near = read_drift(near_path)
        far = read_drift(os.path.join(folder, "MeanDriftMomentum.8"))
        for (omega, beta1, beta2, index), value in sorted(near.items()):
            key = (round(omega, 2), round(beta1, 1), round(beta2, 1), index)
            if index in (1, 2, 6) and key in drift:
                scale = math.hypot(near.get((omega, beta1, beta2, 1), 0.0), near.get((omega, beta1, beta2, 2), 0.0))
                rows.append(dict(nemoh_case=case.name, hams_case=partner.name, velocity=velocity, panels=panels_of(case), omega=omega,
                                 beta1=beta1, beta2=beta2, index=index, hams_near="%.6g" % value,
                                 hams_far="%.6g" % far.get((omega, beta1, beta2, index), float("nan")),
                                 nemoh="%.6g" % drift[key], near_against_nemoh_pct=percent(value, drift[key], scale)))
    partner = hams_case_for_nemoh(case, cases, "differencing")
    if partner is not None and os.path.exists(os.path.join(root, partner.name, "output", "Wamit_format", "AmssDamp.1")):
        a, b, x, rao = read_hams_first_order(os.path.join(root, partner.name, "output", "Wamit_format"))
        added = read_nemoh_matrix(os.path.join(base, "results", "CM.dat"))
        damping = read_nemoh_matrix(os.path.join(base, "results", "CA.dat"))
        excitation = read_nemoh_tec(os.path.join(base, "results", "ExcitationForce.tec"))
        motion = read_nemoh_tec(os.path.join(base, "Motion", "RAO.tec"))
        for omega in sorted(set(key[0] for key in a)):
            if omega not in added:
                continue
            for i, j in ((1, 1), (3, 3), (5, 5), (1, 5)):
                first_order.append(dict(quantity="A%d%d" % (i, j), panels=panels_of(case), omega=omega,
                                        hams=a[(omega, i, j)] * make_mesh.RHO, nemoh=added[omega][i - 1][j - 1]))
                first_order.append(dict(quantity="B%d%d" % (i, j), panels=panels_of(case), omega=omega,
                                        hams=b[(omega, i, j)] * make_mesh.RHO * omega, nemoh=damping[omega][i - 1][j - 1]))
            for i in (1, 3, 5):
                if (omega, 0.0, i) in x and omega in excitation:
                    first_order.append(dict(quantity="X%d" % i, panels=panels_of(case), omega=omega,
                                            hams=x[(omega, 0.0, i)] * make_mesh.RHO * make_mesh.G,
                                            nemoh=excitation[omega][2 * (i - 1)]))
                if (omega, 0.0, i) in rao and omega in motion and case.motion == "free":
                    # Nemoh gives rotations in degrees per metre
                    first_order.append(dict(quantity="RAO%d" % i, panels=panels_of(case), omega=omega,
                                            hams=rao[(omega, 0.0, i)],
                                            nemoh=motion[omega][i - 1] * (math.pi / 180.0 if i > 3 else 1.0)))
        for row in first_order:
            row["nemoh_case"] = case.name
            row["difference_pct"] = percent(row["hams"], row["nemoh"], 1e-12)
    return rows, first_order


def write_csv(path, rows):
    columns = []
    for row in rows:
        columns += [key for key in row if key not in columns]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def gate_lines(cases, rows):
    """Lines of the control-surface gates: per case and load, the worst error over the tolerance, with where it is.
    Horizontal loads against far-field (Phase A), vertical loads against near-field of the run and at m4 (Phase B)."""
    worst = {}
    for row in rows:
        case = cases[row["case"]]
        if not case.control:
            continue
        route = "green" if case.control == "green" else "src"
        checks = [(label, index, "far", label + "_far") for label, index in (("Fx", 1), ("Fy", 2), ("Mz", 6))]
        checks += [(label, index, ref, label + ("_near" if ref == "near" else "_near_m4"))
                   for label, index in VERTICAL for ref in ("near", "near_m4")]
        for label, index, ref, ref_column in checks:
            column = "%s_ctrl_%s" % (label, route)
            if column not in row or ref_column not in row:
                continue
            value, reference = float(row[column]), float(row[ref_column])
            if math.isnan(value) or math.isnan(reference):
                continue
            ratio = gate_ratio(value, reference, index)
            key = (row["case"], label, ref)
            count, best = worst.get(key, (0, -1.0))[:2]
            if ratio > best:
                worst[key] = (count + 1, ratio, row["omega"], row["beta1"], reference, value)
            else:
                worst[key] = (count + 1,) + worst[key][1:]
    lines = ["Case load reference | n | worst err/tol | omega beta reference control"]
    for (name, label, ref), (count, ratio, omega, beta, reference, value) in sorted(worst.items()):
        lines.append("%s %s %s | %d | %.3f | %g %g %.6g %.6g %s" % (name, label, ref, count, ratio, omega, beta, reference, value,
                                                                    "FAIL" if ratio > 1.0 else ""))
    # Phase B gate: the rounded cylinder at m3, green route, F_z, M_x, M_y within the tolerance of near-field of the same run
    gate = [(key, item) for key, item in worst.items() if key[2] == "near" and cases[key[0]].hull == "round"
            and cases[key[0]].multiple == 3 and cases[key[0]].control == "green"]
    lines.append("")
    lines.append("Phase B gate (F_z, M_x, M_y of control against near-field, rounded cylinder, m3, green route):")
    if not gate:
        lines.append("  no case: not checked")
    for name in sorted(set(key[0] for key, _ in gate)):
        ratio = max(item[1] for key, item in gate if key[0] == name)
        lines.append("  %s: worst %.3f of the tolerance: %s" % (name, ratio, "pass" if ratio <= 1.0 else "FAIL"))
    return lines


def collect(args):
    """Write the table of drift values and errors for each finished case, and against Nemoh. Prints time and memory against the
    estimates, and with --control the lines of the control-surface gates."""
    cases = {case.name: case for case in cases_of(args.control, study_of(args))}
    hams_rows, nemoh_rows, first_rows = [], [], []
    print("| Case | Code | Panels | Wall (s) | Estimate (s) | Peak (GB) | Estimate (GB) |")
    print("| :--- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for name in sorted(os.listdir(args.root)):
        path = os.path.join(args.root, name, "result.json")
        if not os.path.exists(path):
            continue
        result = json.load(open(path))
        print("| %s | %d | %d | %.1f | %.0f | %.3f | %.3f |" % (
            name, result["code"], result["panels"], result["wall_s"], result["estimate_hours"] * 3600.0,
            (result.get("peak_bytes") or 0.0) / 1e9, result["estimate_gb"]))
        if result["code"] or name not in cases:
            continue
        if cases[name].solver == "nemoh":
            rows, first = collect_nemoh(cases[name], args.root, cases, result)
            nemoh_rows += rows
            first_rows += first
        else:
            hams_rows += collect_hams(cases[name], args.root, result, (args.reference,) if args.reference else ())
    out = args.out or os.path.join(args.root, "drift_table.csv")
    write_csv(out, hams_rows)
    print("%d rows in %s" % (len(hams_rows), out))
    for label, rows in (("nemoh", nemoh_rows), ("first_order", first_rows)):
        if rows:
            other = out.replace(".csv", "_%s.csv" % label)
            write_csv(other, rows)
            print("%d rows in %s" % (len(rows), other))
    if args.control:
        lines = gate_lines(cases, hams_rows)
        gate_path = out.replace(".csv", "_gate.txt")
        with open(gate_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print("\n".join(lines))
        print("Gate lines in %s" % gate_path)


def fortran_files_in_order(folders):
    """Fortran files of the folders, ordered so that every module is compiled before the files that use it. Files holding a
    program are left out."""
    files, provides, needs = [], {}, {}
    for folder in folders:
        for name in sorted(os.listdir(folder)):
            if not name.lower().endswith((".f90",)):
                continue
            path = os.path.join(folder, name)
            text = open(path, errors="replace").read().lower()
            if re.search(r"^\s*program\s+\w+", text, re.M):
                continue
            files.append(name)
            for module in re.findall(r"^\s*module\s+(?!procedure)(\w+)", text, re.M):
                provides[module] = name
            needs[name] = set(re.findall(r"^\s*use\s*(?:,\s*\w+\s*::\s*)?(\w+)", text, re.M))
    ordered, done = [], set()

    def visit(name, stack=()):
        if name in done:
            return
        if name in stack:
            raise ValueError("files use each other: " + name)
        for module in sorted(needs[name]):
            if module in provides and provides[module] != name:
                visit(provides[module], stack + (name,))
        done.add(name)
        ordered.append(name)

    for name in files:
        visit(name)
    return ordered


SBATCH_HEAD = """#!/bin/bash
#SBATCH -J {name}
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task {cpus}
#SBATCH --mem={memory}G
#SBATCH -t {limit}
#SBATCH -o {logs}/{log}
{extra}module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS={cpus}
export OMP_STACKSIZE=1G
"""


def time_limit(hours):
    """A Slurm time limit of twice the estimate, at least 30 minutes."""
    minutes = int(max(30.0, 120.0 * hours))
    return "%d:%02d:00" % (minutes // 60, minutes % 60)


def write_tests_script(args, out, root, logs, src, head):
    with open(os.path.join(out, "tests.sbatch"), "w", newline="\n") as f:
        f.write(SBATCH_HEAD.format(name="drift_tests", limit="1:00:00", log="tests_%j.log", **head))
        f.write("# Unit and integration tests under ifx on Linux, for the thin margins of the fixtures and the suite time\n")
        f.write("# TESTDRIVE is the folder of the test-drive source files, testdrive.F90 and testdrive_version.f90\n")
        f.write("TESTDRIVE=${TESTDRIVE:?set TESTDRIVE}\nFLAGS=\"%s\"\n" % IFX_FLAGS)
        test_folder = os.path.join(os.path.dirname(src), "test")
        files = [name for name in fortran_files_in_order([src, test_folder])]
        # The testers find test/data three folders above themselves, as under fpm
        f.write("d=%s/build/ifx/test; rm -rf $d; mkdir -p $d; cd $d\n" % args.repo)
        f.write("cp %s/src/*.f90 %s/test/*.f90 $TESTDRIVE/testdrive_version.f90 $TESTDRIVE/testdrive.F90 .\n" % (args.repo, args.repo))
        f.write("OBJ=\"\"\nfor file in testdrive_version.f90 testdrive.F90 %s; do\n" % " ".join(files))
        f.write("  ifx $FLAGS -fpp -c $file || exit 1; OBJ=\"$OBJ ${file%.*}.o\"\ndone\n")
        f.write("ifx $FLAGS -o unit-tester unit_tester.f90 $OBJ\nifx $FLAGS -o integration-tester integration_tester.f90 $OBJ\n")
        f.write("cd %s\n" % args.repo)
        f.write("time timeout 600 $d/unit-tester > %s/unit.log 2>&1; echo unit rc=$?\n" % logs)
        f.write("time timeout 1800 $d/integration-tester %s/hamsmrel_grn4 > %s/integration.log 2>&1; echo integration rc=$?\n" % (root, logs))


def write_slurm(args):
    """Write the scripts of the study into <root>/jobs. Nothing is submitted."""
    root = args.root.rstrip("/")
    jobs = "%s/jobs" % root
    logs = "%s/logs" % root
    out = args.out or jobs
    os.makedirs(out, exist_ok=True)
    script = "%s/pre_postprocessing_python/SecondOrder/drift_convergence.py" % args.repo
    extra = "#SBATCH -A %s\n" % args.account if args.account else ""
    cases = [case for case in cases_of(args.control, study_of(args)) if in_tier(case, args.tier) and case.where == "snellius"]
    head = dict(partition=PARTITION, cpus=CPUS, memory=MEMORY_GB, logs=logs, extra=extra)

    src = os.path.join(args.source or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."), "src")
    order = " ".join(fortran_files_in_order([src]))
    with open(os.path.join(out, "build.sbatch"), "w", newline="\n") as f:
        f.write(SBATCH_HEAD.format(name="drift_build", limit="0:30:00", log="build_%j.log", **head))
        f.write("set -e\nFLAGS=\"%s\"\nLIST=\"%s HAMS_Prog.f90\"\n" % (IFX_FLAGS, order))
        f.write("# GRN_KIND is a compile-time constant, so the second build changes it in a copy of the source\n")
        f.write("for kind in 4 8; do\n  d=%s/build_grn$kind; rm -rf $d; mkdir -p $d; cp %s/src/*.f90 $d/; cd $d\n" % (root, args.repo))
        f.write("  sed -i \"s/GRN_KIND = 4/GRN_KIND = $kind/\" WavDynMods.f90\n")
        f.write("  grep -q \"GRN_KIND = $kind\" WavDynMods.f90\n")
        f.write("  ifx $FLAGS -o %s/hamsmrel_grn$kind $LIST\ndone\n" % root)

    if not args.control:   # the control study needs the build only: its tests run on the laptop
        write_tests_script(args, out, root, logs, src, head)

    # Arrays by the size of the run, since an array has one time limit
    classes = (("short", 0.0, 0.5), ("medium", 0.5, 2.0), ("long", 2.0, 1e9))
    for label, low, high in classes:
        chosen = [case for case in cases if case.solver == "hams" and low <= wall_hours(case) < high]
        if not chosen:
            continue
        with open(os.path.join(out, label + ".txt"), "w", newline="\n") as f:
            f.write("\n".join(case.name for case in chosen) + "\n")
        limit = time_limit(max(wall_hours(case) for case in chosen))
        with open(os.path.join(out, label + ".sbatch"), "w", newline="\n") as f:
            f.write(SBATCH_HEAD.format(name="drift_" + label, limit=limit, log="%A_%a.log",
                                       **dict(head, extra="#SBATCH --array=1-%d%%8\n%s" % (len(chosen), extra))))
            f.write("CASE=$(sed -n \"${SLURM_ARRAY_TASK_ID}p\" %s/%s.txt)\n" % (jobs, label))
            f.write("python3 %s run %s --case \"$CASE\" --hams %s/hamsmrel_grn4 --hams-grn8 %s/hamsmrel_grn8 --threads %d%s\n" % (
                script, root, root, root, CPUS,
                " --control " + args.control if args.control else qtf_flag(study_of(args))))

    # Nemoh is run with one thread per case, as many at a time as there are cases, in one allocation per batch
    # (nemoh_batches: nemoh.txt, or nemoh_<hull>_m<k>.txt for a hull with batches of its own)
    for label, chosen in nemoh_batches([case for case in cases if case.solver == "nemoh"]):
        if len(chosen) > NEMOH_CPUS:
            raise ValueError("%s has %d Nemoh runs for %d CPUs" % (label, len(chosen), NEMOH_CPUS))
        with open(os.path.join(out, label + ".txt"), "w", newline="\n") as f:
            f.write("\n".join(case.name for case in chosen) + "\n")
        limit = time_limit(max(wall_hours(case) for case in chosen))
        with open(os.path.join(out, label + ".sbatch"), "w", newline="\n") as f:
            f.write(SBATCH_HEAD.format(name="drift_" + label, limit=limit, log=label + "_%j.log",
                                       **dict(head, cpus=NEMOH_CPUS, memory=NEMOH_MEMORY_GB)).replace(
                "export OMP_NUM_THREADS=%d" % NEMOH_CPUS, "export OMP_NUM_THREADS=1"))
            f.write("# Needs a build of Nemoh on Snellius, kept outside the repository. Set NEMOH_BUILD to its build folder.\n")
            f.write("NEMOH_BUILD=${NEMOH_BUILD:?set NEMOH_BUILD}\n")
            f.write("xargs -P %d -I{} python3 %s run %s --case {} --hams none --nemoh-build $NEMOH_BUILD --threads 1%s < %s/%s.txt\n" % (
                len(chosen), script, root, qtf_flag(study_of(args)), jobs, label))
    print("Scripts in %s, for %s on the cluster, with logs in %s. Nothing was submitted." % (out, jobs, logs))


def nemoh_batches(cases):
    """[(label, cases)] of the Nemoh allocations: one "nemoh" for the hulls without batches of their own, and one
    "nemoh_<hull>_m<k>" per multiple of a hull with HULLS nemoh_batches (OC4: m1 and m2 are approved and run apart)."""
    groups = {}
    for case in cases:
        own = HULLS[case.hull].get("nemoh_batches")
        groups.setdefault("nemoh_%s_%s" % (case.hull, tag(case.multiple)) if own else "nemoh", []).append(case)
    return sorted(groups.items())


def job_table(tier, markdown, control=None, qtf=False):
    cases = [case for case in cases_of(control, qtf) if in_tier(case, tier)]
    head = ("Case", "Group", "Q", "Hull", "Motion", "Velocity", "Setting", "Panels", "Freq", "Head", "Wall (h)", "Mem (GB)", "Runs on")
    rows = []
    for case in cases:
        setting = []
        if case.grn == 8:
            setting.append("GRN_KIND 8")
        if case.isolv == 2:
            setting.append("ISOLV 2")
        if case.irsp:
            setting.append("IRSP 1")
        if case.waterline == "field":
            setting.append("field")
        if case.mean_drift == 2:
            setting.append("pairs")
        if case.control:
            setting.append("ctrl-" + case.control)
        if case.qtf:
            setting.append("qtf %d" % case.qtf)
        if case.qtf_terms:
            setting.append("qtf_terms %s, %d Omega" % (case.qtf_terms, len(difference_frequencies(case))))
            if case.solver == "hams" and HULLS[case.hull].get("cf_points"):
                setting.append("C+F %d" % HULLS[case.hull]["cf_points"][case.multiple])
        if case.solver == "nemoh":
            setting.append("Nemoh" + (" DUOK+HASBO" if case.nemoh_contrib == 2 else ""))
        count = case.nemoh_range[2] if case.solver == "nemoh" else len(case.frequencies)
        rows.append((case.name, case.group, case.question, case.hull, case.motion,
                     "-" if case.solver == "nemoh" else case.velocity, ", ".join(setting) or "-",
                     "%d" % (panels_of(case) + lid_panels_of(case)), "%d" % count, "%d" % len(case.headings),
                     "%.2f" % wall_hours(case), "%.2f" % memory_gb(case), case.where))
    if markdown:
        print("| " + " | ".join(head) + " |")
        print("| " + " | ".join(":---" if i < 7 or i == 12 else "---:" for i in range(len(head))) + " |")
        for row in rows:
            print("| " + " | ".join(row) + " |")
    else:
        for row in rows:
            print("  ".join(row))
    print()
    for label, chosen in (("HAMS on Snellius, 32 threads each", [c for c in cases if c.solver == "hams"]),
                          ("Nemoh on Snellius, 1 thread each", [c for c in cases if c.solver == "nemoh" and c.where == "snellius"]),
                          ("Nemoh on the laptop", [c for c in cases if c.where == "laptop"])):
        hours = sum(wall_hours(c) for c in chosen)
        longest = max([wall_hours(c) for c in chosen] or [0.0])
        print("%s: %d runs, %.1f wall h in sum, longest %.2f h, largest memory %.1f GB" % (
            label, len(chosen), hours, longest, max([memory_gb(c) for c in chosen] or [0.0])))
    hams_hours = sum(wall_hours(c) for c in cases if c.solver == "hams")
    print("Charged: HAMS %.1f job hours, %.0f CPU-hours, %.0f SBU" % (hams_hours, hams_hours * CHARGED_CPUS, hams_hours * SBU_PER_JOB_HOUR))
    for label, chosen in nemoh_batches([c for c in cases if c.solver == "nemoh" and c.where == "snellius"]):
        hours = max(wall_hours(c) for c in chosen)
        print("Charged: Nemoh allocation %s (%d runs, %d CPUs) %.1f job hours, %.0f CPU-hours, %.0f SBU" % (
            label, len(chosen), NEMOH_CPUS, hours, hours * NEMOH_CPUS, hours * NEMOH_CPUS * SBU_PER_JOB_HOUR / CHARGED_CPUS))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")
    sub.required = True

    control_help = "the control-surface study: all frequencies, or `cheap` with three at the m4 cases of the green route"
    qtf_help = "the QTF study: quadratic QTFs on the Nemoh grid (qtf_cases)"
    potential_help = "the QTF potential study: QTF_Indirect.11d with qtf_terms = body, and Nemoh DUOK + HASBO (qtf_potential_cases)"

    p = sub.add_parser("list")
    p.add_argument("--tier", choices=("full", "small"), default="full")
    p.add_argument("--markdown", action="store_true")
    p.add_argument("--control", nargs="?", const="full", choices=("full", "cheap"), help=control_help)
    p.add_argument("--qtf", action="store_true", help=qtf_help)
    p.add_argument("--qtf-potential", action="store_true", help=potential_help)

    p = sub.add_parser("generate")
    p.add_argument("root")
    p.add_argument("--case")
    p.add_argument("--group")
    p.add_argument("--control", nargs="?", const="full", choices=("full", "cheap"), help=control_help)
    p.add_argument("--qtf", action="store_true", help=qtf_help)
    p.add_argument("--qtf-potential", action="store_true", help=potential_help)

    for command in ("run", "calibrate"):
        p = sub.add_parser(command)
        p.add_argument("root")
        p.add_argument("--hams", required=True)
        p.add_argument("--threads", type=int, default=os.cpu_count())
        p.add_argument("--hams-grn8")
        p.add_argument("--nemoh-build")
        if command == "run":
            p.add_argument("--case", required=True)
            p.add_argument("--control", nargs="?", const="full", choices=("full", "cheap"), help=control_help)
            p.add_argument("--qtf", action="store_true", help=qtf_help)
            p.add_argument("--qtf-potential", action="store_true", help=potential_help)

    p = sub.add_parser("slurm")
    p.add_argument("root", help="folder of the study on the cluster")
    p.add_argument("--repo", required=True, help="checkout of this repository on the cluster")
    p.add_argument("--source", help="local checkout to read the source list from, by default this one")
    p.add_argument("--out", help="local folder for the scripts, by default <root>/jobs")
    p.add_argument("--tier", choices=("full", "small"), default="full")
    p.add_argument("--account", help="Slurm account, by default that of the user")
    p.add_argument("--control", nargs="?", const="full", choices=("full", "cheap"), help=control_help)
    p.add_argument("--qtf", action="store_true", help=qtf_help)
    p.add_argument("--qtf-potential", action="store_true", help=potential_help)

    p = sub.add_parser("collect")
    p.add_argument("root")
    p.add_argument("--out")
    p.add_argument("--reference", help="results of the first study, for near-field at m4 where <root> does not hold it")
    p.add_argument("--control", nargs="?", const="full", choices=("full", "cheap"), help=control_help)
    p.add_argument("--qtf", action="store_true", help=qtf_help)
    p.add_argument("--qtf-potential", action="store_true", help=potential_help)

    args = parser.parse_args()
    if args.command == "list":
        job_table(args.tier, args.markdown, args.control, study_of(args))
    elif args.command == "generate":
        for case in cases_of(args.control, study_of(args)):
            if (args.case and case.name != args.case) or (args.group and case.group != args.group):
                continue
            if case.solver == "nemoh":
                write_nemoh_case(case, args.root)
            else:
                write_hams_case(case, args.root, CPUS)
    elif args.command == "run":
        found = [case for case in cases_of(args.control, study_of(args)) if case.name == args.case]
        if not found:
            raise SystemExit("no case named " + args.case)
        sys.exit(1 if run_case(found[0], args)["code"] else 0)
    elif args.command == "calibrate":
        calibrate(args)
    elif args.command == "slurm":
        write_slurm(args)
    elif args.command == "collect":
        collect(args)


if __name__ == "__main__":
    main()
