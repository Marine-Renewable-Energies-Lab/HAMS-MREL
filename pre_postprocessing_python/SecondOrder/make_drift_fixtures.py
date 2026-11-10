"""
Regenerate the mean drift integration fixtures test/data/drift_*/input and test/data/drift_cyl/field_panels.txt.

Usage:
    python make_drift_fixtures.py [path to test/data]

The cases are described in test/data/drift_cyl/README.md. The default path is that of this repository.
"""

import math
import os
import shutil
import sys

import make_mesh

DRIFT_CYL_FREQUENCIES = (0.70, 0.95, 1.10)
DRIFT_CYL_HEADINGS = (0.0, 30.0)
DRIFT_CYL_MESH = dict(radius=10.0, draft=10.0, ntheta=32, nz=8, nr=6, graded=True)
# The hull of drift_cyl with a rounded bilge: near-field drift converges slowly at the sharp bilge (see the README)
ROUND_MESH = dict(radius=10.0, draft=10.0, bilge_radius=3.0, ntheta=32, nz=7, narc=4, nr=4)
TWO_ROUND_MESH = dict(radius=5.0, draft=5.0, bilge_radius=1.5, ntheta=20, nz=5, narc=3, nr=3)
ELLIPSOID_MESH = dict(ax=20.0, ay=10.0, draft=5.0, ntheta=48, nphi=12)
# On external springs: surge at about 0.43 rad/s, heave moved from 0.80 to about 0.86 rad/s, between the fixture frequencies
EXTERNAL_RESTORING = {0: 1.0e6, 2: 0.5e6}
# Removal of irregular frequencies: the last one is the first irregular frequency of the cylinder, J0(ka) = 0
IRSP_FREQUENCIES = (0.95, 1.10, 1.548)
TWO_CYL_MESH = dict(radius=5.0, draft=5.0, ntheta=20, nz=5, nr=4, graded=True)
# 50 wavelengths at this frequency, deep water wave length 2 pi g / omega^2
FAR_APART_FREQUENCY = 1.35
FAR_APART_SEPARATION = 50.0 * 2.0 * math.pi * 9.80665 / FAR_APART_FREQUENCY ** 2
DRIFT_BOX_MESH = dict(lx=20.0, ly=10.0, draft=5.0, nx=20, ny=10, nz=4)
OBLIQUE_FREQUENCIES = (0.8, 1.0)   # of the box and the ellipsoid
ROTATED_BOX_MESH = dict(lx=20.0, ly=10.0, draft=5.0, nx=10, ny=5, nz=3)
ROTATED_DEG = 30.0
# Bottom-mounted: the wall reaches the sea bed, with no bottom panels, so MacCamy-Fuchs is exact
BOTTOM_MOUNTED_MESH = dict(radius=5.0, draft=10.0, ntheta=64, nz=16, nr=0)
BOTTOM_MOUNTED_FREQUENCIES = (1.0, 1.6)
# Fractions of the waterline segment length
WATERLINE_OFFSETS = (0.005, 0.01, 0.05, 0.1)
# Hull panels whose centroids are the field points of drift_cyl: side, bottom, bilge, lee and weather
FIELD_PANELS = (1, 4, 8, 9, 14, 40, 85, 130, 133, 200, 270, 448)
# Control-surface fixtures: threads, the frequency of the cases that need one, and the finite-difference steps besides 0.02
CONTROL_THREADS = 4
CONTROL_FREQUENCY = 0.95
CONTROL_STEPS = (0.005, 0.01, 0.05, 0.1)
# Main grid of qtf_round_low: differences 0.05 to 0.20 rad/s
QTF_LOW_FREQUENCIES = (0.40, 0.45, 0.50, 0.55, 0.60)
# QTF potential part (qtf_terms = body): the pair of qtf_round_body at 0.15 rad/s for IRSP = 1 and the half mesh, and two
# small rounded hulls 2000 m apart against one alone, at 0.20 rad/s. At 0.80 and 1.00 rad/s (ka 0.33 to 0.51) the wave
# scattered by the other hull changes .10d by 0.8% at 2000 m; at 1.35 and 1.55 rad/s it was 5.9%, too much for a 2% check
QTF_BODY_FREQUENCIES = (0.95, 1.10)
QTF_DISTANT_FREQUENCIES = (0.80, 1.00)
QTF_DISTANT_SEPARATION = 2000.0
# Field points of QTF_FieldPartial.txt on the rounded hull (radius 10, draft 10): free-surface points and points below z = 0,
# on the weather side, the lee side, ahead and astern, and two with y < 0 for the half mesh
QTF_FIELD_POINTS = ((30.0, 0.0, 0.0), (0.0, 40.0, 0.0), (-50.0, 20.0, 0.0), (30.0, 0.0, -5.0), (0.0, -40.0, -10.0),
                    (60.0, 60.0, -20.0))
# A hull of radius and draft 0.05 m in 20 m of water at 0.9 rad/s has ka of about 0.004: the field is that of the incident
# wave alone, and the Stokes wave of the incident pair is the reference
QTF_TINY_MESH = dict(radius=0.05, draft=0.05, ntheta=8, nz=2, nr=2)
QTF_TINY_FREQUENCY = 0.9
QTF_TINY_DEPTH = 20.0
QTF_TINY_POINTS = ((30.0, 0.0, 0.0), (30.0, 10.0, 0.0), (30.0, 0.0, -3.0))


def real12(value):
    """A real for a 12 character field of ReadControlFile. The reader takes all the digits that fit in the field."""
    return "%12.3f" % value if abs(value - round(value, 3)) < 1e-12 else ("%.9f" % value)[:12].rjust(12)


def write_control(path, depth=100.0, frequencies=(1.0,), headings=(0.0,), rotation_centre=(0.0, 0.0, 0.0), refl="1.0D0",
                  sybo=0, solver=None, field_points=(), body_origins=None, lcs_angle=0.0, irsp=0, threads=1):
    """Write ControlFile.in, in the fixed layout read by ReadControlFile. body_origins lists (x, y) per body of a
    multi-body case, which is both the LCS origin and the rotation centre. lcs_angle is the rotation of every body's
    LCS in radians, which is what ReadBodyMeshMulti expects."""
    with open(path, "w") as f:
        f.write("   --------------HAMS Control file---------------\n\n")
        f.write("   Waterdepth  %.1fD0\n\n" % depth)
        f.write("   #Start Definition of Wave Frequencies\n")
        f.write("    0_inf_frequency_limits      %d\n" % sybo)
        f.write("    Input_frequency_type        3\n    Output_frequency_type       3\n")
        f.write("    Number_of_frequencies      %d\n" % len(frequencies))
        f.write("  " + " ".join("%.4f" % x for x in frequencies) + "\n")
        f.write("   #End Definition of Wave Frequencies\n\n")
        if body_origins:
            f.write("   # Start of definition of bodies for multi-body interaction\n")
            f.write("    Number_of_bodies          %d\n" % len(body_origins))
            for i, (x, y) in enumerate(body_origins, start=1):
                f.write(("    LCS_%d" % i).ljust(26) + real12(x) + real12(y) + real12(0.0) + real12(lcs_angle) + "\n")
            f.write("   # End of definition of bodies for multi-body interaction\n\n")
        f.write("   #Start Definition of Wave Headings\n")
        f.write("    Number_of_headings          %d\n" % len(headings))
        f.write("  " + " ".join("%.4f" % x for x in headings) + "\n")
        f.write("   #End Definition of Wave Headings\n\n")
        if body_origins:
            for i, (x, y) in enumerate(body_origins, start=1):
                f.write(("    Reference_body_center_%d" % i).ljust(28) + real12(x) + real12(y) + real12(0.0) + "\n")
        else:
            f.write("    Reference_body_center".ljust(28) + "%12.3f%12.3f%12.3f\n" % tuple(rotation_centre))
        f.write("    Reference_body_length   %s\n" % refl)
        f.write("    Wave_diffrac_solution    1\n    If_remove_irr_freq       %d\n" % irsp)
        f.write("    Number of threads       %d\n" % threads)
        if solver:
            f.write("    Solver_type              %d\n" % solver)
        f.write("\n   #Start Definition of Pressure and/or Elevation (PE)\n")
        f.write("    Number_of_field_points     %d\n" % len(field_points))
        for point in field_points:
            f.write("    %.17g  %.17g  %.17g\n" % tuple(point))
        f.write("   #End Definition of Pressure and/or Elevation\n\n")
        f.write("   ----------End HAMS Control file---------------\n\n")
        f.write("   Note: generated by pre_postprocessing_python/SecondOrder/make_drift_fixtures.py.\n")


def write_second_order(path, **options):
    with open(path, "w") as f:
        f.write("! Mean drift options of this integration test\n")
        for key, value in options.items():
            f.write("%s = %s\n" % (key, value))


def case_dir(data, name):
    path = os.path.join(data, name, "input")
    os.makedirs(path, exist_ok=True)
    return path


def panel_centroid(mesh, number):
    """Centroid as CalPanelCentre computes it from the node coordinates as written to the mesh file, so that a field
    point written at full precision lies on the panel as HAMS sees it."""
    nodes = mesh.panels[number - 1]
    total = [0.0, 0.0, 0.0]
    for node in nodes:
        for k in range(3):
            total[k] += float(make_mesh.NODE_FORMAT % mesh.nodes[node - 1][k])
    return [value / len(nodes) for value in total]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "..", "test", "data")

    both = "far near"

    # Single cylinder with field points on the hull
    path = case_dir(data, "drift_cyl")
    mesh = make_mesh.write_cylinder(os.path.join(path, "HullMesh.pnl"), os.path.join(path, "Hydrostatic.in"),
                                    **DRIFT_CYL_MESH)
    write_control(os.path.join(path, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES, headings=DRIFT_CYL_HEADINGS,
                  field_points=[panel_centroid(mesh, p) for p in FIELD_PANELS])
    write_second_order(os.path.join(path, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                       write_terms=2)
    with open(os.path.join(data, "drift_cyl", "field_panels.txt"), "w") as f:
        f.write("! field point number in ControlFile.in, hull panel whose centroid it is\n")
        for point, panel in enumerate(FIELD_PANELS, start=1):
            f.write("%d %d\n" % (point, panel))

    # The same case with the limit frequencies and all heading pairs, with the iterative solver, and held fixed
    for name, control, options in (("drift_cyl_sybo", dict(sybo=1), dict(drift_methods=both, motion_state="free", write_terms=0,
                                                                         mean_drift=2)),
                                   ("drift_cyl_gmres", dict(solver=2), dict(drift_methods="far", motion_state="free", write_terms=0,
                                                                            mean_drift=1)),
                                   ("drift_cyl_fixed", dict(), dict(mean_drift=1, drift_methods=both, write_terms=1,
                                                                    motion_state="fixed"))):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(path, file), os.path.join(target, file))
        write_control(os.path.join(target, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES,
                      headings=DRIFT_CYL_HEADINGS, **control)
        write_second_order(os.path.join(target, "SecondOrder.in"), **options)

    # The cylinder with a rounded bilge: free, held fixed, and free on external surge and heave springs.
    for name, motion, external in (("drift_round", "free", None),
                                   ("drift_round_fixed", "fixed", None),
                                   ("drift_round_kstf", "free", EXTERNAL_RESTORING)):
        target = case_dir(data, name)
        make_mesh.write_rounded_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                         external_restoring=external, **ROUND_MESH)
        write_control(os.path.join(target, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES,
                      headings=DRIFT_CYL_HEADINGS)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state=motion,
                           write_terms=1)

    # The same cylinder as a half mesh with a y symmetry plane
    target = case_dir(data, "drift_cyl_half")
    make_mesh.write_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"), half_y=True,
                             **DRIFT_CYL_MESH)
    write_control(os.path.join(target, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES, headings=DRIFT_CYL_HEADINGS)
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                       write_terms=0)

    # That cylinder with and without irregular frequency removal, at two frequencies of drift_cyl and at the first irregular
    # frequency. The lid is uniform.
    for name, irsp in (("drift_irsp", 1), ("drift_irsp_ref", 0)):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(path, file), os.path.join(target, file))
        if irsp == 1:
            make_mesh.write_waterplane(os.path.join(target, "WaterPlaneMesh.pnl"), DRIFT_CYL_MESH["radius"],
                                       DRIFT_CYL_MESH["ntheta"], DRIFT_CYL_MESH["nr"])
        write_control(os.path.join(target, "ControlFile.in"), frequencies=IRSP_FREQUENCIES, headings=(0.0,), irsp=irsp)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                           write_terms=0)

    # A small cylinder off the origin, so that M_z is not zero, for the REFL scaling and the near-field M_z
    for name, refl, methods in (("drift_offset", "1.0D0", both), ("drift_offset_refl2", "2.0D0", "far")):
        target = case_dir(data, name)
        make_mesh.write_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                 radius=5.0, draft=5.0, ntheta=16, nz=3, nr=3, x0=6.0, y0=8.0)
        write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.2,), headings=(30.0,),
                      rotation_centre=(6.0, 8.0, 0.0), refl=refl)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=methods,
                           motion_state="fixed", write_terms=0)

    # A box at an oblique heading has a yaw moment of its own
    target = case_dir(data, "drift_box")
    make_mesh.write_box(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"), **DRIFT_BOX_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=OBLIQUE_FREQUENCIES, headings=(30.0,))
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="fixed",
                       write_terms=1)

    # So does a half ellipsoid, which has no edges
    target = case_dir(data, "drift_ellipsoid")
    make_mesh.write_half_ellipsoid(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                   **ELLIPSOID_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=OBLIQUE_FREQUENCIES, headings=(30.0,))
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="fixed",
                       write_terms=1)

    # Two cylinders, multi-body path
    target = case_dir(data, "drift_two_cyl")
    for body in (1, 2):
        make_mesh.write_cylinder(os.path.join(target, "HullMesh_%d.pnl" % body),
                                 os.path.join(target, "Hydrostatic_%d.in" % body), **TWO_CYL_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.35, 1.55), headings=(0.0, 30.0),
                  body_origins=[(0.0, 0.0), (0.0, 30.0)])
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=2, drift_methods=both, motion_state="free",
                       write_terms=2)

    # The same with and without irregular frequency removal, uniform lids
    for name, irsp in (("drift_two_irsp", 1), ("drift_two_irsp_ref", 0)):
        target = case_dir(data, name)
        for body in (1, 2):
            make_mesh.write_cylinder(os.path.join(target, "HullMesh_%d.pnl" % body),
                                     os.path.join(target, "Hydrostatic_%d.in" % body), **TWO_CYL_MESH)
            if irsp == 1:
                make_mesh.write_waterplane(os.path.join(target, "WaterPlaneMesh_%d.pnl" % body), TWO_CYL_MESH["radius"],
                                           TWO_CYL_MESH["ntheta"], TWO_CYL_MESH["nr"])
        write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.35,), headings=(0.0, 30.0),
                      body_origins=[(0.0, 0.0), (0.0, 30.0)], irsp=irsp)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                           write_terms=0)

    # Two cylinders with rounded bilges
    target = case_dir(data, "drift_two_round")
    for body in (1, 2):
        make_mesh.write_rounded_cylinder(os.path.join(target, "HullMesh_%d.pnl" % body),
                                         os.path.join(target, "Hydrostatic_%d.in" % body), **TWO_ROUND_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.35, 1.55), headings=(0.0, 30.0),
                  body_origins=[(0.0, 0.0), (0.0, 30.0)])
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=2, drift_methods=both, motion_state="free",
                       write_terms=0)

    # Two of the sharp cylinders far apart, and one alone
    target = case_dir(data, "drift_far_apart")
    for body in (1, 2):
        make_mesh.write_cylinder(os.path.join(target, "HullMesh_%d.pnl" % body),
                                 os.path.join(target, "Hydrostatic_%d.in" % body), **TWO_CYL_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(FAR_APART_FREQUENCY,), headings=(0.0,),
                  body_origins=[(0.0, 0.0), (0.0, FAR_APART_SEPARATION)])
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                       write_terms=0)
    target = case_dir(data, "drift_isolated")
    make_mesh.write_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"), **TWO_CYL_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(FAR_APART_FREQUENCY,), headings=(0.0,))
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods="near", motion_state="free",
                       write_terms=0)

    # A surface-piercing and a submerged body
    target = case_dir(data, "drift_mixed")
    make_mesh.write_cylinder(os.path.join(target, "HullMesh_1.pnl"), os.path.join(target, "Hydrostatic_1.in"),
                             **TWO_CYL_MESH)
    make_mesh.write_submerged_cylinder(os.path.join(target, "HullMesh_2.pnl"), os.path.join(target, "Hydrostatic_2.in"),
                                       radius=3.0, top=4.0, height=4.0, ntheta=16, nz=3, nr=3)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.35,), headings=(0.0,),
                  body_origins=[(0.0, 0.0), (0.0, 20.0)])
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                       write_terms=1)

    # Two fixed boxes, not in line with the wave, and the whole layout rotated about the origin through the LCS, with the
    # wave
    angle = math.radians(ROTATED_DEG)
    x2, y2 = 30.0, 25.0
    for name, origins, lcs_angle, heading in (
            ("drift_unrotated", [(0.0, 0.0), (x2, y2)], 0.0, 0.0),
            ("drift_rotated", [(0.0, 0.0), (x2 * math.cos(angle) - y2 * math.sin(angle),
                                            x2 * math.sin(angle) + y2 * math.cos(angle))], angle, ROTATED_DEG)):
        target = case_dir(data, name)
        for body in (1, 2):
            make_mesh.write_box(os.path.join(target, "HullMesh_%d.pnl" % body),
                                os.path.join(target, "Hydrostatic_%d.in" % body), **ROTATED_BOX_MESH)
        write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.0,), headings=(heading,),
                      body_origins=origins, lcs_angle=lcs_angle)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="fixed",
                           write_terms=0)

    # The waterline potential by field evaluation, on the bottom-mounted cylinder at each offset
    for offset in WATERLINE_OFFSETS:
        target = case_dir(data, "drift_field_%s" % str(offset).replace(".", "p"))
        make_mesh.write_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                 **BOTTOM_MOUNTED_MESH)
        write_control(os.path.join(target, "ControlFile.in"), depth=BOTTOM_MOUNTED_MESH["draft"],
                      frequencies=BOTTOM_MOUNTED_FREQUENCIES, headings=(0.0,))
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods="near", motion_state="fixed",
                           write_terms=1, waterline_method="field", waterline_offset=offset)

    # The same method with and without irregular frequency removal
    for name, irsp in (("drift_field_irsp", 1), ("drift_field_irsp_ref", 0)):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(path, file), os.path.join(target, file))
        if irsp == 1:
            make_mesh.write_waterplane(os.path.join(target, "WaterPlaneMesh.pnl"), DRIFT_CYL_MESH["radius"],
                                       DRIFT_CYL_MESH["ntheta"], DRIFT_CYL_MESH["nr"])
        write_control(os.path.join(target, "ControlFile.in"), frequencies=IRSP_FREQUENCIES[:2], headings=(0.0,), irsp=irsp)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods="near", motion_state="free",
                           write_terms=1, waterline_method="field")

    # The same method through the multi-body path, on the hulls of drift_two_cyl
    target = case_dir(data, "drift_field_two")
    for body in (1, 2):
        make_mesh.write_cylinder(os.path.join(target, "HullMesh_%d.pnl" % body),
                                 os.path.join(target, "Hydrostatic_%d.in" % body), **TWO_CYL_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=(1.35,), headings=(0.0, 30.0),
                  body_origins=[(0.0, 0.0), (0.0, 30.0)])
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods="near", motion_state="free",
                       write_terms=1, waterline_method="field")

    # The velocity by source solve, on the hulls above: sharp and rounded, free and fixed, half mesh, deep water, two bodies
    # Frequencies are cut to keep the suite short, with every mesh as in the hulls above: the sharp hull at 0.95 rad/s (free) and
    # at 0.70 rad/s (fixed), the fixed rounded hull at 0.70 and 1.10, deep water at 1.10, and the free rounded hull at all three
    for name, source, motion, write_terms, depth, frequencies in (
            ("drift_source_cyl", "drift_cyl", "free", 2, 100.0, (0.95,)),
            ("drift_source_cyl_half", "drift_cyl_half", "free", 0, 100.0, (0.95,)),
            ("drift_source_cyl_fixed", "drift_cyl_fixed", "fixed", 0, 100.0, (0.70,)),
            ("drift_source_round", "drift_round", "free", 2, 100.0, DRIFT_CYL_FREQUENCIES),
            ("drift_source_round_fixed", "drift_round_fixed", "fixed", 0, 100.0, (0.70, 1.10)),
            ("drift_source_deep", "drift_round_fixed", "fixed", 0, -1.0, (1.10,))):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(data, source, "input", file), os.path.join(target, file))
        write_control(os.path.join(target, "ControlFile.in"), depth=depth, frequencies=frequencies,
                      headings=DRIFT_CYL_HEADINGS)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state=motion,
                           write_terms=write_terms, velocity_method="source")
    target = case_dir(data, "drift_source_two_round")
    for body in (1, 2):
        for file in ("HullMesh_%d.pnl" % body, "Hydrostatic_%d.in" % body):
            shutil.copyfile(os.path.join(data, "drift_two_round", "input", file), os.path.join(target, file))
    shutil.copyfile(os.path.join(data, "drift_two_round", "input", "ControlFile.in"), os.path.join(target, "ControlFile.in"))
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=2, drift_methods=both, motion_state="free",
                       write_terms=0, velocity_method="source")

    # Quadratic QTFs: the free rounded hull with the limit frequencies (which must stay out of the QTF grid) and a
    # qtf_max_dw that drops the pair (0.70, 1.10) from .10d, and the two rounded hulls with every heading pair
    target = case_dir(data, "qtf_round")
    for file in ("HullMesh.pnl", "Hydrostatic.in"):
        shutil.copyfile(os.path.join(data, "drift_round", "input", file), os.path.join(target, file))
    write_control(os.path.join(target, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES, headings=DRIFT_CYL_HEADINGS,
                  sybo=1)
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                       write_terms=0, velocity_method="source", qtf=3, qtf_max_dw=0.3)
    # The radiation potentials of the indirect method (qtf_terms = body): the qtf_round case with the check file, and the same
    # hull on a low-frequency grid whose differences go down to 0.05 rad/s, with the zero frequency of SYBO = 1 to compare
    for name, frequencies, options in (("qtf_round_body", DRIFT_CYL_FREQUENCIES, dict(qtf=3, qtf_max_dw=0.3)),
                                       ("qtf_round_low", QTF_LOW_FREQUENCIES, dict(qtf=1))):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(data, "drift_round", "input", file), os.path.join(target, file))
        write_control(os.path.join(target, "ControlFile.in"), frequencies=frequencies, headings=DRIFT_CYL_HEADINGS, sybo=1)
        write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=1, drift_methods=both, motion_state="free",
                           write_terms=1, velocity_method="source", qtf_terms="body", **options)
    target = case_dir(data, "qtf_two_round")
    for file in ("HullMesh_1.pnl", "HullMesh_2.pnl", "Hydrostatic_1.in", "Hydrostatic_2.in", "ControlFile.in"):
        shutil.copyfile(os.path.join(data, "drift_two_round", "input", file), os.path.join(target, file))
    write_second_order(os.path.join(target, "SecondOrder.in"), mean_drift=2, drift_methods=both, motion_state="free",
                       write_terms=0, velocity_method="source", qtf=3)
    write_qtf_body_fixtures(data, both)
    write_qtf_field_fixtures(data, both)

    # The control-surface method on the hull of drift_cyl: held fixed, and with the limit frequencies and all heading pairs
    for name, control, options in (("control_cyl", dict(), dict(mean_drift=1, drift_methods="far control", write_terms=1,
                                                                motion_state="fixed")),
                                   ("control_cyl_sybo", dict(sybo=1), dict(drift_methods="far near control",
                                                                           motion_state="free", write_terms=0, mean_drift=2))):
        target = case_dir(data, name)
        for file in ("HullMesh.pnl", "Hydrostatic.in"):
            shutil.copyfile(os.path.join(path, file), os.path.join(target, file))
        write_control(os.path.join(target, "ControlFile.in"), frequencies=DRIFT_CYL_FREQUENCIES,
                      headings=DRIFT_CYL_HEADINGS, **control)
        write_second_order(os.path.join(target, "SecondOrder.in"), **options)

    write_control_physics(data)


def write_qtf_body_fixtures(data, both):
    """Fixtures of QTF_Indirect.11d (qtf_terms = body), besides qtf_round_body and qtf_round_low. Kept small: two
    frequencies where one pair is enough."""
    body = dict(mean_drift=1, drift_methods=both, motion_state="free", write_terms=1, qtf_terms="body", qtf=1)
    # The rounded hull of qtf_round_body with irregular frequency removal, as a half mesh, and the whole mesh without
    # removal as their reference. Velocity by differencing in all three, since the source solve does not take IRSP = 1.
    for name, irsp, half in (("qtf_round_body_ref", 0, False), ("qtf_round_body_irsp", 1, False),
                             ("qtf_half_body", 0, True)):
        target = case_dir(data, name)
        make_mesh.write_rounded_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                         half_y=half, **ROUND_MESH)
        if irsp == 1:
            make_mesh.write_waterplane(os.path.join(target, "WaterPlaneMesh.pnl"), ROUND_MESH["radius"],
                                       ROUND_MESH["ntheta"], ROUND_MESH["nr"])
        write_control(os.path.join(target, "ControlFile.in"), frequencies=QTF_BODY_FREQUENCIES, headings=DRIFT_CYL_HEADINGS,
                      irsp=irsp)
        write_second_order(os.path.join(target, "SecondOrder.in"), **body)
    # Two of the small rounded hulls far apart, with every heading pair, and one alone
    for name, origins, pairs in (("qtf_two_round_body", [(0.0, 0.0), (0.0, QTF_DISTANT_SEPARATION)], 2),
                                 ("qtf_isolated_body", None, 1)):
        target = case_dir(data, name)
        for b in ((1, 2) if origins else (0,)):
            suffix = "_%d" % b if b else ""
            make_mesh.write_rounded_cylinder(os.path.join(target, "HullMesh%s.pnl" % suffix),
                                             os.path.join(target, "Hydrostatic%s.in" % suffix), **TWO_ROUND_MESH)
        write_control(os.path.join(target, "ControlFile.in"), depth=50.0, frequencies=QTF_DISTANT_FREQUENCIES,
                      headings=DRIFT_CYL_HEADINGS, body_origins=origins)
        write_second_order(os.path.join(target, "SecondOrder.in"), **dict(body, mean_drift=pairs, velocity_method="source"))
    # Every second frequency of qtf_round_low: psi is solved at each difference frequency, so the grid spacing must not
    # change the common pairs
    target = case_dir(data, "qtf_round_low_spacing")
    for file in ("HullMesh.pnl", "Hydrostatic.in"):
        shutil.copyfile(os.path.join(data, "drift_round", "input", file), os.path.join(target, file))
    write_control(os.path.join(target, "ControlFile.in"), frequencies=QTF_LOW_FREQUENCIES[::2], headings=DRIFT_CYL_HEADINGS)
    write_second_order(os.path.join(target, "SecondOrder.in"), **dict(body, velocity_method="source"))


def write_qtf_field_fixtures(data, both):
    """Fixtures of QTF_FieldPartial.txt (qtf_field = 1): the rounded hull with the limit frequencies and every heading pair,
    its half mesh against the whole mesh (the set-up of qtf_half_body), a hull too small to scatter, and the two ways the
    field points can stop the run. The reader needs near in drift_methods for every qtf."""
    field = dict(mean_drift=2, drift_methods=both, motion_state="free", write_terms=0, qtf=3, qtf_max_dw=0.3, qtf_field=1)

    def rounded_case(name, frequencies, points, half=False, sybo=0, **options):
        target = case_dir(data, name)
        make_mesh.write_rounded_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"),
                                         half_y=half, **ROUND_MESH)
        write_control(os.path.join(target, "ControlFile.in"), frequencies=frequencies, headings=DRIFT_CYL_HEADINGS,
                      sybo=sybo, field_points=points)
        write_second_order(os.path.join(target, "SecondOrder.in"), **dict(field, **options))

    rounded_case("qtf_field", DRIFT_CYL_FREQUENCIES, QTF_FIELD_POINTS, sybo=1)
    # The half mesh (ISY = 1) and the whole mesh of the same hull, at the frequencies and headings of qtf_half_body
    rounded_case("qtf_field_half", QTF_BODY_FREQUENCIES, QTF_FIELD_POINTS, half=True)
    rounded_case("qtf_field_full", QTF_BODY_FREQUENCIES, QTF_FIELD_POINTS)
    # A point inside the hull, and no points at all
    rounded_case("qtf_field_inside", DRIFT_CYL_FREQUENCIES, QTF_FIELD_POINTS + ((0.0, 0.0, -5.0),), sybo=1)
    rounded_case("qtf_field_no_points", DRIFT_CYL_FREQUENCIES, (), sybo=1)
    target = case_dir(data, "qtf_field_tiny")
    make_mesh.write_cylinder(os.path.join(target, "HullMesh.pnl"), os.path.join(target, "Hydrostatic.in"), **QTF_TINY_MESH)
    write_control(os.path.join(target, "ControlFile.in"), depth=QTF_TINY_DEPTH, frequencies=(QTF_TINY_FREQUENCY,),
                  headings=(0.0,), field_points=QTF_TINY_POINTS)
    write_second_order(os.path.join(target, "SecondOrder.in"), **dict(field, mean_drift=1))


def copy_inputs(data, source, target, files):
    for file in files:
        shutil.copyfile(os.path.join(data, source, "input", file), os.path.join(target, file))


def write_control_physics(data):
    """The control-surface method against far-field, exact results and itself. These run on CONTROL_THREADS threads, with as
    few frequencies and headings as each test needs: the field rows cost about 50 s per frequency on one thread at 448 hull
    panels."""
    hull = ("HullMesh.pnl", "Hydrostatic.in")
    two_hulls = ("HullMesh_1.pnl", "Hydrostatic_1.in", "HullMesh_2.pnl", "Hydrostatic_2.in")
    one = (CONTROL_FREQUENCY,)
    fixed = dict(mean_drift=1, drift_methods="far control", write_terms=0, motion_state="fixed")
    free = dict(fixed, motion_state="free")
    # (name, source fixture, files, ControlFile.in options, SecondOrder.in options)
    cases = [
        ("control_cyl_free", "drift_cyl", hull, dict(frequencies=DRIFT_CYL_FREQUENCIES), free),
        ("control_cyl_wide", "drift_cyl", hull, dict(frequencies=one), dict(fixed, cs_margin=30)),
        ("control_cyl_box", "drift_cyl", hull, dict(frequencies=one), dict(fixed, cs_shape="box")),
        ("control_bottom", "drift_field_0p01", hull, dict(depth=BOTTOM_MOUNTED_MESH["draft"],
                                                          frequencies=BOTTOM_MOUNTED_FREQUENCIES, headings=(0.0,)), fixed),
        ("control_box", "drift_box", hull, dict(depth=50.0, frequencies=OBLIQUE_FREQUENCIES, headings=(30.0,)), fixed),
        # C at half the default panel size: with the default, the sum over the bodies is 2.28% from far-field, with this 0.58%
        ("control_two_round", "drift_two_round", two_hulls,
         dict(depth=50.0, frequencies=(1.35,), body_origins=[(0.0, 0.0), (0.0, 30.0)]),
         dict(free, drift_methods="far near control", cs_panel_size=0.75)),
        ("control_far_apart", "drift_far_apart", two_hulls,
         dict(depth=50.0, frequencies=(FAR_APART_FREQUENCY,), headings=(0.0,),
              body_origins=[(0.0, 0.0), (0.0, FAR_APART_SEPARATION)]), free),
        ("control_isolated", "drift_isolated", hull, dict(depth=50.0, frequencies=(FAR_APART_FREQUENCY,), headings=(0.0,)),
         free),
        ("control_irsp", "drift_irsp", hull + ("WaterPlaneMesh.pnl",),
         dict(frequencies=IRSP_FREQUENCIES, headings=(0.0,), irsp=1), free),
        ("control_irsp_ref", "drift_irsp_ref", hull, dict(frequencies=IRSP_FREQUENCIES, headings=(0.0,)), free),
        ("control_half", "drift_cyl_half", hull, dict(frequencies=one), free),
        ("control_gmres", "drift_cyl", hull, dict(frequencies=one, solver=2), free),
        ("control_pairs", "drift_cyl", hull, dict(frequencies=one), dict(free, mean_drift=2)),
        # The source route on C, with near-field by differencing: the driver builds the sources on a scratch copy of the
        # state. control_source_ref is the same case without control, for the near-field rows
        ("control_source", "drift_cyl", hull, dict(frequencies=DRIFT_CYL_FREQUENCIES),
         dict(free, drift_methods="far near control", velocity_method="differencing", cs_velocity="source")),
        ("control_source_ref", "drift_cyl", hull, dict(frequencies=DRIFT_CYL_FREQUENCIES),
         dict(free, drift_methods="far near", velocity_method="differencing")),
        # Deep water: the Green function goes by INFGREEN3D, as drift_source_deep
        ("control_deep", "drift_round_fixed", hull, dict(depth=-1.0, frequencies=(1.10,)), fixed),
    ]
    # The finite-difference step of the green route. 0.02, the default, is control_cyl
    for step in CONTROL_STEPS:
        cases.append(("control_step_%s" % str(step).replace(".", "p"), "drift_cyl", hull, dict(frequencies=one),
                      dict(fixed, cs_fd_step=step)))
    for name, source, files, control, options in cases:
        target = case_dir(data, name)
        copy_inputs(data, source, target, files)
        control = dict(dict(headings=DRIFT_CYL_HEADINGS, threads=CONTROL_THREADS), **control)
        write_control(os.path.join(target, "ControlFile.in"), **control)
        write_second_order(os.path.join(target, "SecondOrder.in"), **options)


if __name__ == "__main__":
    main()
