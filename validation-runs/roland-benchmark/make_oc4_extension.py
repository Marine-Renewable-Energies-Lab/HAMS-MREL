"""OC4 full WAMIT hull (9208 or 5022 panels) over the frequency range of Roland et al. (2026), 0.25 to 2.55 rad/s.

The existing cases (q_oc4wfull_sym_free_L0 and _L1: hull, C, F and the WAMIT-matched Hydrostatic.in) stop at 1.90 rad/s
because velocity_method = source cannot remove irregular frequencies; the first is about 1.98 rad/s for the offset columns
(radius 6 m), the next 2.50 rad/s (j_11) [derived, snellius/COST.md]. This case copies their meshes and matrices unchanged,
and changes only:
  - ControlFile.in: 47 frequencies, 0.25 to 2.55 rad/s in steps of 0.05, If_remove_irr_freq = 1, 24 threads;
  - SecondOrder.in: velocity_method = differencing (source rejects IRSP = 1); every other key as before;
  - WaterPlaneMesh.pnl: a lid on every waterline loop (make_volturnus_case.lid_mesh): a ring on the waterline nodes one
    waterline row high, rings of about 0.5 m inside, and a fan; its nearest centroid must be at least 0.3 waterline row
    heights from the wall (input_files.md, "Lid mesh requirements").
A full rerun is needed: a QTF pair needs the first-order states of both its frequencies in one run.

Usage: python -I make_oc4_extension.py EXISTING_CASE_DIR OUT_CASE_DIR [THREADS]
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_volturnus_case as mv  # noqa: E402

mm = mv.mm
FMAX = 2.55
LID_SPACING = 0.5   # m


def main(argv):
    src, out = os.path.join(argv[1], "input"), os.path.join(argv[2], "input")
    threads = int(argv[3]) if len(argv) > 3 else 24
    os.makedirs(out, exist_ok=True)
    for name in ("HullMesh.pnl", "CSMesh_1.pnl", "FSMesh_1.pnl", "Hydrostatic.in"):
        shutil.copyfile(os.path.join(src, name), os.path.join(out, name))
    old = open(os.path.join(src, "SecondOrder.in")).read()
    if "velocity_method = source" not in old or "qtf_terms = body" not in old:
        raise ValueError("unexpected SecondOrder.in in " + src)
    new = old.replace("velocity_method = source", "velocity_method = differencing")
    new = "! OC4 full WAMIT hull over 0.25-2.55 rad/s (Roland et al. 2026 range); make_oc4_extension.py\n" + new
    open(os.path.join(out, "SecondOrder.in"), "w").write(new)
    depth_line = [l for l in open(os.path.join(src, "ControlFile.in")) if "Waterdepth" in l][0]
    if "200" not in depth_line:
        raise ValueError("expected depth 200 m")
    mv.write_control_file(os.path.join(out, "ControlFile.in"), threads, mv.frequencies(FMAX), 1,
                          note="make_oc4_extension.py")
    hull = mm.read_mesh(os.path.join(out, "HullMesh.pnl"))
    loops = mv.waterline_loops(hull)
    row = mv.waterline_row_height(hull)
    lid, worst = mv.lid_mesh(loops, row, LID_SPACING)
    if worst < 0.3:
        raise ValueError("a lid centroid is %.2f waterline row heights from the wall (at least 0.3 wanted)" % worst)
    mm.write_mesh(os.path.join(out, "WaterPlaneMesh.pnl"), lid, title="Waterplane Mesh File")
    for name in os.listdir(out):
        path = os.path.join(out, name)
        data = open(path, "rb").read().replace(b"\r\n", b"\n")   # read first: opening for writing empties the file
        open(path, "wb").write(data)
    empty = [name for name in os.listdir(out) if os.path.getsize(os.path.join(out, name)) == 0]
    if empty:
        raise ValueError("empty input files: " + ", ".join(empty))
    print("hull: %d panels; %d waterline loops:" % (len(hull.panels), len(loops)))
    for pts, (x, y, r) in loops:
        print("  centre (%.3f, %.3f), mean radius %.3f m, %d nodes" % (x, y, r, len(pts)))
    print("lid: %d panels, waterline row height %.3f m, nearest centroid %.2f row heights from the waterline"
          % (len(lid.panels), row, worst))
    print("frequencies: %d, %.2f to %.2f rad/s" % (len(mv.frequencies(FMAX)), 0.25, FMAX))


if __name__ == "__main__":
    main(sys.argv)
