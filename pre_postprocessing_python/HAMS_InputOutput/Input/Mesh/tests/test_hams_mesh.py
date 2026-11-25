"""Tests of hams_mesh.py: one per shape (volume and waterplane area against the exact values within the panel error, and
the written files read back in the HAMS format), the dipole and hull checks, the lid, the hydrostatics and the multi-body
set. Run from the repository root or this folder: python -m pytest pre_postprocessing_python/HAMS_InputOutput/Input/Mesh"""

import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import hams_mesh as hm  # noqa: E402

RG = hm.RHO * hm.G


def polygon_area(n, r):
    return 0.5 * n * r * r * math.sin(2.0 * math.pi / n)


def sagitta_bound(body, radius, surface):
    """Volume a polyhedron with its nodes on a sphere of the radius can miss: the surface times the largest gap between a
    panel and the sphere, (panel size)^2 / (2 R)."""
    size = max(hm.panel_size(body.mesh.points(p)) for p in body.mesh.panels)
    return surface * size ** 2 / (2.0 * radius)


def write_and_read(body, folder, **kw):
    """Write the input files and read them back as HAMS does; check that the mesh comes back unchanged."""
    results, lines, warnings = hm.write_input(body, str(folder), **kw)
    nodes, panels, dipoles, sym = hm.read_pnl(os.path.join(str(folder), "HullMesh.pnl"))
    assert sym == (0, 0)
    assert panels == body.mesh.panels and dipoles == body.mesh.dipoles
    assert np.allclose(nodes, body.mesh.nodes, atol=1e-9)
    cog, blocks = hm.read_hydrostatic(os.path.join(str(folder), "Hydrostatic.in"))
    assert np.allclose(blocks[0], results[0]["mass_matrix"], rtol=1e-5, atol=1e-6 * max(1.0, abs(blocks[0]).max()))
    assert np.allclose(blocks[3], results[0]["restoring"], rtol=1e-5, atol=1e-6 * max(1.0, abs(blocks[3]).max()))
    return results[0], lines, warnings


def test_spacing():
    for fs, fe in ((0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.2)):
        s = hm.spacing(10, fs, fe)
        assert s[0] == 0.0 and s[-1] == 1.0 and all(b > a for a, b in zip(s, s[1:]))
    s = hm.spacing(10, 0.0, 1.0)
    assert s[-1] - s[-2] < 0.5 * (s[1] - s[0])
    assert hm.spacing(4) == [0.0, 0.25, 0.5, 0.75, 1.0]


def test_cylinder(tmp_path):
    body = hm.cylinder(5.0, 10.0, panel_size=1.0, waterline_grading=0.5, edge_grading=0.5)
    h, lines, _ = write_and_read(body, tmp_path)
    n = body.ntheta
    assert h["volume"] == pytest.approx(polygon_area(n, 5.0) * 10.0, rel=1e-10)
    assert h["waterplane"] == pytest.approx(polygon_area(n, 5.0), rel=1e-10)
    deficit = 1.0 - polygon_area(n, 5.0) / (math.pi * 25.0)
    assert 0.0 < 1.0 - h["volume"] / body.exact["volume"] <= deficit * (1 + 1e-9)
    # Waterplane second moment of the regular polygon: n R^4 sin(2 pi / n) (2 + cos(2 pi / n)) / 24
    t = 2.0 * math.pi / n
    ixx = n * 5.0 ** 4 * math.sin(t) * (2.0 + math.cos(t)) / 24.0
    assert h["restoring"][3, 3] == pytest.approx(RG * ixx, rel=1e-9)
    assert h["cog"][2] == pytest.approx(-5.0, rel=1e-10)
    assert "Reference_body_center" in lines[0]


def test_cylinder_matches_make_mesh():
    """Same uniform cylinder as SecondOrder/make_mesh.py: the same nodes and the same Hydrostatic.in matrices."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "SecondOrder"))
    import make_mesh
    old, polygon = make_mesh.cylinder(10.0, 10.0, 32, 6, 6)
    new = hm.cylinder(10.0, 10.0, ntheta=32, nz=6, nr=6)
    assert len(new.mesh.panels) == len(old.panels)
    assert sorted(np.round(new.mesh.nodes, 9).tolist()) == sorted(np.round(old.nodes, 9).tolist())
    ref = make_mesh.hydrostatics(polygon, 10.0)
    h = hm.hydrostatics(new.mesh, (0.0, 0.0, 0.0))
    assert np.allclose(h["mass_matrix"], ref["mass"], rtol=1e-9, atol=1e-6)
    assert np.allclose(h["restoring"], ref["restoring"], rtol=1e-9, atol=1e-6)


def test_box(tmp_path):
    body = hm.box(20.0, 10.0, 5.0, panel_size=1.5, waterline_grading=0.3, edge_grading=0.5)
    h, _, _ = write_and_read(body, tmp_path)
    assert h["volume"] == pytest.approx(1000.0, rel=1e-12)
    assert h["waterplane"] == pytest.approx(200.0, rel=1e-12)
    assert h["restoring"][3, 3] == pytest.approx(RG * 20.0 * 10.0 ** 3 / 12.0, rel=1e-10)
    assert h["restoring"][4, 4] == pytest.approx(RG * 10.0 * 20.0 ** 3 / 12.0, rel=1e-10)
    m = 1e6
    assert h["mass_matrix"][3, 3] == pytest.approx(m * (10.0 ** 2 + 5.0 ** 2) / 12.0 + m * 2.5 ** 2, rel=1e-10)
    assert h["mass_matrix"][5, 5] == pytest.approx(m * (20.0 ** 2 + 10.0 ** 2) / 12.0, rel=1e-10)


def test_hemisphere(tmp_path):
    body = hm.hemisphere(5.0, panel_size=0.6, waterline_grading=0.5)
    h, _, _ = write_and_read(body, tmp_path)
    miss = body.exact["volume"] - h["volume"]
    assert 0.0 < miss <= sagitta_bound(body, 5.0, 2.0 * math.pi * 25.0)
    assert h["waterplane"] == pytest.approx(body.exact["waterplane_polygon"], rel=1e-10)
    assert miss / body.exact["volume"] < 0.01


def test_sphere(tmp_path):
    body = hm.sphere(3.0, 8.0, panel_size=0.6)
    h, lines, _ = write_and_read(body, tmp_path)
    miss = body.exact["volume"] - h["volume"]
    assert 0.0 < miss <= sagitta_bound(body, 3.0, 4.0 * math.pi * 9.0)
    assert abs(h["waterplane"]) < 1e-9
    assert abs(h["restoring"]).max() < 1e-6 * RG
    assert h["cog"][2] == pytest.approx(-8.0, abs=1e-9)
    with pytest.raises(ValueError):
        hm.lid(body)


def test_cone(tmp_path):
    body = hm.cone(6.0, 3.0, 8.0, panel_size=0.8, edge_grading=0.5)
    h, _, _ = write_and_read(body, tmp_path)
    assert h["volume"] == pytest.approx(body.exact["volume_polygon"], rel=1e-10)
    assert h["waterplane"] == pytest.approx(body.exact["waterplane_polygon"], rel=1e-10)
    assert 0.0 < body.exact["volume"] - h["volume"] <= body.exact["volume"] * (
        1.0 - polygon_area(body.ntheta, 1.0) / math.pi) * (1 + 1e-9)


def test_spar_stepped(tmp_path):
    profile = [(3.25, 0.0), (3.25, -4.0), (4.7, -12.0), (4.7, -40.0), (6.0, -40.0), (6.0, -45.0)]
    body = hm.spar(profile, panel_size=1.5, waterline_grading=0.3, edge_grading=0.5)
    h, _, _ = write_and_read(body, tmp_path)
    k = polygon_area(body.ntheta, 1.0)
    frustum = 8.0 / 3.0 * (3.25 ** 2 + 3.25 * 4.7 + 4.7 ** 2)
    exact_polygon = k * (3.25 ** 2 * 4.0 + frustum + 4.7 ** 2 * 28.0 + 6.0 ** 2 * 5.0)
    assert h["volume"] == pytest.approx(exact_polygon, rel=1e-10)
    assert h["waterplane"] == pytest.approx(k * 3.25 ** 2, rel=1e-10)


def test_heave_plate(tmp_path):
    body = hm.heave_plate(10.0, 9.9, 15.0, panel_size=1.5, edge_grading=1.0)
    h, _, _ = write_and_read(body, tmp_path)
    n = body.ntheta
    assert h["volume"] == pytest.approx(polygon_area(n, 10.0) * 9.9, rel=1e-10)
    sheet = sum(_area(body.mesh.points(p)) for p in body.mesh.dipoles)
    assert sheet == pytest.approx(body.exact["sheet_area_polygon"], rel=1e-10)
    with open(os.path.join(str(tmp_path), "HullMesh.pnl")) as f:
        text = f.read().splitlines()
    assert len(text[3].split()) == 5 and int(text[3].split()[4]) == len(body.mesh.dipoles)
    assert hm.DIPOLE_HEADER in text[6 + len(body.mesh.nodes) + 3 + len(body.mesh.panels)]
    # The annulus shares the nodes of the column's bottom ring
    hull_nodes = {i for p in body.mesh.panels for i in p}
    sheet_nodes = {i for p in body.mesh.dipoles for i in p}
    assert len(hull_nodes & sheet_nodes) == n


def _area(points):
    p = np.array(points, float)
    s = np.cross(p[1] - p[0], p[2] - p[0])
    if len(p) == 4:
        s = s + np.cross(p[2] - p[0], p[3] - p[0])
    return 0.5 * np.linalg.norm(s)


def test_disc_and_plate_on_a_hull(tmp_path):
    """Sheets are added to a hull body; they add no volume, so the hydrostatics stay those of the hull."""
    hull = hm.cylinder(3.0, 6.0, panel_size=0.5)
    alone, _, _ = write_and_read(hull, tmp_path / "hull")
    sheet_disc = hm.disc(2.0, 12.0, ntheta=32, nrings=12, rim_grading=1.0)
    sheet_plate = hm.plate(4.0, 2.0, 10.0, panel_size=0.25, rim_grading=0.5, centre=(8.0, 0.0))
    assert sum(_area(sheet_disc.mesh.points(p)) for p in sheet_disc.mesh.dipoles) == pytest.approx(
        polygon_area(32, 2.0), rel=1e-10)
    assert sum(_area(sheet_plate.mesh.points(p)) for p in sheet_plate.mesh.dipoles) == pytest.approx(8.0, rel=1e-12)
    body = hm.add_sheets(hm.cylinder(3.0, 6.0, panel_size=0.5), sheet_disc, sheet_plate)
    assert len(body.mesh.dipoles) == 384 + 16 * 8 and len(body.mesh.panels) == len(hull.mesh.panels)
    h, _, _ = write_and_read(body, tmp_path / "both")
    assert np.allclose(h["mass_matrix"], alone["mass_matrix"]) and np.allclose(h["restoring"], alone["restoring"])
    xs = [body.mesh.nodes[i - 1][0] for p in body.mesh.dipoles[384:] for i in p]
    assert min(xs) == pytest.approx(6.0) and max(xs) == pytest.approx(10.0)
    h, _, _ = write_and_read(body, tmp_path / "mass", mass=5e4, cog=(0.0, 0.0, -4.0), gyration=(2.0, 2.0, 1.0))
    assert h["mass"] == 5e4 and h["restoring"][2, 2] == pytest.approx(alone["restoring"][2, 2])


def test_sheet_alone_is_refused(tmp_path):
    for sheet in (hm.disc(1.0, 50.0, panel_size=0.2), hm.plate(4.0, 2.0, 10.0, panel_size=0.5)):
        with pytest.raises(ValueError, match="NELEM = 0"):
            hm.write_input(sheet, str(tmp_path))
        assert not os.listdir(str(tmp_path))


def test_dipole_checks():
    with pytest.raises(ValueError, match="panel size"):
        hm.check_dipoles(hm.disc(5.0, 0.5, panel_size=1.0).mesh)
    # A sheet over the column's own bottom disc overlaps coplanar hull panels
    body = hm.cylinder(5.0, 6.0, panel_size=1.0)
    hm.revolve(body.mesh, [(4.0, -6.0), (2.0, -6.0)], body.ntheta, dipole=True)
    with pytest.raises(ValueError, match="overlaps|coincides"):
        hm.check_dipoles(body.mesh)
    hm.check_dipoles(hm.heave_plate(5.0, 6.0, 8.0, panel_size=1.0).mesh)


def test_hull_checks():
    body = hm.cylinder(5.0, 6.0, panel_size=1.0)
    body.mesh.panels[3].reverse()
    with pytest.raises(ValueError, match="same direction"):
        hm.check_hull(body.mesh)
    body = hm.box(4.0, 4.0, 2.0, panel_size=1.0)
    body.mesh.panels.pop()
    with pytest.raises(ValueError, match="not closed"):
        hm.check_hull(body.mesh)
    body = hm.box(4.0, 4.0, 2.0, panel_size=1.0)
    body.mesh.panels = [p[::-1] for p in body.mesh.panels]
    with pytest.raises(ValueError, match="into the fluid"):
        hm.check_hull(body.mesh)


def test_lid(tmp_path):
    for body in (hm.cylinder(5.0, 10.0, panel_size=1.0), hm.box(20.0, 10.0, 5.0, panel_size=1.5),
                 hm.hemisphere(5.0, panel_size=0.8)):
        h, lines, warnings = write_and_read(body, tmp_path / body.name, irsp=True)
        nodes, panels, _, _ = hm.read_pnl(os.path.join(str(tmp_path / body.name), "WaterPlaneMesh.pnl"))
        assert all(z == 0.0 for _, _, z in nodes)
        area = sum(_area([nodes[i - 1] for i in p]) for p in panels)
        assert area == pytest.approx(h["waterplane"], rel=1e-10)
        assert not warnings and lines[-1].split()[-1] == "1"
        mesh = hm.Mesh()
        mesh.nodes, mesh.panels = nodes, panels
        assert all(hm.hams_normal(mesh.points(p))[2] > 0.0 for p in panels)


def test_hydrostatics_centre_and_mass():
    body = hm.box(20.0, 10.0, 5.0, panel_size=2.0)
    h = hm.hydrostatics(body.mesh, (0.0, 0.0, -2.0))
    m = 1e6
    assert h["mass_matrix"][0, 4] == pytest.approx(m * (-2.5 + 2.0), rel=1e-10)
    assert h["restoring"][3, 3] == pytest.approx(RG * 20.0 * 10.0 ** 3 / 12.0, rel=1e-10)
    h = hm.hydrostatics(body.mesh, (1.0, 0.0, 0.0), mass=8e5, cog=(0.0, 0.0, -1.0), gyration=(3.0, 6.0, 7.0))
    assert h["mass_matrix"][3, 3] == pytest.approx(8e5 * (9.0 + 1.0), rel=1e-12)
    assert h["mass_matrix"][4, 4] == pytest.approx(8e5 * (36.0 + 1.0 + 1.0), rel=1e-12)
    assert h["mass_matrix"][1, 5] == pytest.approx(8e5 * -1.0, rel=1e-12)
    # C44 = rho g (I_xx + V z_b) - m g z_g, z_b = -2.5, z_g = -1
    assert h["restoring"][3, 3] == pytest.approx(RG * (20.0 * 1000.0 / 12.0 - 1000.0 * 2.5) + 8e5 * hm.G, rel=1e-10)
    # C35 = -rho g int (x - xr) dA = rho g A xr; C46 = -rho g V (xb - xr) + m g (xg - xr)
    assert h["restoring"][2, 4] == pytest.approx(RG * 200.0, rel=1e-10)
    assert h["restoring"][3, 5] == pytest.approx(RG * 1000.0 - 8e5 * hm.G, rel=1e-10)


def test_multi_body(tmp_path):
    body = hm.cylinder(2.0, 4.0, panel_size=0.5)
    positions = [(0.0, 0.0, 0.0), (20.0, 0.0, 0.0), (10.0, 15.0, math.radians(30.0))]
    results, lines, warnings = hm.write_input(body, str(tmp_path), positions, irsp=True)
    assert not warnings
    for i, (x, y, a) in enumerate(positions, start=1):
        for name in ("HullMesh_%d.pnl", "WaterPlaneMesh_%d.pnl", "Hydrostatic_%d.in"):
            assert os.path.exists(os.path.join(str(tmp_path), name % i))
        nodes, panels, _, _ = hm.read_pnl(os.path.join(str(tmp_path), "HullMesh_%d.pnl" % i))
        assert np.allclose(nodes, body.mesh.nodes, atol=1e-9)          # local frame, placed by the LCS line
        cog, blocks = hm.read_hydrostatic(os.path.join(str(tmp_path), "Hydrostatic_%d.in" % i))
        assert cog[:2] == pytest.approx([x, y], abs=1e-9)
        assert blocks[3][2, 2] == pytest.approx(results[0]["restoring"][2, 2], rel=1e-5)
        assert blocks[3][3, 3] == pytest.approx(results[0]["restoring"][3, 3], rel=1e-5)
    # The lines in the fixed formats of ReadControlFile
    assert int(lines[1][24:40]) == 3
    lcs = [[float(l[26 + 12 * k:38 + 12 * k]) for k in range(4)] for l in lines[2:5]]
    assert lcs[2] == pytest.approx([10.0, 15.0, 0.0, 0.524])
    xr = [[float(l[28 + 12 * k:40 + 12 * k]) for k in range(3)] for l in lines[6:9]]
    assert xr[1] == pytest.approx([20.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="one body"):
        hm.write_input(hm.heave_plate(5.0, 6.0, 8.0, panel_size=1.0), str(tmp_path / "x"), positions[:2])


def test_command_line(tmp_path, capsys):
    hm.main(["heave-plate", "--column-radius", "5", "--draft", "6", "--plate-radius", "8", "--panel-size", "1",
             "--edge-grading", "1", "--out-dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert "dipole panels" in out and "ISOLV = 1" in out
    assert os.path.exists(os.path.join(str(tmp_path), "HullMesh.pnl"))
    hm.main(["cylinder", "--radius", "3", "--draft", "6", "--panel-size", "0.5", "--disc", "2", "12",
             "--plate", "4", "2", "10", "8", "0", "--out-dir", str(tmp_path / "sheets")])
    _, panels, dipoles, _ = hm.read_pnl(os.path.join(str(tmp_path / "sheets"), "HullMesh.pnl"))
    assert panels and len(dipoles) > 0
    with pytest.raises(SystemExit):
        hm.main(["disc", "--radius", "1", "--depth", "50"])
