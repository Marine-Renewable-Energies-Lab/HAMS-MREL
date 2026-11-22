"""Readers for the comparison of HAMS-MREL with the HydroStar reference files that ship with the Nemoh v3.0.3 test cases.

HydroStar files (TestCases/<case>/reference/hydrostar/<mesh>/), as read from their headers:
- AddedMass_ij.dat, Damping_ij.dat, f?f1st.rao / m?f1st.rao, surge.rao ... yaw.rao: '#' header lines, then one row per
  frequency: omega (rad/s), the amplitude of each heading, then the phase (degrees) of each heading. Units in the header
  (#UNIT): kg, kg m^2, N s/m ..., N/m, N m/m, m/m, deg/m. ADDED-MASS and DAMPING are real: the "phase" columns are 0.
- QTF files (QTF?? .rao, MQTF??.rao, SPG??.rao, Drift??.rao): #QTFMODE Diff or Sum, #NBHEADING, #HEADING list, #DIFF list
  (the offsets dw >= 0), then blocks per heading pair ("#Heading1 = b1", "#Heading2 = b2", or "#Heading = b" for one
  heading), each with one row per omega: omega, then (Re, Im) for each dw of #DIFF. The row (omega, dw) is the pair
  (omega, omega + dw); pairs beyond the frequency list are written as 0. Dimensional: N/m^2 (force) or N m/m^2 (moment)
  per unit amplitude pair; #COORD is the moment reference point.
HAMS-MREL files (output/Wamit_format), as output_files.md describes them: non-dimensional with rho = 1000, g = 9.80665 and
L = 1 (REFL): A/rho, B/(rho omega), X/(rho g A), QTF/(rho g A^2).
"""
import os
import re

import numpy as np

RHO_HS, G_HS = 1025.0, 9.81
RHO_H, G_H = 1000.0, 9.80665


def _header(path):
    hdr, rows = [], []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#"):
                hdr.append(ln.rstrip("\n"))
            elif ln.strip():
                rows.append(ln)
    return hdr, rows


def hs_rao(path):
    """(omega, amplitude[nfreq, nhead], phase_deg[nfreq, nhead], headings, unit)."""
    hdr, rows = _header(path)
    heads = []
    unit = ""
    for h in hdr:
        if h.startswith("#HEADING"):
            heads = [float(v) for v in h.split()[1:]]
        if h.startswith("#UNIT"):
            unit = h.split(":", 1)[1].strip()
    data = np.array([[float(v) for v in r.split()] for r in rows])
    n = len(heads)
    return data[:, 0], data[:, 1:1 + n], data[:, 1 + n:1 + 2 * n], heads, unit


def hs_qtf(path):
    """dict with mode ('Diff'/'Sum'), headings, dw (offsets), omega (rows) and blocks {(b1, b2): complex[nomega, ndw]}."""
    hdr_mode, heads, dws = None, [], []
    blocks = {}
    cur = None
    b1 = None
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("#QTFMODE"):
                hdr_mode = s.split()[1]
            elif s.startswith("#NBHEADING"):
                pass
            elif s.startswith("#HEADING"):
                heads = [float(v) for v in s.split()[1:]]
            elif s.startswith("#DIFF"):
                dws = [float(v) for v in s.split()[1:]]
            elif re.match(r"^##?Heading1", s):
                b1 = float(s.split("=")[1])
            elif s.startswith("#Heading2"):
                cur = (b1, float(s.split("=")[1]))
                blocks[cur] = []
            elif s.startswith("#Heading ") or s.startswith("#Heading="):
                b = float(s.split("=")[1])
                cur = (b, b)
                blocks[cur] = []
            elif s.startswith("#"):
                continue
            else:
                v = [float(x) for x in s.split()]
                blocks[cur].append(v)
    out = {}
    omega = None
    for k, rows in blocks.items():
        a = np.array(rows)
        omega = a[:, 0]
        out[k] = a[:, 1::2] + 1j * a[:, 2::2]
    return dict(mode=hdr_mode, headings=heads, dw=np.array(dws), omega=omega, blocks=out)


def hs_pair_matrix(q, block, grid):
    """The reference as a square matrix T[i, j] on grid (omega_i, omega_j) for omega_j >= omega_i, from the row
    (omega_i, dw = omega_j - omega_i); nan where the file has no value."""
    n = len(grid)
    t = np.full((n, n), np.nan + 0j)
    b = q["blocks"][block]
    for i, w in enumerate(grid):
        r = np.argmin(abs(q["omega"] - w))
        if abs(q["omega"][r] - w) > 1e-3:
            continue
        for j in range(i, n):
            c = np.argmin(abs(q["dw"] - (grid[j] - w)))
            if abs(q["dw"][c] - (grid[j] - w)) < 1e-3:
                t[i, j] = b[r, c]
    return t


def hams_qtf(path, mode_index, b1, b2, grid):
    """T[i, j] of a HAMS .1?d/.1?s file for mode I (1..6), headings (b1, b2), on grid; nan where absent."""
    n = len(grid)
    t = np.full((n, n), np.nan + 0j)
    a = np.loadtxt(path)
    sel = (a[:, 4] == mode_index) & (abs(a[:, 2] - b1) < 1e-6) & (abs(a[:, 3] - b2) < 1e-6)
    a = a[sel]
    gi = {round(w, 4): k for k, w in enumerate(grid)}
    for row in a:
        i, j = gi.get(round(row[0], 4)), gi.get(round(row[1], 4))
        if i is not None and j is not None:
            t[i, j] = row[7] + 1j * row[8]
    return t


def hams_first(path_dir):
    """A/rho, B/(rho omega) as dict[(i, j)] -> (omega, A, B); X and motion as dict[(beta, i)] -> (omega, complex)."""
    out = {}
    a = np.loadtxt(os.path.join(path_dir, "AmssDamp.1"))
    am = {}
    for row in a:
        am.setdefault((int(row[1]), int(row[2])), []).append((row[0], row[3], row[4]))
    out["AB"] = {k: np.array(v) for k, v in am.items()}
    for name, key in (("ExcForce.3", "X"), ("Motion.4", "RAO")):
        a = np.loadtxt(os.path.join(path_dir, name))
        d = {}
        for row in a:
            d.setdefault((round(row[1], 3), int(row[2])), []).append((row[0], row[5] + 1j * row[6]))
        out[key] = {k: (np.array([w for w, _ in v]), np.array([c for _, c in v])) for k, v in d.items()}
    return out
