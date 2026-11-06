# Mean wave drift

HAMS-MREL computes the mean (time-averaged) second-order wave loads on bodies in regular waves, from the first-order solution. It has two methods:

- **Far-field** (momentum conservation): F_x, F_y and M_z of the whole system, from the far-field amplitude of the waves the bodies send out.
- **Near-field** (pressure integration over the mean wetted surface): all six loads of each body, as the sum of six terms (waterline, flare, velocity, displacement, rotation, hydrostatic).

With `mean_drift = 2`, both give the loads of every pair of headings at the same frequency, for bidirectional seas.

This page covers usage, the choice of method, mesh guidance, limits, and the evidence behind them. The design is in [the spec](superpowers/specs/2026-11-04-second-order-waves-design.md), decisions D1 to D3.

## Usage

Put a `SecondOrder.in` in the input directory. Without it, HAMS-MREL runs exactly as before. All keys and the checks on them are in [Input Files](input_files.md#secondorderin-optional). A typical file:

```
! Mean drift with both methods, default velocity and waterline methods
mean_drift      = 1
drift_methods   = far near
motion_state    = free
```

Requirements: `ISOL = 1` in `ControlFile.in`, and depth `H` at most 500 m (use `H < 0` for deep water). The near-field method does not run with generalized mode files.

Outputs, in `Wamit_format/`, are described in [Output Files](output_files.md):

- `MeanDriftMomentum.8`: far-field F_x, F_y, M_z of the system, M_z about the origin.
- `MeanDriftPressure.9`: near-field loads of each body, moments about that body's `XR`.
- `MeanDriftTerms.txt`: the six near-field terms (with `write_terms` 1 or 2).

When both methods run, the log prints their relative difference for F_x, F_y and M_z (M_z moved to the origin). Use it as a mesh check: a large difference means the near-field result has not converged.

## Methods and defaults

| Key | Default | Other value |
| :--- | :--- | :--- |
| `velocity_method` | `differencing`: tangential velocity from the potential on neighbouring panels | `source`: velocity from an extra source solve |
| `waterline_method` | `extrapolate`: waterline potential from the adjacent panel and its gradient | `field`: boundary integral formula at a point just off the waterline |

The defaults were kept after the convergence study below. `differencing` meets the gate of the plan (near-field and far-field agree and converge on a free body), and on rounded hulls it is closer to far-field than `source` at every frequency tested from 0.60 rad/s up. There is no switch of method by frequency.

## Which method to use

- **Rounded hulls** (no sharp edges at the bilge or corners): use the defaults. On the rounded cylinder at 7680 panels, `differencing` is within 0.24% of far-field from 0.70 to 1.80 rad/s, free, and within 0.08% fixed. `source` is off by 0.12% to 2.3% on the same mesh.
- **Hulls with sharp edges**, such as a flat bottom meeting a vertical wall:
  - Near-field with `differencing` converges slowly. On the sharp cylinder at 7168 panels it is still +17.4% at 0.95 rad/s, free.
  - On that cylinder `source` is closer up to about 1.2 rad/s (ka about 1.5), and `differencing` is closer from 1.4 rad/s (ka about 2.0). This crossing rests on one radius. On the sharp box `source` is closer at all three frequencies run, 0.80, 1.00 and 1.40 rad/s. ka is inferred from the deep-water dispersion relation, a = 10 m.
  - `source` costs about twice the run time at 7680 panels: +98% (ifx, 32 threads). The earlier record is +51% (ifort, 8 threads).
  - `source` cannot remove irregular frequencies: `IRSP = 1` is rejected. At the two irregular frequencies tested it is off by 32% to 198% (see below). Keep away from them, or use `differencing` with a lid there.
- **Free sharp hull at low frequency:** use far-field. At 0.70 rad/s, where the net force is small next to its terms, `differencing` has the wrong sign at 7168 panels and `source` is still 10.3% off.
- **Sharp box:** near-field F_y and M_z are not reliable, with either method, on the meshes tested. F_y is +12.5% off with `source` at 7040 panels. Round the vertical corners and the bilge if near-field F_y or M_z are needed. That makes it a different hull: judge it against its own far-field.
- **Waterline:** keep `extrapolate`. `field` is opt-in. It is first order where `extrapolate` is second order, and in finite depth it can jump with the mesh or the offset (see [Limits](#limits)).
- **Settings:** `GRN_KIND = 8` changed no printed digit of F_x. `ISOLV = 2` (GMRES) matches LU to within 0.003 in F_x / (ρg) on the cases tried.

For F_x, F_y and M_z of the whole system, far-field is the reference on any hull. Near-field is needed for loads per body and for F_z, M_x and M_y.

## Mesh guidance

- Refine towards the waterline. The cylinders here use wall rows cosine spaced towards the waterline, so the top rows are thin. Since the free-surface image term is integrated over the panel (see [Developer Documentation](developer_documentation.md#free-surface-image-term)), near-field converges on such rows: on the rounded cylinder on springs at 1.10 rad/s it is +1.36%, +0.40%, +0.21% from far-field at 480, 1920, 4320 panels.
- Refine at corners and at the bilge. Better, round them if the real hull allows. On the box, rounding the vertical corners and the bilge took near-field against far-field from 23.4% (sharp, 7040 panels) to 0.41% (rounded, 8000 panels), for the force vector at 1.00 rad/s.
- Check convergence with both methods on two meshes, panels per side doubled. The rate of fall per doubling is a better check than one difference.
- Irregular frequencies: with `differencing`, use `IRSP = 1` with a lid near them. A lid with thin rim panels did not shift results on the cylinder: over a scan of the rim ratio down to 0.17, `IRSP = 1` moved near-field F_x by at most 0.42% against `IRSP = 0`. This is one cylinder at 0.95 and 1.10 rad/s only.

## Limits

- **Generalized modes:** the near-field method is rejected when a `gen_mod` file is present. Far-field runs, but no integration test covers it.
- **Narrow gaps:** not tested. In the test fixtures, bodies are at least 20 m apart centre to centre. With `waterline_method = field`, no other body may be closer to the field point than `waterline_offset` times the segment length.
- **Multi-body:** far-field gives only the total of all bodies. Near-field per body is checked against the far-field total in the two-body fixtures.
- **`waterline_method = field` in finite depth:** term I can jump with the mesh or the offset. On 32 x 16 rows at depth 10 and offset 0.01, term I is -11.40% at 1.0 rad/s, against -1.93% and -2.06% at offsets 0.005 and 0.05 (measured on this branch, see [the fixtures](../test/data/drift_cyl/README.md)). A separate measurement on commit 36d1733, not reproduced on this branch, replaced the Green function with a direct eigenfunction sum: term I went from -11.40% to -1.96%. In that measurement, first-order results of the repository cases moved by at most 0.25%, and small panels at kh of 50 or more by up to 9%. A fix (Linton's method for small R/h) is planned on a separate branch and is not merged. The run logs a warning when this method is used in finite depth.
- **Not tested:** a correction of the velocity fit for the normal offset of neighbouring panels. It was skipped because the term II offset it was meant to explain went away with the image term fix.
- **Second-order QTFs** (difference and sum frequency) are not implemented. Only the mean drift is.

## Convergence study

Run in November 2026 on Snellius (`fat_genoa`, 32 threads, ifx 2024.2) with [drift_convergence.py](../pre_postprocessing_python/SecondOrder/drift_convergence.py). 80 HAMS runs and 12 Nemoh runs. The per-frequency tables are not in the repository. All values below are measured unless marked inferred.

**Hulls.** Truncated cylinder, radius 10 m, draft 10 m, depth 100 m, with a sharp bilge (graded rows) or a bilge radius of 3 m (cosine rows); 448 to 7168 and 480 to 7680 panels. Bottom-mounted cylinder, radius 5 m, depth 10 m, against MacCamy-Fuchs; 256 to 4096 panels. Box 20 x 10 x 5 m, depth 50 m, sharp, with rounded vertical corners only, and with rounded corners and bilge. Frequencies 0.40 to 1.80 rad/s in steps of 0.2, plus 0.70, 0.95 and 1.10, set before any result. 0.80 rad/s is the heave resonance, and 1.60 is 3% from the first irregular frequency, 1.548 rad/s.

### Near-field against far-field, truncated cylinder

F_x, heading 0, finest mesh, percent of far-field, `differencing` / `source`:

| rad/s | Sharp fixed, 7168 | Sharp free, 7168 | Rounded fixed, 7680 | Rounded free, 7680 |
| ---: | ---: | ---: | ---: | ---: |
| 0.70 | +4.76 / -0.06 | -136 / +10.3 | -0.08 / -0.41 | +0.13 / +1.66 |
| 0.80 | +3.52 / -0.26 | +141 / +5.27 | -0.08 / -0.47 | +0.24 / -1.23 |
| 0.95 | +2.08 / -0.34 | +17.4 / +0.53 | -0.08 / -0.45 | +0.03 / +0.12 |
| 1.00 | +1.74 / -0.34 | +11.3 / +0.05 | -0.07 / -0.45 | +0.04 / -0.33 |
| 1.10 | +1.29 / -0.44 | +5.29 / -0.58 | -0.06 / -0.55 | +0.03 / -0.95 |
| 1.20 | +1.00 / -0.69 | +2.68 / -0.97 | -0.04 / -0.78 | +0.03 / -1.26 |
| 1.40 | +0.53 / -1.59 | +0.89 / -1.84 | -0.02 / -1.45 | +0.03 / -1.86 |
| 1.60 | +0.31 / +2.34 | +0.41 / +3.13 | +0.01 / +1.52 | +0.04 / +2.31 |
| 1.80 | +0.41 / -1.94 | +0.46 / -1.63 | +0.08 / -1.60 | +0.10 / -1.37 |

- The crossing on the sharp hull, `source` closer to 1.20 rad/s and `differencing` from 1.40, also holds at 1792 panels.
- At 0.40 and 0.60 rad/s the free net force is about 3e-6 and 2e-3 of ρg, so percentages there say little. At 0.40 fixed, rounded, the two methods nearly tie (-0.060 / -0.044).
- At 1.60 rad/s `source` is off by 1.5% to 3.1%, with the sign opposite to its neighbours. `differencing` is within 0.5%.

**Sharp free hull at 0.70 rad/s**, F_x / (ρg):

| Panels | Far | Near, `differencing` | Near, `source` |
| ---: | ---: | ---: | ---: |
| 448 | 0.1193 | -0.2829 | 0.1325 |
| 1792 | 0.1127 | -0.1327 | 0.1281 |
| 7168 | 0.1107 | -0.0402 | 0.1222 |

`differencing` keeps falling but has the wrong sign at 7168 panels. `source` is 10% to 14% off on all meshes. Its error, 0.0115 at 7168, is 0.44% of term I.

**Rounded hull, `source`.** Its error falls with refinement, but reaching zero is not shown. Fixed at 0.95 rad/s, near minus far is -0.0478, -0.0369, -0.0290, -0.0238 of ρg at 480, 1920, 4320, 7680 panels. The ratio per doubling goes from 0.77 to 0.65, so the local order rises from about 0.4 to 0.6. A constant offset would make it fall instead. At order 0.6, halving the error needs about 10 times the panels (inferred).

**Bottom-mounted cylinder at 4096 panels, against MacCamy-Fuchs:** `differencing` within 0.12% at all 11 frequencies, far-field within 0.17%. `source` within 0.15% up to 1.20 rad/s, then -0.30%, -0.73%, -1.58% at 1.40, 1.60, 1.80.

### Box, heading 30, 1.00 rad/s, fixed

Near-field against far-field, percent:

| Hull | Panels | Method | F_x | F_y | M_z | Force vector |
| :--- | ---: | :--- | ---: | ---: | ---: | ---: |
| Sharp | 7040 | `differencing` | +2.8 | +33.9 | +20.0 | 23.4 |
| Sharp | 7040 | `source` | -0.09 | +12.5 | +9.3 | 8.6 |
| Rounded vertical corners | 1760 | `differencing` | +9.7 | +35.8 | -4.1 | 25.5 |
| Rounded corners and bilge | 8000 | `differencing` | -0.13 | -0.58 | +0.08 | 0.41 |
| Rounded corners and bilge | 8000 | `source` | -0.04 | +2.1 | +1.14 | 1.4 |

- The M_z excess follows the vertical corners and the F_y excess the bilge, on two meshes each.
- On the sharp box, `source` removes about 60% of the residual, which still falls at about 0.7 per doubling.
- Free sharp box at 1760 panels: F_x is -210% with `differencing` and -55% with `source`. Free rounded box at 2000 panels: +6.0% and -9.8%. No free box above 2000 panels was run.

### Irregular frequencies

Free cylinders, F_x error against far-field, percent. `differencing` with `IRSP` 0 / 1, and `source` with `IRSP = 0`:

| Hull, panels | 1.548 rad/s, diff | 1.548, src | 1.939 rad/s, diff | 1.939, src |
| :--- | ---: | ---: | ---: | ---: |
| Rounded, 1920 | +4.0 / +0.50 | -46.6 | +4.7 / +1.2 | -32.4 |
| Sharp, 1792 | +15.5 / +1.6 | -125.9 | +12.6 / +2.7 | -198.5 |

At 1.10 rad/s, away from them, `IRSP = 1` changes the error by at most 0.03 points. On the rounded hull at 1.548 rad/s, `source` goes from -66.1% at 480 panels to -46.6% at 1920.

### Waterline method

Free cylinders at 1920 and 1792 panels, F_x error against far-field. Rounded hull: `extrapolate` is closer at all 11 frequencies, by 0.14 to 1.12 points from 0.70 rad/s up. Sharp hull: `field` is closer by 0.27 to 0.57 points from 0.95 to 1.80 rad/s, where both carry the sharp-edge error (1.3% to 28% with `extrapolate`).

### Settings

- `GRN_KIND = 8`: identical F_x to all printed digits at all 15 entries, at the same time per frequency, and 0.9 GB more memory at 4320 panels.
- `ISOLV = 2`: F_x differs from LU by at most 0.0009 (near) and 0.0003 (far) of ρg on the rounded hull, and 0.003 on the sharp hull at 0.70 rad/s. Run time within 3% at 4320 panels.

### Run time

On Snellius, 32 threads, 11 frequencies, wall time per run:

| Hull, panels | `differencing` | `source` |
| :--- | ---: | ---: |
| Sharp, 7168 | 884 to 885 s | 1750 s |
| Rounded, 7680 | 1031 s | 2042 to 2045 s |

The largest run, the rounded box with `source` at 8000 panels, used 8.40 GB.

## Comparison with Nemoh

Nemoh v3.0.3 ran as a separate program, built outside this repository with ifort and sequential MKL, on Linux (Snellius) and, for the coarsest meshes, on Windows. No Nemoh source or binary, and no Nemoh output file, is in the repository. Only its documented input and output files were used. Nemoh is a code comparison, not a reference truth: it also uses constant panels and shares error sources with the near-field method.

**Inputs.** The nodes and panels of the HAMS case, written by `make_mesh.py --format nemoh`. No lid, `IRSP = 0` in HAMS, no external loads. Fixed bodies through a stiffness of 1e13 in all six modes. 21 frequencies from 0.40 to 1.40 rad/s, below the first irregular frequency. Heading 0, and headings 0 and 30 for the pair. QTF part with the DUOK term only. Drift is read from the diagonal of `results/QTF/OUT_QTFM_N.dat`, times ρg.

**Rule.** The plan set 5% on a converged F_x. Nemoh did not converge on the meshes that worked (below), so the rule for an unconverged mesh applies: HAMS no further from Nemoh than the Nemoh-HydroSTAR spread. That spread is a median of 9.7% in F_x on Nemoh's own `8b_QTF_Cylinder` case at 634 panels, for ω ≥ 0.9 rad/s.

**Nemoh at 4320 panels is unusable.** Its F_x on the rounded fixed cylinder is -3.85 of ρg at 0.95 rad/s, against 5.35 at 1920 panels and 5.31 from HAMS. Nemoh gave no warning and returned 0. On the same mesh its first-order results agree with HAMS (table below), and so do its fixed F_z and M_y, within 1.2%. The cause is unknown. A guess, untested: the top row is 0.0196 m tall at 4320 panels, against 0.044 m at 1920. The comparison therefore rests on 448 to 1920 panels.

**First order**, worst absolute difference in A, B, X and RAOs from 0.70 to 1.10 rad/s, percent:

| Hull, panels | A | B | X | RAO |
| :--- | ---: | ---: | ---: | ---: |
| Rounded 480 | 3.99 | 5.86 | 1.96 | 3.77 |
| Rounded 1920 | 2.27 | 3.47 | 1.01 | 2.11 |
| Rounded 4320 | 1.54 | 2.15 | 0.79 | 1.44 |
| Sharp 448 | 3.22 | 4.84 | 1.67 | 10.07 |
| Sharp 1792 | 1.84 | 3.46 | 0.84 | 6.25 |

The difference falls with the mesh in all columns. The largest RAO difference on the sharp hull is in surge at 0.70 rad/s.

**Drift**, F_x, heading 0, absolute difference between HAMS near-field and Nemoh, frequencies up to 1.40 rad/s with far-field above 0.5 of ρg:

| Hull, motion | Panels | `differencing` median / worst | `source` median / worst | Spread rule |
| :--- | ---: | ---: | ---: | :--- |
| Rounded, fixed | 1920 | 1.80 / 4.26 | 1.32 / 3.24 | met by both |
| Rounded, free | 1920 | 1.00 / 39 | 0.88 / 30 | median met; worst at 0.80 rad/s |
| Sharp, fixed | 1792 | 3.93 / 7.18 | 1.16 / 1.99 | met by both |
| Sharp, free | 1792 | 14.43 / 176 | 0.80 / 2.48 | missed by `differencing`, met by `source` |

- On the rounded free hull the worst value is at the heave resonance, 0.80 rad/s, where Nemoh itself moves from 0.86 to 1.45 of ρg between 480 and 1920 panels.
- The sharp free hull matches the near-against-far result: `differencing` carries the sharp-edge error, `source` does not.
- Bidirectional pair, rounded free, headings 0 and 30, 1920 panels: within 6.7% at 0.70 rad/s and within 0.85% at 0.95 and 1.10. Informative only. Nemoh's documentation states its QTFs are untested for multi-body problems, and no pair ran at 4320 panels.

## Tests

The integration fixtures and their measured values are described in [test/data/drift_cyl/README.md](../test/data/drift_cyl/README.md). Two limits have thin margins on ifx and Linux (Snellius, ifx 2024.2, measured):

| Test | Measured, ifx / Linux | Measured, ifort / Windows | Limit |
| :--- | ---: | ---: | ---: |
| `source_two_body_total` | 4.69% | 4.69% | 5% |
| `source_near_far_low_frequency` | 4.10% | 4.11% | 5% |

A change in compiler or flags may trip them. Check the printed value before changing a limit.

## Compilers

ifx 2024.2 with `-heap-arrays 0` stopped the far-field method with a segmentation fault on the first frequency. The cause was an array constructor `[s1, s2]` of a type with allocatable components in `far_field_pair`. A 40-line program with no HAMS code fails the same way under ifx 2024.2 with `-heap-arrays 0`, and runs under ifort with that flag and under ifx 2025.1.1. The code now copies the two states into a named array (commit dff7a32). With that change, unit (133) and integration (118) tests all pass under ifx 2024.2 on Linux with `-qmkl -qopenmp -heap-arrays 0` (Snellius job 28933784), with the same two thin margins as above. No other array constructor of a type with allocatable components is in `src/` or `test/`. ifx 2025.0, used by the Linux CI, was not tested here.

## Open items

- A free sharp body at low frequency has no good near-field result. Far-field is the number to use.
- Nemoh at 4320 panels: the cause is unknown, and testing it needs another Nemoh run.
- WAMIT remains an option if a licence becomes available. It is still needed to settle the conjugation and axes of the WAMIT-format output.
