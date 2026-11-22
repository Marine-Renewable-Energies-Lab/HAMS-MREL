# Nemoh-paper cases on Snellius (cluster-side steps)

Placeholders: `<user>` (cluster login), `<login_host>` (login node), `<account>` (Slurm account), `<run_root>` (the run
folder), `<nemoh_root>` (where Nemoh v3.0.3 is built). The filled-in copies of the job files are kept outside the
repository. Submit only after the user's approval of the cost (COST.md). Login by key only; in Git Bash every ssh/scp
line needs the `MSYS_NO_PATHCONV=1` prefix.

## 1. Inputs (laptop)

```
python -I make_cases.py <Nemoh-3.0.3>/TestCases <work>/cases
python -I cost_estimate.py <work>/cases/cases.json <work>/snellius/cost.png <work>/snellius/COST_table.md
```

`<work>` is outside the repository: the meshes are derived from Nemoh's test cases (GPL). Then one checksum list per case:
`(cd <work>/cases/<case>/input && md5sum *) > md5sums_<case>.txt`.

## 2. Staging (copy only)

- `<run_root>/run_nemohpaper/<case>/input/` for hp_m1, hp_m2, sw_fo, cyl8b_diff, sw_diff, sw_full, oc4_bi,
  oc3_diff, and `<run_root>/run_nemohpaper/nemoh_hp_m1/` (Nemoh's input of the heave-plate m1 mesh).
- In `<run_root>/run_nemohpaper/`: `common.sh`, the five `run_J*.sbatch`, `mechanics_to_hydrostatic.py`, the
  `md5sums_<case>.txt` lists, and `hamsmrel`, the binary built from `development_2` at 8f27716 (md5 in
  `md5sums_bin.txt`).
- Each job first checks its inputs against the checksum lists and stops if one differs.

## 3. Submit (only after the approval)

```
cd <run_root>/run_nemohpaper && for j in run_J1_small run_J2_sw_diff run_J3_sw_full run_J4_oc4_bi run_J5_oc3_diff; do sbatch $j.sbatch; done
```

For the options of COST.md (J2b, J5b: 32 frequencies to 2.01 rad/s), cut the frequency list of `ControlFile.in` of
sw_diff or oc3_diff on the laptop, write new checksums, copy again, and set `#SBATCH -t 09:00:00` (J2b) or `15:00:00`
(J5b).

## 4. Outputs

Each case writes `log.txt`, `time.txt`, `result.json` and `output/`; Nemoh writes its `results/` in `nemoh_hp_m1`, and
`nemoh_mech_<case>` keeps the hydrosCal output used for Hydrostatic.in of cyl8b_diff and oc3_diff.
