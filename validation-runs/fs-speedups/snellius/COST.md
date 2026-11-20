# Cost of the fs-speedups Snellius jobs (G1 item 3 and the OC4 comparison)

Nothing here has been submitted. fat_genoa, account tdsei3290, 24 CPUs per job (the 24-CPU billing minimum), 32G (below
the memory that would be billed as more CPUs). Core-hours = wall time x 24 CPUs. Expected times come from the measured
times of the same cases in earlier jobs (sacct and result.json); the maximum is the time limit x 24.

| Job | Runs | Measured basis | Expected wall | Expected core-h | Time limit | Max core-h |
|---|---|---|---|---|---|---|
| build.sbatch | build of the branch, GRN_KIND = 4 | 29453290: 51 s | 1 min | 0.3 | 0:30:00 | 12 |
| g1_item3/run_sharp.sbatch (array of 4) | sharp m2 fixed and free, each with samples and with panels | stage 5, 29405797: HAMS 11436 s (fixed), 11394 s (free), from the "code 0, ... s" lines of logs/29405797_1.log and _2.log in snellius-qtf-directsum | 3.2 h per task | 304 (4 x 76) | 6:00:00 per task | 576 |
| g1_item3/run_round.sbatch (array of 2) | round m2 free, samples and panels | stage 2, 29336170_2: 768 s, before the near-waterline fit (which made sharp m2 28 % slower) | 0.27 h per task | 13 (2 x 6.5) | 1:00:00 per task | 48 |
| oc4_compare/run_oc4.sbatch | OC4 WAMIT model, sum, all items on, qtf_fs_circle = panels | 29453291: 26100 s (174 core-h) without the items; about 2.0 h with them (derived, oc4wfull-sum/profile/findings_speedups.md section 1) | 2.0 h | 48 | 9:00:00 | 216 |
| **Total** | 7 runs + build | | | **366** | | **852** |

Notes
- The sharp expectation gives the branch no credit for its speed-ups (deep-water field rows above k0 h = 10, which is
  about 1 rad/s and up in depth 100 m, and the psi reuse). They should make the runs shorter; by how much was not measured
  for these cases, so the stage-5 times are used as they are. The 6 h limit is 1.9 x the stage-5 time.
- The OC4 expectation (2.0 h) is a derived estimate, not a measurement: deep water from k0 h = 10 (2.3 h), psi reuse
  (2.1 h), the threaded fit (2.0 h), from the profile of 29453291. The circle panels save part of what is left of the circle
  samples. If the items saved nothing, the run would take the 7.25 h of 29453291 (174 core-h), inside the 9 h limit.
- SBU: stage 5 was 153.4 core-h and about 230 SBU on fat_genoa (snellius-qtf-directsum/findings.md section 10), about
  1.5 SBU per core-hour: the total is about 550 SBU expected, 1280 SBU at most.
- Cheaper G1 variant (not prepared): skip the two sharp samples runs (152 core-h) and compare the panels runs with the
  stage-5 outputs (results_wl, commit 47d5911). The difference would then mix the circle route with every change since
  47d5911 (field_deep_kh above all), so it would not isolate item 3. The job set above runs both routes on the same build.
- No Snellius run for item 2 (qtf_fs_inner = source): ruling 1, it is decided on the laptop numbers and stays off.
