# Roland et al. (2026) benchmark

HAMS against the second-order benchmark of Roland, Pribadi, Al-Ghuwaidi and Lataire (2026), J. Phys.: Conf. Ser. 3346,
012006, doi:10.1088/1742-6596/3346/1/012006 (open access, CC BY 4.0). Results: `docs/second_order_validation.md`, sections
"OC4 against the solvers of Roland et al. (2026)" and "VolturnUS-S with the free-surface integral".

Not in the repository: the paper's PDF, the WAMIT files (OpenFAST r-test `marin_semi.12d`; IEA-15-240-RWT v1.1.18
HydroData, Apache-2.0) and the case inputs derived from the IEA mesh. Run Python on them with `-I`.

| File | What it does |
|---|---|
| `extract_roland_curves.py PDF OUT_CSV` | Reads the curves of Figs. 5 and 7 from the vector paths of the PDF |
| `roland2026_fig5_fig7.csv` | Its output: figure, platform, diagonal, solver, omega_1, |QTF| in N/m^2 (derived from the paper, CC BY 4.0) |
| `check_wamit_extraction.py CSV FIG WAMIT_12D` | The WAMIT curves against the WAMIT file (0.034 % for Fig. 5, 0.042 % for Fig. 7) |
| `compare_diagonals.py CSV FIG OUT_PNG OUT_TXT LABEL=HAMS_FILE ...` | HAMS on the two diagonals against the six solver set-ups of the paper |
| `compare_maps.py HAMS_FILE WAMIT_12D OUT_PNG OUT_TXT TITLE` | The whole QTF, magnitude and phase, in the style of Figs. 4 and 6 |
| `make_volturnus_case.py IEA_HYDRODATA_DIR OUT_CASE_DIR` | The VolturnUS-S HAMS case with the free-surface integral |
| `plot_volturnus_case.py`, `plot_cost.py` | Figures of the case meshes and of the run cost |
| `snellius/` | Job file, input checksums, cost and the cluster steps |

OC4 commands (HAMS outputs of the full WAMIT hull at 5022 and 9208 panels):

```
python -I extract_roland_curves.py roland2026.pdf roland2026_fig5_fig7.csv
python -I check_wamit_extraction.py roland2026_fig5_fig7.csv 5 marin_semi.12d
python -I compare_diagonals.py roland2026_fig5_fig7.csv 5 oc4_fig5_hams.png oc4_fig5_hams.txt \
    "HAMS, WAMIT hull, 5022 panels=<L0>/output/Wamit_format/QTF_Indirect.11d" \
    "HAMS, WAMIT hull, 9208 panels=<L1>/output/Wamit_format/QTF_Indirect.11d"
python -I compare_maps.py <L1>/output/Wamit_format/QTF_Indirect.11d marin_semi.12d oc4_fig4_hams_maps.png \
    oc4_fig4_hams_maps.txt "DeepCwind ..." --band 0.25 0.35
```
