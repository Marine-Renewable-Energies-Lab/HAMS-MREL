# J3C and J10: cost (nothing submitted)

fat_genoa, 24 CPUs billed per job (the minimum), 1.5 SBU per core-hour. Each job builds its own binary from the staged
source of the branch fix/field-limes-lid-weight at 62e46e2 (the direct-sum commit 9d32939, then the lid-weighting
commit; a few minutes, included). Source staged as src_62e46e2 with its md5 list.

| Job | Content | Expected wall (h) | Core-h | SBU | Limit (h) | Core-h at limit | SBU at limit |
|---|---|---:|---:|---:|---:|---:|---:|
| J3C sw_full | rerun of J3B (SOFTWIND, QTF⁻ and QTF⁺ with the free-surface integral, 32 frequencies to 2.01 rad/s) with the direct field sum above ν·h = 50 | 17.4 | 418 | 626 | 22 | 528 | 792 |
| J10 lid weight | heave plate p0 (1.09 m lid) and p_lid (0.375 m lid), test case 8a with two lids, each at α = 1, 0.1, 0.01; p0 unweighted on the new binary | 1.1 | 26 | 40 | 3 | 72 | 108 |
| **Both** | | | **444** | **666** | | **600** | **900** |

- **J3C expected wall.** J3B took 17.4 h (62 660 s). Its passes above ν·h = 50 took 25 300 s. In those passes, the pairs that now take the direct sum (R/h from 0.0005 to 0.05) are about 20 % of the Green-function calls on the SOFTWIND inner points (desk timing on the SOFTWIND hull). They cost 8.4 µs each, against 26 µs for PIDCOCK with Limes, so those passes should be no slower; the estimate keeps J3B's 17.4 h.
- **J10 expected wall.** p0 and p_lid took 535 s and 873 s on 17 threads in J9 (about 400 s and 650 s on 24). Three weights make about 3200 s, plus the unweighted p0 check, the 8a cases (seconds) and the build.
