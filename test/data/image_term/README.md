# Image term fixtures

Used by `test/test_integration_image.f90`. All cases are a cylinder of radius 5 m and draft 5 m with a flat bottom, in 50 m of water, at 1.35 rad/s and heading 0, with no external loads. Inputs are committed as they are: they were generated with the cylinder and control file writers of the mean drift work, which are not on this branch.

Rows are at z = -(d/2)(1 - cos(pi i/nz)) with cosine rows, thin at the waterline, or z = -d i/nz with uniform rows. Ring radii on the bottom follow the same choice.

| Case | Rows | ntheta x nz x nr | Panels | Other |
| :--- | :--- | :--- | ---: | :--- |
| `image_cos` | cosine | 20 x 6 x 3 | 180 | every panel centroid is a field point |
| `image_cos_half` | cosine | 20 x 6 x 3 | 90 | half mesh, `ISY = 1` |
| `image_cos_gmres` | cosine | 20 x 6 x 3 | 180 | `Solver_type = 2` |
| `image_cos_lid` | cosine | 20 x 6 x 3 | 180 | `IRSP = 1`, uniform lid 20 x 4, every panel centroid is a field point |
| `image_cos_two` | cosine | 20 x 6 x 3 | 180 each | a second body 5 km away, every panel centroid of the first is a field point |
| `image_uni` | uniform | 20 x 6 x 3 | 180 | |
| `image_cos_fine` | cosine | 30 x 9 x 4 | 390 | field points on the top 5 rows of the first column |
| `image_uni_fine` | uniform | 30 x 9 x 4 | 390 | |

## Checks, before and after the change

Before is `c40ad4d`. All are relative.

| Test | Limit | Before | After |
| :--- | ---: | ---: | ---: |
| `top_rows_smooth` | 0.2 | 0.43, fails | 0.15 |
| `lid_matches_no_lid` | 5e-3 | 8.9e-3, fails | 3.3e-3 |
| `field_points_match_collocation` | 1e-5 | 1.1e-7 | 9.1e-8 |
| `field_points_match_collocation_multi` | 1e-5 | 2.2e-7 | 2.9e-7 |
| `half_matches_full` | 1e-4 | 0 | 0 |
| `gmres_matches_direct` | 1e-3 | 1.2e-5 | 1.7e-5 |
| `multi_matches_single` | 1e-3 | 7.8e-4 | 7.8e-4 |
| `haskind_heave` | 1e-2 | 7.8e-3 | 9.3e-3 |

`multi_matches_single` is limited by the far body, which still radiates into this one. `cosine_and_uniform_approach` only asks that cosine and uniform rows get closer with refinement. It passes before and after, because A, B and X on this hull are dominated by the sharp bilge: A55 differs by 5.9% and 5.4% at the two meshes on both.
