# High-rank fibers

Each curve file identifies a fiber over ℚ with an ICARM database curve. The following direct specializations use the repository's base coordinate `T`.

| ICARM curve | Rank lower bound | Surface | Point on C (u, v) | T |
| --- | ---: | --- | --- | ---: |
| [273](273.gp) | 30 | 1020 | (2, −19) | 587609/377717 |
| [273](273.gp) | 30 | 1020 | (153, 1530) | −6421189/587609 |
| [302](302.gp) | 31 | 1092 | (9/2, 0) | −2261185/7523886 |
| [751](751.gp) | 31 | 1008 | Fixed surface over ℚ(T) | 370273/90650 |
| [819](819.gp) | 30 | 1092 | (−1/8, 37/16) | −500246/38623 |
| [1766](1766.gp) | 30 | 1092 | (75/2, 435) | 1147413/9716170 |

From the repository root, run `gp -q curves/819.gp` to verify one match. From this directory, run `gp -q 819.gp`. Each file loads the parent surface and includes the database equation and all its witness points.

### Coordinates

`fiber()` uses the specialized surface coordinates and `curve()` uses the database coordinates. For a database point `P`, `ellchangepointinv(P, change)` gives its coordinates on the fiber. Going back uses `ellchangepoint(Q, change)`. For example, after loading a match in GP, run

```gp
P = points()[1];
Q = ellchangepointinv(P, change);
[ellisoncurve(fiber(), Q), ellchangepoint(Q, change) == P]
\\ [1, 1]
```

When the database commentary gives a base parameter, `database_t` stores it and `repository_t(database_t)` gives `t_value`. Curve 273 also provides `second_coefficient_point`, `second_t_value`, `second_repository_t`, and `second_change` for its second specialization. Loading the file checks both.

### Search coverage

The database snapshot contains 78 curves with rank lower bound at least 30. The search considered the fixed 1008 surface and the 1020, 1092, and 1152 families. For a given coefficient point, matching used the rational solutions of `j(T) = j_database`, including the fiber at infinity, and checked the change of Weierstrass coordinates over ℚ.

For 1020 and 1092 the coefficient search covered the listed points and all points of the form `nG + D` with `|n| ≤ 10000` and rational 2-torsion `D`, including the identity. The table gives the Weierstrass model used to enumerate these points and the conversion to the repository coordinates.

| Surface | Weierstrass coefficients | G | (u, v) from an affine point (X, Y) |
| --- | --- | --- | --- |
| 1020 | `[1,1,1,-421,-3157]` | (−13, 24) | ((X − 23)/2, (Y + X − 11)/4) |
| 1092 | `[1,0,1,-137,380]` | (−6, 34) | ((X − 3)/2, (Y + 3X − 7)/2) |

Reduction modulo primes excluded impossible pairs before solving for rational `T`. For 1152 the search covered the four listed points, the four affine points with `u = 0`, and the four points at infinity. These are finite searches in the coefficient parameter.

The equations, rational witness points, and rank lower bounds come from the ICARM Elliptic Curve Rank Leaderboard, maintained under NSF grant DMS 2425401.
