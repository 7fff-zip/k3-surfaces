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

### Rank 32 on X1152 by changing the fibration

All three records start from the repository equation `surface(T,4,4)`. Specializing the coefficient curve at `(4,4)` fixes the K3 surface. On this surface a new rational function of `(T,x,y)` gives a different elliptic fibration. Its published base parameter is denoted by `t_value` in the curve file. The original `T` varies along the record curve.

| ICARM curve | Fibrations visited after the original | Final base parameter |
| --- | --- | ---: |
| [4239](4239.gp) | #4239 | −2868768/362209 |
| [6534](6534.gp) | #4239 | −248188/1981293 |
| [6782](6782.gp) | #4239 → #3900 → #6782 | 431951/282772 |

Each arrow, including the one from the original fibration to #4239, is a 2-neighbor step. A general old fiber and a general new fiber meet in two points. In the first step, write the reduced old coordinates as `(z,x,y)` and the chosen section as `Q = (X/k²,Y/k³)`. Define `c0 = −Y/X mod k²`, `c1 = z*c0 mod k²`, `d0 = (Y+c0*X)/k²`, and `d1 = (z*Y+c1*X)/k²`. The pencil parameter is

```text
s = (z*k*y - c1*x + d1)/(k*y - c0*x + d0).
```

If `b` denotes the published #4239 parameter, then `s = (4098*b − 5568)/(8795*b − 4750)`. Consequently #4239 is the curve cut out by `s = 6886495488/13475653655` on the original surface in reduced coordinates. The first two functions in the [map implementation](neighbors/1152.gp) give the conversion between reduced and repository coordinates. The [pencil data](neighbors/1152-data.gp) give the exact section polynomials and the parameter changes for all three steps.

From the repository root, start GP with `gp -q -s 64M`. The larger stack accommodates the parent 1152 equation. Paste the following to display the original coordinates of a database point and verify the return map.

```gp
read("curves/6782.gp");
P = points()[1];
Q = from_database(P);                    \\ [original T, original x, original y]
Q[1]                                     \\ the original T for this point
to_record_fibration(Q)[1]                \\ 431951/282772
to_database(Q) == P                      \\ 1
```

`to_record_fibration` and `from_record_fibration` transport triples `[base,x,y]` in either direction on the surface. They use affine charts where the displayed denominators are nonzero. Loading each file checks the section equations, the quartic's rational origin, the generic Weierstrass coefficient identities at each step, and all 32 database points on the original surface and back. The construction uses the 2-neighbor method of Elkies and Kumar. [Method and derivation](https://msp.org/ant/2014/8-10/ant-v8-n10-p01-p.pdf)

### Search coverage

The database snapshot contains 78 curves with rank lower bound at least 30. The search considered the fixed 1008 surface and the 1020, 1092, and 1152 families. For a given coefficient point, matching used the rational solutions of `j(T) = j_database`, including the fiber at infinity, and checked the change of Weierstrass coordinates over ℚ.

For 1020 and 1092 the coefficient search covered the listed points and all points of the form `nG + D` with `|n| ≤ 10000` and rational 2-torsion `D`, including the identity. The table gives the Weierstrass model used to enumerate these points and the conversion to the repository coordinates.

| Surface | Weierstrass coefficients | G | (u, v) from an affine point (X, Y) |
| --- | --- | --- | --- |
| 1020 | `[1,1,1,-421,-3157]` | (−13, 24) | ((X − 23)/2, (Y + X − 11)/4) |
| 1092 | `[1,0,1,-137,380]` | (−6, 34) | ((X − 3)/2, (Y + 3X − 7)/2) |

Reduction modulo primes excluded impossible pairs before solving for rational `T`. For 1152 the search covered the four listed points, the four affine points with `u = 0`, and the four points at infinity. These are finite searches in the coefficient parameter.

The equations, rational witness points, and rank lower bounds come from the ICARM Elliptic Curve Rank Leaderboard, maintained under NSF grant DMS 2425401.
