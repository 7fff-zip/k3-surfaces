# High-rank fibers

Each numbered GP script reproduces an ICARM database curve from a K3 surface in this repository.

## Direct specializations on 1008, 1020 and 1092

Here `T` is the repository's base coordinate.

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

## Curves on X27, X440, X552, X588 and X808

| ICARM curve | Exact rank over ℚ | Parent model | T | Construction |
| --- | ---: | --- | ---: | --- |
| [7640](7640.gp) | 6 | [27](../27/QT/27.surface.gp) | 2401/25 | Direct specialization |
| [7704](7704.gp) | 18 | [440](../440/QT/440.surface.gp) | 7043/517 | Direct specialization |
| [7705](7705.gp) | 18 | [440](../440/QT/440.surface.gp) | 7043/517 | 2-isogeny quotient |
| [7702](7702.gp) | 17 | [440](../440/QT/440.surface.gp) | 7115/2756 | Direct specialization |
| [7681](7681.gp) | 13 | [552](../552/QT/552.surface.gp) | −14/15 | Direct specialization |
| [7683](7683.gp) | 13 | [552](../552/QT/552.surface.gp) | 44/273 | Direct specialization |
| [7682](7682.gp) | 13 | [552](../552/QT/552.surface.gp) | −14/15 | 2-isogeny quotient |
| [7693](7693.gp) | 16 | [588](../588/QT/588.surface.gp) | −2/505 | Direct specialization |
| [7699](7699.gp) | 16 | [588](../588/QT/588.surface.gp) | −2/505 | 2-isogeny quotient |
| [7685](7685.gp) | 15 | [588](../588/QT/588.surface.gp) | 179/4 | Direct specialization |
| [7697](7697.gp) | 17 | [808.9-sections-model-01](../808/QT/808.9-sections-model-01.surface.gp) | −254/887 | Direct specialization |
| [7696](7696.gp) | 17 | [808.9-sections-model-01](../808/QT/808.9-sections-model-01.surface.gp) | −254/887 | 2-isogeny quotient |
| [7700](7700.gp) | 16 | [808.9-sections-model-02](../808/QT/808.9-sections-model-02.surface.gp) | −298/515 | Direct specialization |

From the repository root, run `gp -q curves/7640.gp` or substitute any of the other filenames above. From this directory, pass the filename alone. Each script runs offline and includes the exact database coefficients and rational witnesses submitted by Andrii Prokhorenkov. [ICARM data format](https://elliptic-rank.icarm.cloud/api)

In each script, `fiber()` gives the specialized parent equation, `curve()` gives the database model, and `points()` gives the rational witnesses. For a direct specialization, `to_database(Q)` and `from_database(P)` are the two directions of the coordinate change. Each quotient uses the parent fiber's point of order 2 at `(0, 0)`. Its `quotient()` gives the codomain before the final coordinate change, and `to_database(Q)` applies the isogeny and then that change.

The scripts check equality of the five database coefficients and verify every witness on its database model. They use a prime modulus $\ell$, either 5 or 7, and good primes $p$ with exactly one factor of $\ell$ in $\#E(\mathbf F_p)$. Each row consists of discrete logs of the witnesses after reduction and multiplication by $\#E(\mathbf F_p)/\ell$. The rows have full column rank over $\mathbf F_\ell$. A separate good prime with group order prime to $\ell$ excludes rational $\ell$-torsion. PARI/GP's 2-descent confirms the matching rank upper bound.

## Direct specializations on X₈

The coefficient curve for X₈ is the rational line with coordinate $u$. Each numbered file below gives the exact rational values `coefficient_point = [u]`, `t_value` ($T$), and the PARI coordinate change `change`. All three load the same [X₈ surface and section](../8/QCT/8.gp).

| ICARM curve | Rank lower bound over ℚ | Database point on the marked section |
| --- | ---: | --- |
| [6665](6665.gp) | 13 | `points()[1]` |
| [7631](7631.gp) | 13 | `points()[1]` |
| [7632](7632.gp) | 18 | `points()[2]` |

From this directory, run `gp -q 7632.gp` to check the equation, the witness points and the K3 surface certificates. The script constructs `fiber()` with `ellinit(surface(t_value, coefficient_point[1]))` and provides the point maps `to_database(Q)` and `from_database(P)`. At the GP prompt, check the marked point and the database equation with

```gp
Q = sections(t_value, coefficient_point[1])[1];
to_database(Q) == points()[2]                           \\ 1
ellchangecurve(fiber(), change)[1..5] == curve()[1..5]  \\ 1
```

### Recover the X₈ parameters from one point

Write the database equation as $y^2=x^3+a_2x^2+a_4x+a_6$. For the marked witness $P=(x,y)$, let

$$
m=\frac{3x^2+2a_2 x+a_4}{2y},\qquad n=3x+a_2-m^2,\qquad k=2m-\frac{2y}{n}.
$$

Putting $x_{\rm db}=k^2 x_c+x$ and $y_{\rm db}=k^3 y_c+mk^2 x_c+y$ moves the marked point to $(0,0)$ and its tangent to $y_c=0$. The resulting coefficients are $[2m/k,n/k^2,2y/k^3,0,0]$. The compact surface has coefficients $[1+w,-t,-wt,0,0]$ with $t=T/216$, $v=54u$ and $w=vt^2$. Matching coefficients gives the catalogue parameters and scale

$$
u=y\left(\frac{k}{3n}\right)^3,\qquad T=-\frac{216n}{k^2},\qquad h=\frac{6}{k}.
$$

The PARI change from the catalogue fiber to the database model is $[h,h^2a_2/3,0,0]$. For #6665 and #7631, $a_2=0$ and scaling alone recovers the database equation. For #7632, $a_2=1$ and the change includes the translation $h^2/3$. Each script checks its saved parameters against this tangent calculation.

### Verify the fiber rank bounds and the K3 discriminants

For the elliptic curve, each prime in `independence_primes` gives a homomorphism from the rational points into a cyclic group of prime order $\ell$. Multiply the reduction of each point by $\#E(\mathbf F_p)/\ell$ and take its discrete logarithm in that subgroup. The resulting matrix has column rank 13 modulo 3 for #6665 and #7631, and rank 18 modulo 5 for #7632. A separate good prime excludes rational $\ell$-torsion, so these ranks certify independence over ℤ.

For the K3 surface, the compact equation has split fibers $I_7$ at zero and $I_{11}$ at infinity, together with six simple nodal fibers. The marked section $P$ is disjoint from the zero section. In compact coordinates $x(2P)=t$ has order one at zero and $x(4P)/t^4=(vt^3+t^2)/t^4$ has order one at infinity. Thus $P$ meets component $\pm3$ of each reducible fiber and has height $4-12/7-24/11=8/77$. The fiber, zero section, components and marked section give 19 rational divisor classes and absolute determinant $7\cdot11\cdot(8/77)=8$. [Height and determinant formulas, §§11.8–11.10](https://arxiv.org/pdf/0907.0298#page=51)

Each script counts points on its K3 surface at two good primes. For each surface the reductions give distinct Néron–Severi squareclasses.

| Record | Prime $p$ | Points over $\mathbf F_p$ | Points over $\mathbf F_{p^2}$ | Reduction NS squareclass |
|---|---:|---:|---:|---:|
| 6665, 7631 | 5 | 122 | 1112 | −1 |
| 6665 | 11 | 306 | 17016 | −2 |
| 7631 | 17 | 632 | 88728 | −2 |
| 7632 | 17 | 621 | 88805 | −43 |
| 7632 | 29 | 1418 | 725336 | −7 |

After removing the 19 known classes, Frobenius has polynomial $(X-\varepsilon p)(X^2-sX+p^2)$. The two counts determine $\varepsilon$ and $s$ here. The test $p\nmid s$ rules out roots of unity in the normalized quadratic factor. Over $\mathbf F_{p^2}$ the geometric Picard rank is 20 and the Artin–Tate formula gives squareclass $s^2-4p^2$. Distinct squareclasses at two good primes exclude rank 20 in characteristic zero. [Two-prime comparison](https://arxiv.org/abs/math/0506416)

An index-two enlargement is the only possible proper enlargement of the rank-19 lattice above. The 19 rational classes span it over ℚ, so Galois fixes the full geometric Néron–Severi group and every section is rational. A prime-order section can meet a nonidentity component of at most one of the two fibers. Its maximum single height correction is $30/11<4$, so the torsion is trivial. An index-two enlargement would therefore halve the marked section. Reduce the marked database point modulo 5 for #6665 and #7631, or modulo 23 for #7632. Multiplying it by half the group exponent gives a nonzero point, whereas the same multiple kills every double. This proves saturation and gives $|\det\mathrm{NS}|=8$ for each of the three specialized K3 surfaces.

## Rank 32 on X1152 by changing the fibration

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

## Search coverage for rank lower bounds of at least 30

The database snapshot contains 78 curves with rank lower bound at least 30. The search considered the fixed 1008 surface and the 1020, 1092, and 1152 families. For a given coefficient point, matching used the rational solutions of `j(T) = j_database`, including the fiber at infinity, and checked the change of Weierstrass coordinates over ℚ.

For 1020 and 1092 the coefficient search covered the listed points and all points of the form `nG + D` with `|n| ≤ 10000` and rational 2-torsion `D`, including the identity. The table gives the Weierstrass model used to enumerate these points and the conversion to the repository coordinates.

| Surface | Weierstrass coefficients | G | (u, v) from an affine point (X, Y) |
| --- | --- | --- | --- |
| 1020 | `[1,1,1,-421,-3157]` | (−13, 24) | ((X − 23)/2, (Y + X − 11)/4) |
| 1092 | `[1,0,1,-137,380]` | (−6, 34) | ((X − 3)/2, (Y + 3X − 7)/2) |

Reduction modulo primes excluded impossible pairs before solving for rational `T`. For 1152 the search covered the four listed points, the four affine points with `u = 0`, and the four points at infinity. These are finite searches in the coefficient parameter.

The equations, rational witness points, and rank lower bounds come from the ICARM Elliptic Curve Rank Leaderboard, maintained under NSF grant DMS 2425401.
