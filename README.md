# K3 surfaces in PARI/GP

Start GP in `1020/QCT/`, `1092/QCT/`, or `1152/QCT/` and load the matching `1020.gp`, `1092.gp`, or `1152.gp`.

For a fixed surface over Q(T), start GP in the directory listed below and run its load command. Use `surface()` for the coefficients and `sections()` for the saved sections.

| Discriminant | Coefficient-curve point | Directory | GP load command | Saved sections |
| --- | --- | --- | --- | ---: |
| 27 | — | `27/QT/` | `read("27.gp")` | 1 |
| 27 | — | `27/QT/` | `read("27.1-section-cm03-27-base-change-1.gp")` | 1 |
| 27 | — | `27/QT/` | `read("27.1-section-schuett-torsion-3.gp")` | 2 |
| 440 | — | `440/QT/` | `read("440.gp")` | 10 |
| 552 | — | `552/QT/` | `read("552.gp")` | 10 |
| 588 | — | `588/QT/` | `read("588.gp")` | 10 |
| 808 | — | `808/QT/` | `read("808.gp")` | 10 |
| 808 | — | `808/QT/` | `read("808.9-sections-model-01.gp")` | 52 |
| 808 | — | `808/QT/` | `read("808.9-sections-model-02.gp")` | 54 |
| 1008 | — | `1008/QT/` | `read("1008.gp")` | 17 |
| 1152 | `[4, 4]` | `1152/QT/` | `read("1152.17-u4-v4.gp")` | 17 |
| 1152 | `[4, -7]` | `1152/QT/` | `read("1152.17-u4-vm7.gp")` | 17 |

The models for discriminants 27, 440, 552, 588, and 808 use the names and coordinates from the `weierstrass-surfaces` catalogue.

| File | Contents |
| --- | --- |
| `*.surface.gp` | `surface(...)` returns `[a1,a2,a3,a4,a6]` for `y² + a1·x·y + a3·y = x³ + a2·x² + a4·x + a6`. |
| `*.sections.gp` | `sections(...)` returns the saved affine points `[x,y]` in the same coordinates. Where provided, `independent_section_count` records the number of independent sections. |
| `*.curve.gp` | `parameter_points` lists saved rational points on the coefficient curve defined by `parameter_relation(u,v) = 0`. |

## Specialize a family over Q(C)(T)

For discriminant 1152 and others, the coefficient curve has coordinates `u,v`. Substituting its rational point `(4,4)` leaves `T` free.

```gp
read("1152.gp");

P = parameter_points[1];            \\ [4, 4]
parameter_relation(P[1], P[2])      \\ 0

E = ellinit(surface('T, P[1], P[2]));
S = sections('T, P[1], P[2]);       \\ all sections on that surface

print([#S, vecmin(vector(#S, i, ellisoncurve(E, S[i])))])
\\ [17, 1]

h(Q) = if(#Q == 1, 0, max(poldegree(numerator(Q[1])), 4 + poldegree(denominator(Q[1]))));
G = matrix(#S, #S, i, j, (h(S[i]) + h(S[j]) - h(ellsub(E, S[i], S[j]))) / 2);
[matrank(G), matdet(G)]
\\ [17, 1152]
```

## Check the section lattice for discriminant 1008

```gp
read("1008.gp");

E = ellinit(surface('T));
S = sections('T);

print([#S, vecmin(vector(#S, i, ellisoncurve(E, S[i])))])
\\ [17, 1]

h(Q) = if(#Q == 1, 0, max(poldegree(numerator(Q[1])), 4 + poldegree(denominator(Q[1]))));
G = matrix(#S, #S, i, j, (h(S[i]) + h(S[j]) - h(ellsub(E, S[i], S[j]))) / 2);
[matrank(G), matdet(G)]
\\ [17, 1008]
```
