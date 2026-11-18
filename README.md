# K3 surfaces in PARI/GP

Start GP in `1020/QCT/`, `1092/QCT/`, or `1152/QCT/` and load the matching `1020.gp`, `1092.gp`, or `1152.gp`.

For the fixed surface over Q(T), start GP in `1008/QT/` and run `read("1008.gp")`. Use `surface()` for the coefficients and `sections()` for its 17 sections.

| File | Contents |
| --- | --- |
| `*.surface.gp` | `surface(...)` returns `[a1,a2,a3,a4,a6]` for `y² + a1·x·y + a3·y = x³ + a2·x² + a4·x + a6`. |
| `*.sections.gp` | `sections(...)` returns the saved affine points `[x,y]` in the same coordinates. |
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

## Specialize a family over Q(T)

For discriminant 1008 and other, there is no available coefficient curve.

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
