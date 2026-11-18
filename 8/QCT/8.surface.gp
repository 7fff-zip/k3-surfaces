\\ Compact coordinates t = T/216 and v = 54*u.
compact_surface(t, v) = {
  my(w = v*t^2);
  [1+w, -t, -w*t, 0, 0]
};

\\ Return the short equation in the original catalogue coordinates T,u.
surface(T='T, u='u) = {
  my(a = compact_surface(T/216, 54*u));
  my(b2 = a[1]^2 + 4*a[2], b4 = a[1]*a[3], b6 = a[3]^2);
  [0, 0, 0, -27*(b2^2-24*b4), 54*(b2^3-36*b2*b4+216*b6)]
};

\\ Count the resolved K3 over F_p or F_(p^2) in compact coordinates.
surface_point_count(v, p, degree) = {
  my(q = p^degree, g = ffgen(ffinit(p, degree)));
  my(elements = vector(q, i, ((i-1)%p) + ((i-1)\p)*g));
  my(total = 18*q, a, linear, rhs, d);
  \\ The two split multiplicative fibers contribute 7q + 11q.
  for(i = 2, q,
    a = compact_surface(elements[i], v);
    total++;                              \\ point at infinity on this cubic
    for(j = 1, q,
      my(x = elements[j]);
      linear = a[1]*x + a[3];
      rhs = x^3 + a[2]*x^2;
      d = linear^2 + 4*rhs;
      total += if(d == 0, 1, if(issquare(d), 2, 0))));
  total
};

\\ Return [p, N_p, N_(p^2), epsilon, quadratic trace, NS squareclass].
picard_reduction(v, p) = {
  if(!isprime(p) || p < 5 || valuation(v, p) != 0, error("bad prime or parameter at ", p));
  my(z = lift(Mod(v,p)));
  my(generic = ellinit(compact_surface('t, z)));
  my(residual = Mod(1,p) * (generic.disc / 't^7));
  if(poldegree(residual) != 6 || polcoef(residual,0) == 0 || poldisc(residual) == 0,
    error("the six nodal fibers do not stay distinct at ", p));
  my(n1 = surface_point_count(z, p, 1), n2 = surface_point_count(z, p, 2));
  my(a = n1 - 1 - p^2 - 19*p, b = n2 - 1 - p^4 - 19*p^2);
  forstep(eps = -1, 1, 2,
    my(s = a - eps*p);
    if(b == s^2 - p^2 && s%p != 0 && s^2 < 4*p^2,
      return([p, n1, n2, eps, s, core(s^2 - 4*p^2)])));
  error("the counts do not certify reduction Picard rank 20 at ", p)
};
