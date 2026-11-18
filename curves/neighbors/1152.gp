\\ Birational maps from the repository's X1152 at (4,4) to the record fibrations.
\\ Each step passes through a quartic W^2 = F(t,s), where s is the new base.
read_repo("1152/QCT/1152.curve.gp");
read_repo("1152/QCT/1152.surface.gp");
read_repo("curves/neighbors/1152-data.gp");

\\ The reduced source uses z as base. These give the repository T and scale.
original_t(z) = (-1498959*z - 65454543)/(64490776*z - 36734758);
original_scale(z) = 35997971672139965871075/(9954091672436*(64490776*z - 36734758)^2);
inverse_mobius(f, var, value) = {
  my(n = numerator(f), d = denominator(f));
  (polcoef(n,0,var) - value*polcoef(d,0,var))/(value*polcoef(d,1,var) - polcoef(n,1,var))
};
short_model(a, k = 1) = {
  my(e = ellinit(a));
  [0, 0, 0, -k^4*e.c4/48, -k^6*e.c6/864]
};
short_to_long(a) = [1, (4*a[2]+a[1]^2)/12, a[1]/2, a[3]/2];
published_fiber(t) = if(neighbor_steps == 1, published_4239(t), published_6782(t));
fiber() = ellinit(published_fiber(t_value));

\\ Convert between the quartic with a chosen rational origin and its short Weierstrass model.
neighbor_forward(pt,pencil,base,scale,quartic_data,t,s,h,k,v) = {
  my(u = pt[1], x = pt[2], raw = pencil(u,x,pt[3]), nt = inverse_mobius(base(t),t,raw));
  my(w = subst(subst((2*h^2*x-v)/k,t,u),s,raw), coeff = subst(quartic_data,s,raw));
  my(z = u-coeff[1], q = coeff[2], xx = (2*q*(w+q)+coeff[6]*z)/z^2);
  my(yy = ((xx^2-4*q^2*coeff[3])*z-coeff[6]*xx-2*q^2*coeff[4])/(2*q));
  [nt,9*(xx+coeff[5]/3)/scale(nt)^2,27*yy/scale(nt)^3]
};
neighbor_backward(pt,base,scale,quartic_data,t,s,pencil_coeffs) = {
  my(nt = pt[1], raw = base(nt), coeff = subst(quartic_data,s,raw), q = coeff[2]);
  my(xx = pt[2]*scale(nt)^2/9-coeff[5]/3, yy = pt[3]*scale(nt)^3/27);
  my(z = (2*q*yy+coeff[6]*xx+2*q^2*coeff[4])/(xx^2-4*q^2*coeff[3]));
  my(w = (xx*z^2-coeff[6]*z)/(2*q)-q, u = z+coeff[1]);
  my(values = subst(subst(pencil_coeffs,t,u),s,raw));
  my(x = (values[1]*w+values[5])/(2*values[2]^2));
  [u,x,(values[3]*x-values[4])/(values[1]*values[2])]
};

\\ Intersect the old cubic with the pencil line. The fixed intersection is -Q.
\\ Completing the square in x for the other two intersections gives W^2 = F(t,s).
build_step(index, source, target) = {
  my(t = 'T, s = varlower("pencil_s", 'T));
  my(section = [neighbor_section_1,neighbor_section_2,neighbor_section_3][index](t));
  my(k = section[1], K = k^2, X = section[2], Y = section[3]);
  if(Y^2 != X^3 + source[4]*X*K^2 + source[5]*K^3, error("invalid pencil section ", index));
  my(c0 = lift(-Mod(Y,K)/Mod(X,K)), d0 = (Y+c0*X)/K);
  my(h,c,d,c1,d1);
  if(index == 3,
    h = 1; c = c0+s*K; d = d0+s*X,
    c1 = lift(Mod(t*c0,K)); d1 = (t*Y+c1*X)/K;
    h = t-s; c = c1-s*c0; d = d1-s*d0);
  my(v = (c^2-h^2*X)/K, e = (X*v-2*c*d)/K-h^2*source[4]);
  my(F = (v^2+4*h^2*e)/K);
  if(type(F) != "t_POL" || poldegree(F,t) != 4, error("invalid neighbor quartic ",index));
  my(a4 = polcoef(F,4,t), a3 = polcoef(F,3,t), a2 = polcoef(F,2,t), a1 = polcoef(F,1,t), a0 = polcoef(F,0,t));
  my(ii = 12*a4*a0-3*a3*a1+a2^2, jj = 72*a4*a2*a0+9*a3*a2*a1-27*a4*a1^2-27*a3^2*a0-2*a2^3);
  my(base = [neighbor_base_1,neighbor_base_2,neighbor_base_3][index]);
  my(scale = [neighbor_scale_1,neighbor_scale_2,neighbor_scale_3][index]);
  if(subst(-27*ii,s,base(t)) != scale(t)^4*target[4] || subst(-27*jj,s,base(t)) != scale(t)^6*target[5],
    error("generic Weierstrass identity failed at step ",index));
  my(origin = if(index == 3, neighbor_origin_3(s), [s,subst(v/k,t,s)]));
  if(origin[2]^2 != subst(F,t,origin[1]), error("invalid quartic origin at step ",index));
  my(shifted = subst(F,t,t+origin[1]), b = polcoef(shifted,3,t), cc = polcoef(shifted,2,t), dd = polcoef(shifted,1,t));
  my(quartic_data = [origin[1],origin[2],a4,b,cc,dd]);
  my(pencil = if(index == 3,
    (u,x,y) -> subst((k*y-c0*x+d0)/(K*x-X),t,u),
    (u,x,y) -> subst((t*k*y-c1*x+d1)/(k*y-c0*x+d0),t,u)));
  my(forward = pt -> neighbor_forward(pt,pencil,base,scale,quartic_data,t,s,h,k,v));
  my(backward = pt -> neighbor_backward(pt,base,scale,quartic_data,t,s,[k,h,c,d,v]));
  [forward, backward]
};

build_chain(depth, c) = {
  if(c != [4,4] || parameter_relation(c[1],c[2]) || (depth != 1 && depth != 3), error("unknown X1152 route"));
  my(t = 'T, a = surface(original_t(t),c[1],c[2]), unit = original_scale(t));
  my(source = ellchangecurve(ellinit(a),[unit,0,0,0])[1..5]);
  my(targets = [short_model(published_4239(t)),short_model(published_3900(t),6),short_model(published_6782(t))]);
  my(chain = vector(depth));
  for(i=1,depth, chain[i] = build_step(i,source,targets[i]); source = targets[i]);
  chain
};

\\ Affine maps between triples [base,x,y] on the two presentations of the K3 surface.
from_record_fibration(p) = {
  my(q = concat([p[1]],ellchangepointinv(p[2..3],short_to_long(published_fiber(p[1])))));
  forstep(i=#neighbor_chain,1,-1,q = neighbor_chain[i][2](q));
  my(z = q[1], k = original_scale(z));
  [original_t(z), k^2*q[2], k^3*q[3]]
};
to_record_fibration(p) = {
  my(z = inverse_mobius(original_t('T),'T,p[1]), k = original_scale(z));
  my(q = [z,p[2]/k^2,p[3]/k^3]);
  for(i=1,#neighbor_chain,q = neighbor_chain[i][1](q));
  concat([q[1]],ellchangepoint(q[2..3],short_to_long(published_fiber(q[1]))))
};

\\ The final elliptic fiber has the usual PARI [k,r,s,h] change to the database model.
from_database(p) = from_record_fibration(concat([t_value],ellchangepointinv(p,change)));
to_database(p) = {
  my(q = to_record_fibration(p));
  if(q[1] != t_value,error("point lies on another fiber of the final fibration"));
  ellchangepoint(q[2..3],change)
};

check() = {
  my(db = curve(), p = points(), c = coefficient_point, a = surface('T,c[1],c[2]), q);
  if(ellchangecurve(fiber(),change)[1..5] != db[1..5],error("database model mismatch"));
  for(i=1,#p,
    if(!ellisoncurve(db,p[i]),error("invalid database point ",i));
    q = from_database(p[i]);
    if(!ellisoncurve(ellinit(subst(a,'T,q[1])),q[2..3]) || to_database(q) != p[i],
      error("point transport failed for point ",i)));
  print("Verified ICARM #",curve_id,", C = ",c,", final parameter = ",t_value,", neighbor steps = ",#neighbor_chain,", points checked = ",#p);
  1
};
