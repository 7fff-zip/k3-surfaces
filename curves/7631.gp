\\ ICARM #7631, submitted by NDElkies. https://elliptic-rank.icarm.cloud/curve/7631
\\ From ns8-icarm run gp -q curves/7631.gp. ICARM / NSF DMS 2425401.

curve_id = 7631;
rank_lower_bound = 13;
marked_point_index = 1;
surface_label = "8/QCT";

\\ Inputs to surface(T,u) in 8.surface.gp. The coefficient curve has one coordinate u.
{coefficient_point = [
  -221048628537384033482033575782236545162648987748047074772997857748118740666749979107418188525191122528245054765308932968434414175368310442583445380736
  / 199939377930896413438341282159880863312975059413131893895162454465520651836864192201426519924432882515412854083755655173645514960324745265573740719717409
];}
{t_value =
  254519858654703872900529070365993360856147208203276789716756093060299193015246
  / 5529824136519172019485832909692605475042293639253373931217954370494842574609;
}

\\ PARI change from the X_8 fiber to the database model.
\\ x_fiber = scale^2*x_db + change[2],  y_fiber = scale^3*y_db.
{scale =
  23384551439543606669198997266982
  / 74362787310046225966202495031601674647;
}
change = [scale, 0, 0, 0];

record_fiber = [coefficient_point[1], t_value, change];

read_repo(path) = iferr(read(path), err, read(Str("../", path)), errname(err) == "e_FILE");
read_repo("8/QCT/8.curve.gp");
read_repo("8/QCT/8.surface.gp");
read_repo("8/QCT/8.sections.gp");

fiber(t = t_value, c = coefficient_point) = ellinit(surface(t, c[1]));
to_database(Q) = ellchangepoint(Q, change);
from_database(P) = ellchangepointinv(P, change);
\\ ellchangecurve(fiber(), change)[1..5] gives the five database coefficients.

curve() = ellinit([0, 0, 0, -17687486497420251833361117, 0]);
points() = {
  [
    [4315106393697, 2006118853068097482],
    [-2071013618607, 5267660417172039426],
    [-2045317491279, 5255501926276140498],
    [-373201407887, 2559104459647016974],
    [-45645786169431/25, 639917804032132164378/125],
    [1187114715838011/121, 36952847800135652297058/1331],
    [400316020798771/81, 4206672894249390635378/729],
    [-5988617793879/4, 38470826474508207507/8],
    [26668833149827, 135999560819024706782],
    [-759635171988, 3605231298830104482],
    [4624198413312, 4133976683793971232],
    [-29451654948672/49, 1106885883464978678424/343],
    [2393449554503907, 117094298010931080007218]
  ]
};

nondivision_prime = 5;
independence_modulus = 3;
picard_primes = [5, 17];
independence_primes = [23, 47, 59, 131, 167, 191, 227, 239, 263, 311, 347, 383, 419];

\\ Recover the same fiber from the database model (all three have a1=a3=0).
parameters_from_point(db, P) = {
  my(x = P[1], y = P[2], m = (3*x^2 + 2*db.a2*x + db.a4)/(2*y));
  specialization(x, y, m, db.a2)
};

\\ Map each witness to the cyclic subgroup of order ell in E(F_p).
witness_row(db, P, p, ell) = {
  if(!isprime(p) || p == ell || valuation(db.disc,p) != 0 || denominator(P)%p == 0,
    error("unusable witness prime ", p));
  my(e = ellinit(db[1..5], p), size = ellcard(e));
  if(valuation(size,ell) != 1, error("expected one factor of ", ell, " at ", p));
  my(images = vector(#P, i, ellmul(e, P[i], size/ell)), generator = [0]);
  for(i = 1, #P, if(images[i] != [0], generator = images[i]; break));
  vector(#P, i, if(images[i] == [0], 0, elllog(e,images[i],generator,ell)))
};

check() = {
  my(db = curve(), e = fiber(), P = points(), u = coefficient_point[1]);
  my(marked = P[marked_point_index]);
  if(parameters_from_point(db,marked) != [u,t_value,change],
    error("marked point gives different specialization parameters"));
  if(ellchangecurve(e,change)[1..5] != db[1..5],
    error("specialized equation does not reproduce ICARM coefficients"));
  if(to_database(sections(t_value,u)[1]) != marked,
    error("the saved section does not give the marked ICARM point"));
  for(i = 1, #P,
    my(Q = from_database(P[i]));
    if(!ellisoncurve(db,P[i]) || !ellisoncurve(e,Q) || to_database(Q) != P[i],
      error("witness point transport failed at ", i)));
  my(p = nondivision_prime, ell = independence_modulus);
  if(!isprime(p) || !isprime(ell) || p == ell || valuation(db.disc,p) != 0,
    error("invalid torsion-test prime"));
  my(ep = ellinit(db[1..5],p), exponent = ellgroup(ep)[1]);
  if(ellcard(ep)%ell == 0 || exponent%2 != 0 || ellmul(ep,marked,exponent/2) == [0],
    error("failed the torsion or nondivision certificate at ", p));
  my(rows = vector(#independence_primes, i, witness_row(db,P,independence_primes[i],ell)));
  if(#P != rank_lower_bound || matrank(Mod(Mat(rows~),ell)) != #P,
    error("the witness points failed the independence certificate"));
  my(red = vector(#picard_primes, i, picard_reduction(54*u,picard_primes[i])));
  if(#Set(vector(#red, i, red[i][6])) < 2, error("need distinct NS squareclasses"));
  print("Verified ICARM #", curve_id, ", T = ", t_value, ", points checked = ", #P,
        ", K3 Picard rank 19, |disc NS| = 8");
  1
};

iferr(check(), err, print(err); quit(1));
