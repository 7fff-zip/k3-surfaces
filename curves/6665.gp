\\ ICARM #6665, submitted by NDElkies. https://elliptic-rank.icarm.cloud/curve/6665
\\ From ns8-icarm run gp -q curves/6665.gp. ICARM / NSF DMS 2425401.

curve_id = 6665;
rank_lower_bound = 13;
marked_point_index = 1;
surface_label = "8/QCT";

\\ Inputs to surface(T,u) in 8.surface.gp. The coefficient curve has one coordinate u.
{coefficient_point = [
  1712659481719971044539726809931343046715947074669711502096092761045611566568742743481738296714087418055281057355406920376598860293802090725774304
  / 1825810174935033398589847280889207374903438250807666364663172775137400343865025111680055785926416632499869841647844489853458258084954517921247670049
];}
{t_value =
  -1153696501480194030830909298442344806265123479065141659395811286163570564811
  / 316642531665071164669823691718169745955412744677069849304981372649291601;
}

\\ PARI change from the X_8 fiber to the database model.
\\ x_fiber = scale^2*x_db + change[2],  y_fiber = scale^3*y_db.
{scale =
  39039395160715901832544745455863
  / 1125419977901709659658109328439589102;
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

curve() = ellinit([0, 0, 0, -57993414599174271251532, 0]);
points() = {
  [
    [760466147382, 629032755220242588],
    [-240275722032, 7918967200776984],
    [-56136529968, 55485542698226088],
    [282795191382, 78840060066927588],
    [-120332633091, 72360836522272221],
    [-233185987202, 29044780446500084],
    [-152033460578, 72820407730835212],
    [100098581716772164, 31669549538347975909338236],
    [8355502534564/9, 23325627346112161484/27],
    [-86834602146075/529, 868607321664539336145/12167],
    [241675979908, 10001057371096484],
    [1745660114710374/3721, 62593071464191557648756/226981],
    [-4486587072, 16127684012617584]
  ]
};

nondivision_prime = 5;
independence_modulus = 3;
picard_primes = [5, 11];
independence_primes = [11, 47, 59, 83, 131, 167, 191, 227, 239, 263, 311, 347, 419];

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
