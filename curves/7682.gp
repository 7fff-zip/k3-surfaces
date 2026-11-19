\\ ICARM #7682, submitted by Andrii Prokhorenkov. https://elliptic-rank.icarm.cloud/curve/7682
\\ ICARM Elliptic Curve Rank Leaderboard, NSF DMS 2425401.

curve_id = 7682;
rank_lower_bound = 13;
surface_stem = "552/QT/552";
t_value = -14/15;

kernel_point = [0, 0];

\\ PARI change from the isogeny quotient to the database model.
{change = [
  8/225, -6397887953792/50625, -1623952/75,
  6233922881312218624/2278125
];}

curve() = {
  ellinit([
    1, 0, 0,
    -23692879778853666532675,
    3255070248668104198542185468530625
  ])
};

points() = {
  [
    [-157181542150, 55640177475772325],
    [89072606750, 43027614562128425],
    [-170422338400, -48406096634190175],
    [3859451909590/9, 7242517442928225755/27],
    [3466547449853198/289, -204084694238746297400999/4913],
    [82302977456750/5329, -20923306392397266819775/389017],
    [5890283779150/9, 13979789828984472275/27],
    [-20086067912050/121, -68047180782680309525/1331],
    [-199203761650, 8363695386734825],
    [57369525332, 45657819991574441],
    [-854833290475/16, 4230012734309026925/64],
    [-12304125025552/17689, 134564877745234292074067/2352637],
    [-43392580078516450/657721, 35905519832332479793232075/533411731]
  ]
};

independence_modulus = 5;
torsion_test_prime = 13;
{independence_primes = [
  31, 47, 61, 103, 223, 227, 307, 313, 337, 349, 397, 443,
  509
];}

iferr(read("curves/fixed-fiber.gp"), err, read("fixed-fiber.gp"), errname(err) == "e_FILE");

to_database(Q) = ellchangepoint(ellisogenyapply(fiber_isogeny()[2], Q), change);

check() = {
  load_parent();
  check_isogenous_fiber();
  check_rank();
};

iferr(check(), err, print(err); quit(1));
