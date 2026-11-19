\\ ICARM #7683, submitted by Andrii Prokhorenkov. https://elliptic-rank.icarm.cloud/curve/7683
\\ ICARM Elliptic Curve Rank Leaderboard, NSF DMS 2425401.

curve_id = 7683;
rank_lower_bound = 13;
surface_stem = "552/QT/552";
t_value = 44/273;

\\ PARI change from the specialized fiber to the database model.
{change = [
  80/74529, 6179184602907200/5554571841, -19055774/24843,
  -50463954722746022115200/59139526391127
];}

curve() = {
  ellinit([
    1, 0, 0,
    -753476147642486726339363,
    172543550574399322515349393679561217
  ])
};

points() = {
  [
    [246984631462, -38896553291689031],
    [3342163728142, 5914934243858218729],
    [778927238812, -241324120863087431],
    [965905783642, 588150473668848229],
    [591366234374038/289, 13181744866725162587777/4913],
    [720648211852, 61719320760777049],
    [962841296062, -582821627614719431],
    [1492513780942, -1540354541348540471],
    [2223868974742, 3081437292318621529],
    [237745367758, -82741492245641783],
    [260474725051822/361, 461429567258451634291/6859],
    [3123708262012, -5319649460205850631],
    [931222894702, -527659001366819351]
  ]
};

independence_modulus = 5;
torsion_test_prime = 23;
{independence_primes = [
  13, 67, 163, 197, 227, 239, 241, 251, 313, 337, 389, 433,
  443
];}

iferr(read("curves/fixed-fiber.gp"), err, read("fixed-fiber.gp"), errname(err) == "e_FILE");

to_database(Q) = ellchangepoint(Q, change);
from_database(P) = ellchangepointinv(P, change);

check() = {
  load_parent();
  check_direct_fiber();
  check_rank();
};

iferr(check(), err, print(err); quit(1));
