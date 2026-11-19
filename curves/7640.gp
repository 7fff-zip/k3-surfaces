\\ ICARM #7640, submitted by Andrii Prokhorenkov. https://elliptic-rank.icarm.cloud/curve/7640
\\ ICARM Elliptic Curve Rank Leaderboard, NSF DMS 2425401.

curve_id = 7640;
rank_lower_bound = 6;
surface_stem = "27/QT/27";
t_value = 2401/25;

\\ PARI change from the specialized fiber to the database model.
{change = [
  147/125, 0, 0,
  0
];}

curve() = {
  ellinit([
    0, 0, 0,
    0,
    -110868517847561067
  ])
};

points() = {
  [
    [2437566487/1369, -119159156431070/50653],
    [3818607059623/808201, -7458115233655421090/726572699],
    [90659727118423/90231001, -814679128068728811470/857104278499],
    [695423033121529/148596100, -18328992932242260188333/1811386459000],
    [13131991/25, 23071206914/125],
    [3545584240068814123/94109401, 6676236494759972727702013220/912955299101]
  ]
};

independence_modulus = 5;
torsion_test_prime = 13;
{independence_primes = [
  179, 239, 269, 359, 389, 479
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
