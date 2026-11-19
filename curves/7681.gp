\\ ICARM #7681, submitted by Andrii Prokhorenkov. https://elliptic-rank.icarm.cloud/curve/7681
\\ ICARM Elliptic Curve Rank Leaderboard, NSF DMS 2425401.

curve_id = 7681;
rank_lower_bound = 13;
surface_stem = "552/QT/552";
t_value = -14/15;

\\ PARI change from the specialized fiber to the database model.
{change = [
  8/225, -6397887953792/50625, -1623952/75,
  6233922881312218624/2278125
];}

curve() = {
  ellinit([
    1, 0, 0,
    -31552033597777759696195,
    2155151814867626298158228643583937
  ])
};

points() = {
  [
    [150953370014, 28844826437572793],
    [82242822494, -10793754609773527],
    [-2926618066, -47407461025808887],
    [206037974774, 66339071653159673],
    [176154505454, 45423093573643913],
    [105289030244, 528917318273033],
    [34472052974, 33293410872539273],
    [-139309387306, -62024557084644367],
    [-5732620306, -48330518056126327],
    [195020947694, 58473051762241673],
    [1036225748549/16, 1253088485758477127/64],
    [2417969520175574, 118898459424811587642353],
    [16881077578125902/101761, 1251424886973648809064407/32461759]
  ]
};

independence_modulus = 7;
torsion_test_prime = 13;
{independence_primes = [
  19, 61, 139, 211, 233, 281, 331, 347, 353, 499, 541, 563,
  571
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
