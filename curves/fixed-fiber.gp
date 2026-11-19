read_repo(path) = iferr(read(path), err, read(Str("../", path)), errname(err) == "e_FILE");

fiber() = ellinit(surface(t_value));

load_parent() = {
  read_repo(Str(surface_stem, ".surface.gp"));
  read_repo(Str(surface_stem, ".sections.gp"));
};

check_direct_fiber() = {
  my(e = fiber(), db = curve(), witnesses = points());
  if(ellchangecurve(e, change)[1..5] != db[1..5],
    error("specialized equation does not match ICARM #", curve_id));
  for(i = 1, #witnesses,
    my(Q = from_database(witnesses[i]));
    if(!ellisoncurve(e, Q) || to_database(Q) != witnesses[i],
      error("point transport failed for witness ", i)));
};

fiber_isogeny() = ellisogeny(fiber(), kernel_point);
quotient() = ellinit(fiber_isogeny()[1]);

check_isogenous_fiber() = {
  my(e = fiber(), db = curve(), saved = sections(t_value));
  if(#kernel_point != 2 || !ellisoncurve(e, kernel_point) || ellmul(e, kernel_point, 2) != [0],
    error("expected a nonzero point of order two in the kernel"));
  if(poldegree(fiber_isogeny()[2][3]) != 1 || ellchangecurve(quotient(), change)[1..5] != db[1..5],
    error("the 2-isogeny quotient does not match the database equation"));
  if(to_database(kernel_point) != [0], error("the kernel point did not map to infinity"));
  for(i = 1, #saved,
    if(!ellisoncurve(e, saved[i]) || !ellisoncurve(db, to_database(saved[i])),
      error("isogeny failed on saved section ", i)));
};

\\ Project the reductions onto the order-ell subgroup and take discrete logs.
reduction_row(db, witnesses, p, ell) = {
  if(!isprime(p) || p == ell || db.disc % p == 0 || denominator(witnesses) % p == 0,
    error("invalid independence prime ", p));
  my(e = ellinit(db[1..5], p), size = ellcard(e));
  if(valuation(size, ell) != 1, error("expected one factor of ", ell, " at ", p));
  my(images = vector(#witnesses, i, ellmul(e, witnesses[i], size/ell)), generator = [0]);
  for(i = 1, #images, if(images[i] != [0], generator = images[i]; break));
  vector(#images, i, if(images[i] == [0], 0, elllog(e, images[i], generator, ell)))
};

check_rank() = {
  my(db = curve(), witnesses = points(), ell = independence_modulus, p = torsion_test_prime);
  if(#witnesses != rank_lower_bound || !vecmin(vector(#witnesses, i, ellisoncurve(db, witnesses[i]))),
    error("invalid database witnesses"));
  if(!isprime(ell) || !isprime(p) || p == ell || db.disc % p == 0 || ellcard(ellinit(db[1..5], p)) % ell == 0,
    error("invalid torsion-test prime ", p));
  my(rows = vector(#independence_primes, i, reduction_row(db, witnesses, independence_primes[i], ell)));
  if(matrank(Mod(Mat(rows~), ell)) != #witnesses,
    error("witness independence certificate failed"));

  my(seed = getrand(), bounds);
  setrand(1);
  bounds = ellrank(db, 0, witnesses)[1..2];
  setrand(seed);
  if(bounds != [rank_lower_bound, rank_lower_bound], error("unexpected descent bounds ", bounds));
  print("Verified ICARM #", curve_id, ", T = ", t_value, ", rank = ", bounds[1],
        ", independent witnesses = ", #witnesses);
};
