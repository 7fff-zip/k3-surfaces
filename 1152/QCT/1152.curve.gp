parameter_points = [[4, 4], [4, -7], [1/4, -1], [1/4, 7/4]];
parameter_relation(u = 'u, v = 'v) =
{
  -v^4 + (-2*u + 2)*v^3 + (u^2 + 2*u + 1)*v^2 + (2*u^3 - 2*u^2 + 2*u - 2)*v
    + (-2*u^3 - 2*u^2 - 2*u)
}
