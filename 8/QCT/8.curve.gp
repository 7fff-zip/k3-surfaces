\\ C is the rational line with coordinate u. Each row gives [u, T, change].
\\ The exact inputs are the marked x,y coordinates, tangent slope, and a2.
specialization(x, y, slope, a2) = {
  my(n = 3*x + a2 - slope^2, k = 2*slope - 2*y/n, h = 6/k);
  [y*(k/(3*n))^3, -216*n/k^2, [h, h^2*a2/3, 0, 0]]
};

{record_fibers = [
  specialization(760466147382, 629032755220242588, 4961560707220/3722253, 0),
  specialization(4315106393697, 2006118853068097482, 4423175188115/464906, 0),
  specialization(6873059167553728443935438715689/5447620864,
                 17145815613579789528084408428786749282444601387/402078000730112,
                 25129778539729853952149109723087078805267733590204408309/465606724910358946453999713134336596400739872, 1)
];}
parameter_points = vector(#record_fibers, i, [record_fibers[i][1]]);
parameter_relation(u = 'u) = 0;
