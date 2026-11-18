independent_section_count = 1;
sections(T='T, u='u) = {
  my(t = T/216, w = 54*u*t^2);
  [[3*((1+w)^2 - 4*t), -108*w*t]]
};
