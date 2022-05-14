#include "Zoytendeik.h"
int main()
{
  int counter = 0;
  int iterations = 0;
  auto f = [&](vec x) -> double {
    counter++;
    return x[0] * x[0] + 2 * x[1] * x[1] + x[2] * x[2] + exp(x[0] * x[0] + x[1] * x[1]);
  };

  auto df = [&](vec x) -> vec {
    vec res(3);
    res[0] = 2 * x[0] * (exp(x[0] * x[0] + x[1] * x[1]) + 1);
    res[1] = 2 * x[1] * (exp(x[0] * x[0] + x[1] * x[1]) + 2);
    res[2] = 2 * x[2];
    return res;
  };

  // on border

  // inequality
  std::vector<std::function<double(vec)>> in_eq_restrictions = {
    [&](vec x) -> double {return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 4;},
    [&](vec x) -> double {return x[0] * x[0] + 2 * x[1] * x[1] - 2; },
    [&](vec x) -> double {return x[0] * x[0] + x[2] * x[2] - 2;},
  };

  std::vector<std::function<vec(vec)>> d_in_eq_restrictions = {
    [&](vec x) -> vec {
      vec res(3);
      res[0] = 2 * x[0];
      res[1] = 2 * x[1];
      res[2] = 2 * x[2];
      return res; },
    [&](vec x) -> vec {
      vec res(3);
      res[0] = 2 * x[0];
      res[1] = 4 * x[1];
      res[2] = 0;
      return res; },
    [&](vec x) -> vec {
      vec res(3);
      res[0] = 2 * x[0];
      res[1] = 0;
      res[2] = 2 * x[2];
      return res; },
  };

  //equality
  matr A(1, 3); // (0, 0, 0)
  A[0][2] = 1; // (0, 0, 1)
  // that means A * x = 0 <=> x[3] = 0

  vec x0 = vec(3);
  x0[0] = 0.7;
  x0[1] = 0.3;
  x0[2] = 0;

  Zoytendeik z(x0, f, df, in_eq_restrictions, d_in_eq_restrictions, A);

  vec sol = z.solve();

  sol.print();
}