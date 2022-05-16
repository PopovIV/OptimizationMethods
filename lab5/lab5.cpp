#include "Zoytendeik.h"

vec FindFirstPoint(std::function<double(vec)> target_f,
  std::function<vec(vec)> target_df,
  std::vector<std::function<double(vec)>> ineq_restr,
  std::vector<std::function<vec(vec)>> d_ineq_restr, 
  matr A) {

  // must be a solution of Ax = b
  vec newx0(3);
  newx0[0] = 4;
  newx0[1] = 6;
  newx0[2] = 0;

  // min(nu)
  auto new_f = [&](vec x) -> double { return x[x.size() - 1]; };
  auto new_df = [&](vec x) -> vec {
    vec res(x.size()); // all others are 0
    res[x.size() - 1] = 1;
    return res;
  };

  std::vector<std::function<double(vec)>> new_ineq_restr;
  std::vector<std::function<vec(vec)>> new_d_ineq_restr;


  // phi_i(x) <= nu <=> phi_i(x) - nu <= 0
  for (int i = 0; i < ineq_restr.size(); i++) {
    new_ineq_restr.push_back(
      [&, i](vec x) -> double {
      return ineq_restr[i](x) - x[x.size() - 1];
    }
    );
    new_d_ineq_restr.push_back(
      [&, i](vec x) -> vec {
      vec res = d_ineq_restr[i](x);
      res[x.size() - 1] = -1;
      return res;
    });
  }


  // -10 - nu <= 0 <=> nu >= -10
  new_ineq_restr.push_back(
    [&](vec x) -> double {
    return -10 - x[x.size() - 1];
  });
  // -10 - nu <= 0 <=> nu >= -10
  new_d_ineq_restr.push_back(
    [&](vec x) -> vec {
    vec res(x.size());
    res[x.size() - 1] = -1;
    return res;
  });

  double nu = -INFINITY;
  for (auto& ineq : ineq_restr)
  {
    double val = ineq(newx0);
    nu = nu > val ? nu : val;
  }

  newx0.append(nu);

  matr newA = A;
  matr sideA(A.sizeH(), 1);
  newA.concatinateRight(sideA);

  Zoytendeik tmp(newx0, new_f, new_df, new_ineq_restr, new_d_ineq_restr, newA);

  vec sol = tmp.solve(true);

  nu = sol.popend();
  if (nu < 0)
    printf("first point founded\n");

  return sol;
}

int main()
{
  int counter = 0;
  int iterations = 0;
  auto f = [&](vec x) -> double {
    counter++;
    return x[0] * x[0] + 2 * x[1] * x[1] + 2 * x[2] * x[2] + 8 * x[2] + exp(x[0] * x[0] + x[1] * x[1]);
  };

  auto df = [&](vec x) -> vec {
    vec res(x.size()); // !! IMPORTANT
    res[0] = 2 * x[0] * (exp(x[0] * x[0] + x[1] * x[1]) + 1);
    res[1] = 2 * x[1] * (exp(x[0] * x[0] + x[1] * x[1]) + 2);
    res[2] = 4 * x[2] + 8;
    return res;
  };

  // on border

  // inequality
  std::vector<std::function<double(vec)>> in_eq_restrictions = {
    [&](vec x) -> double {return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;},
    [&](vec x) -> double {return 2 * x[0] * x[0] + x[1] * x[1] - 0.5;},
    [&](vec x) -> double {return (x[0] - 1) * (x[0] - 1) + x[1] * x[1] - 1;},// - 1 for inside
  };

  std::vector<std::function<vec(vec)>> d_in_eq_restrictions = {
    [&](vec x) -> vec {
      vec res(x.size()); //// !! IMPORTANT
      res[0] = 2 * x[0];
      res[1] = 2 * x[1];
      res[2] = 2 * x[2] - 0;// - 0 for inside
      return res; }, 
    [&](vec x) -> vec {
      vec res(x.size()); // !! IMPORTANT
      res[0] = 4 * x[0];
      res[1] = 2 * x[1];
      res[2] = 0;
      return res; },
    [&](vec x) -> vec {
      vec res(x.size()); // !! IMPORTANT
      res[0] = 2 * (x[0] - 1);
      res[1] = 2 * x[1];
      res[2] = 0;
      return res; },
  };

  //equality
  matr A(1, 3); // (0, 0, 0)
  A[0][2] = 1; // (0, 0, 1)
  // that means A * x = 0 <=> x[3] = 0
  /*
  vec x0 = vec(3);
  x0[0] = 0.7;
  x0[1] = 0.3;
  x0[2] = 0;
  */
  // find first point
  std::cout <<"FIND FISRT POINT" << std::endl;
  vec x0 = FindFirstPoint(f, df, in_eq_restrictions, d_in_eq_restrictions, A);
  std::cout << "First point f=" << f(x0) << " x =";
  x0.print();
  std::cout << "FIND SOLUTION" << std::endl;

  Zoytendeik z(x0, f, df, in_eq_restrictions, d_in_eq_restrictions, A);

  vec sol = z.solve();

  std::cout << "Solution point f=" << f(sol) << " x =";
  sol.print();
}