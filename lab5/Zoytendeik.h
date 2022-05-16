#pragma once
#include "simplex/math.h"
#include <iostream>
#include <functional>
#include "simplex/TaskLoader.h"
#pragma warning(disable: 4996)

#define EPS 1e-4
#define GOLDEN_RATIO 1.6180339887498948482045868343656381177203091798057628621

class Zoytendeik {
  vec x0;
  vec Xi = vec(1);
  std::function<double(vec)> target_f;
  std::function<vec(vec)> target_df;
  std::vector<std::function<double(vec)>> ineq_restr;
  std::vector<std::function<vec(vec)>> d_ineq_restr;
  matr A;
  double lambda;
  double delta;

  // function to find almost Active
  std::vector<int> findAlmostActive()
  {
    std::vector<int> res;
    for (int i = 0; i < ineq_restr.size(); i++) 
    {
      double val = ineq_restr[i](x0);
      if (-delta <= val && val <= 0)
        res.push_back(i);
    }
    return res;
  }

  // function to find almost Active
  std::vector<int> findNonAlmostActive()
  {
    std::vector<int> res;
    for (int i = 0; i < ineq_restr.size(); i++)
    {
      double val = ineq_restr[i](x0);
      if (!(-delta <= val && val <= 0))
        res.push_back(i);
    }
    return res;
  }

  // function for finding direction
  vec makeAndDoSimplex(std::vector<int> almostActive) 
  {
    // lets do stupid thins yahoooooo
    FILE *F = fopen("tmp_simplex.txt", "w");
    
    std::vector<vec> ineq_xk;
    for (auto& i : almostActive)
      ineq_xk.push_back(d_ineq_restr[i](x0));
    
    // write our func
    vec df = target_df(x0);
    for (int el = 0; el < df.size(); el++)
      fprintf(F, "%17g ", df[el]);

    fprintf(F, "%17g <= 0", -1.0);
    fprintf(F, "\n");


    // write all ineq
    for (int line = 0; line < ineq_xk.size(); line++)
    {
      for (int el = 0; el < ineq_xk[line].size(); el++)
        fprintf(F, "%17g ", ineq_xk[line][el]);
      
      fprintf(F, "%17g <= 0", -Xi[line]);
      fprintf(F, "\n");
    }

    // write all eq
    for (int line = 0; line < A.sizeH(); line++)
    {
      for (int el = 0; el < A.sizeW(); el++)
        fprintf(F, "%17g ", A[line][el]);

      fprintf(F, "0 == 0");
      fprintf(F, "\n");
    }

    // Add bounds
    for (int line = 0; line < x0.size(); line++)
    {
      for (int el = 0; el < x0.size(); el++)
        fprintf(F, "%d ", el == line ? 1 : 0);

      fprintf(F, "0 <= 1");
      fprintf(F, "\n");

      for (int el = 0; el < x0.size(); el++)
        fprintf(F, "%d ", el == line ? 1 : 0);

      fprintf(F, "0 >= -1");
      fprintf(F, "\n");
    }

    // no params with sign
    fprintf(F, "\n");
    fprintf(F, "\n");

    // target function
    for (int el = 0; el < x0.size(); el++)
      fprintf(F, "0 ");
    fprintf(F, "1\n");
    
    fclose(F);

    TaskLoader t;
    LinearProgTask task = t.load("tmp_simplex.txt").first;
    vec sol(0);

    // YAHOOO
    sol = task.tableSimplexMethod();
    if (task.checkSolution(sol))
    {
      sol = task.retrieveCorrectAnswer(sol);
      sol = t.retrieveCorrectAnswer(sol);

      return sol;
    }
    else // in case of bad simplex solution
    {
      sol = task.retrieveCorrectAnswer(sol);
      sol = t.retrieveCorrectAnswer(sol);

      if (sol == vec(sol.size()))
      {
        sol = df * -1.0;
        for (int i = 0; i < ineq_xk.size(); i++)
          sol = sol + ineq_xk[i] * -0.1;
      
        double len2 = sol * sol;
        sol = sol * (1.0 / sqrt(len2));

        sol.append(0);
      }
      return sol;
    }
  }

  double getAlpha(double eta, vec s) {
    double alpha = 1;
    while (true)
    {
      vec xk = x0 + (s * alpha);
      bool first = target_f(xk) <= target_f(x0) + 0.5 * eta * alpha;
      bool second = true;

      for (int i = 0; i < ineq_restr.size(); i++)
      {
        double val = ineq_restr[i](xk);
        second = second && val <= 0.0;
      }

      if (first && second)
      {
        return alpha;
      }
      alpha *= lambda;
    }
  }



public:
  // x0 - first point
  // target_f - target function
  // target_df - derivative of target function
  // ineq_restr array of inequality restrictions
  // d_ineq_restr array of derivatives of inequality restrictions
  // A - matrix of equality restrictions (Ax = 0)
  // lambda - division parameter
  // delta - for determine "Almost passive"
  Zoytendeik(
    vec x0,
    std::function<double(vec)> target_f,
    std::function<vec(vec)> target_df,
    std::vector<std::function<double(vec)>> ineq_restr,
    std::vector<std::function<vec(vec)>> d_ineq_restr,
    matr A,
    double lambda = 1/ GOLDEN_RATIO,
    double delta = 1.1) : 
    x0(x0), 
    target_f(target_f), 
    target_df(target_df), 
    ineq_restr(ineq_restr), 
    d_ineq_restr(d_ineq_restr),
    A(A),
    lambda(lambda), 
    delta(delta)
  {
    //these "Speedup" parameters
    Xi = vec(ineq_restr.size());
    for (int i = 0; i < Xi.size(); i++)
      Xi[i] = 1;
  }

  // solve task
  vec solve(bool isFirstX = false)
  {
    int k = 0;
    //if (isFirstX == false)
    //   x0[1] = 1;
    while (true) 
    {
      std::cout << "Iteration " << k << std::endl;
      std::cout << "xk= ";
      x0.print();
      std::cout << "f_xk = " << target_f(x0) << std::endl;
      std::cout << "delta = " << delta << std::endl;


      std::vector<int> almostActive = findAlmostActive();
      std::vector<int> nonAlmostActive = findNonAlmostActive();

      // returns direction concatinated with eta
      vec simplex_solution = makeAndDoSimplex(almostActive);

      // unpack this data
      vec s = vec(simplex_solution.size() - 1);
      double eta = simplex_solution[simplex_solution.size() - 1];
      for (int i = 0; i < s.size(); i++)
        s[i] = simplex_solution[i];


      std::cout << "eta = " << eta << std::endl;
      std::cout << "s = ";
      s.print();

      if (eta < -delta || eta == 0) 
      {
        double alpha = getAlpha(eta, s);
        vec x1 = x0 + s * alpha;
        double prev_f = target_f(x0);
        double cur_f = target_f(x1);
        if (x1 == x0)
        {
          // our step is too small for computer
          std::cout << "WARNING!!! we have so small alpha, that PC cannot proccess"
            "operations with it without \"Arbitrary-precision arithmetic\" due to machine zero problem.\n";
          x0 = x1;
          break;
        }
        x0 = x1;
      }
      else
      {
        delta *= lambda;
      }

      // calc delta_0
      double delta_0 = -delta;
      for (int i = 0; i <ineq_restr.size(); i++)
      {
          double val = ineq_restr[i](x0);
          if(abs(val) > EPS)
              delta_0 = delta_0 < val ? val : delta_0;
      }
      delta_0 *= -1;
      if (isFirstX)
      {
        if (x0[x0.size() - 1] < 0)
          break;
      }
      else
      {
        if (abs(eta) < EPS && delta <= delta_0)
          break;
      }
     
      k++;
    }
    return x0;
  }
};


