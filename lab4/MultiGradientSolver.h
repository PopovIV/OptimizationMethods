#ifndef MULTIGRADIENTSOLVER_H_INCLUDED__
#define MULTIGRADIENTSOLVER_H_INCLUDED__

#include <limits>
#include <functional>
#include <algorithm>
#include <vector>
#include "math.h"

// x0 - first solution
// f - our function
// df - gradient of our function
// eps - our error
// iterations - number of iterations
// alpha_0 - step on zero iteration
// lambda - our coefficient of step shrinking
// delta - our criteria of shinking
vec GradMethodSplitStep(vec x0, std::function<double(vec)> f, std::function<vec(vec)> df, double eps, int& iterations, double alpha_0 = 0.5, double lambda = 0.6, double delta = 0.5)
{
  vec xk = x0;
  vec grad_xk = df(xk);
  double alpha_k = alpha_0;
  iterations = 0;
  std::cout << std::endl << std::endl;
  while (grad_xk.lenSquared() >= eps * eps)
  {
    //for printing results 
    //
    std::cout << "k:" << iterations << std::endl;
    std::cout << "xk:";
    xk.print();
    

    //std::cout << "->";
    //
    double newF = f(xk - grad_xk * alpha_k);
    double oldF = f(xk);
    bool was_nan = false;
    double alpha_prev = alpha_k;
    alpha_k = alpha_0;

    if (isnan(newF) || newF == INFINITY)
      was_nan = true;

    int ik = 0;
    while (isnan(newF) || newF == INFINITY || newF - oldF > -delta * alpha_k * grad_xk.lenSquared()) {
      alpha_k = alpha_k * lambda;
      ik++;
      newF = f(xk - grad_xk * alpha_k);
    }

    std::cout << "ik: " << ik << std::endl;
    std::cout << "ak: " << alpha_k << std::endl;

    double qk = (xk - grad_xk * alpha_k).len();
    qk /= (xk).len();
    std::cout << "qk: " << qk << std::endl;

    xk = xk - grad_xk * alpha_k;

    if (was_nan)
      alpha_k = alpha_prev;

    grad_xk = df(xk);
    iterations++;
  }
  //std::cout << std::endl << std::endl;
  return xk;
}

double min(double a, double b) {
  return a > b ? b : a;
}

double FindAlpha(vec xk, vec pk, std::function<double(vec)> f, std::function<vec(vec)> df, double eps, int version = 0)
{
  auto phi = [&](double alpha) -> double {
    return f(xk + pk * alpha);
  };

  double b = 50; // max possible alpha value
  double a = 0; // min possible alpha value

  // we use 1d minimization method of Dihotomy but with some improvments to fight with inf-s
  // Ughhhh
  while (phi((b+a) / 2.0) > 1e+5 || fabs(b - a) >= eps)
  {
    double m = (b + a) / 2.0;
    double delta = (b - a) / 100.0;
    double y1 = phi(m - delta);
    double y2 = phi(m + delta);
    if (y1 > y2)
      a = m;
    else
      b = m;
  }
  double res = phi((b + a) / 2.0);
  return (b + a) / 2.0;
}

// x0 - first solution
// f - our function
// df - gradient of our function
// eps - our error
// n - update moments coefficient
// iterations - number of iterations
vec BFGSMethod(vec x0, std::function<double(vec)> f, std::function<vec(vec)> df, double eps, int n, int& iterations)
{
  vec xk = x0;
  vec wk = df(xk) * -1;
  matr Ak(xk.size(), xk.size());
  
  int k = 1;

  for (int i = 0; i < xk.size(); i++)
    Ak[i][i] = 1;
  iterations = 0;
  std::cout << std::endl << std::endl;
  while (wk.lenSquared() >= eps * eps)
  {
    vec pk = Ak * wk;
    
    std::cout << "k:" << iterations << std::endl;
    std::cout << "xk:";
    xk.print();


    // trying to find optimal alpha_k. Ughhhhhh
    double alpha_k = FindAlpha(xk, pk, f, df, eps);

    std::cout << "ak: " << alpha_k << std::endl;


    vec delta_x_k = pk * alpha_k;
    vec delta_w_k = wk;

    //
    //xk.print();
    //std::cout << "->";
    //

    double qk = (xk + pk * alpha_k).len();
    qk /= (xk).len();

    std::cout << "qk: " << qk << std::endl;

    xk = xk + pk * alpha_k;
    wk = df(xk) * -1;
    delta_w_k = wk - delta_w_k;
    if (k % n == 0) 
    {
      k = 1;
      for (int i = 0; i < xk.size(); i++)
        Ak[i][i] = 1;
    }
    else 
    {
      k++;
      // BFGS formula... 
      // Oh no, senpai it is too huge, ahhhhhh
      // I will take it from here https://habr.com/ru/post/333356/

      matr I(Ak.sizeH(), Ak.sizeW());
      double rho = 1.0 / (delta_x_k * delta_w_k);
      for (int i = 0; i < xk.size(); i++)
        I[i][i] = 1;

      matr m1(Ak.sizeH(), Ak.sizeW());
      matr m2(Ak.sizeH(), Ak.sizeW());
      matr m3(Ak.sizeH(), Ak.sizeW());
      m1.fromVectors(delta_x_k, delta_w_k);
      m2.fromVectors(delta_w_k, delta_x_k);
      m3.fromVectors(delta_x_k, delta_x_k);
      m1 = m1 * rho;
      m2 = m2 * rho;

      Ak = (I - m1) * Ak * (I - m2) + m3;
    }
    iterations++;
  }
  std::cout << std::endl << std::endl;
  return xk;
}


#endif