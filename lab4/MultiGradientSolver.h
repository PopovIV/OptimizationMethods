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
// alpha_0 - step on zero iteration
// lambda - our coefficient of step shrinking
// delta - our criteria of shinking
vec GradMethodSplitStep(vec x0, std::function<double(vec)> f, std::function<vec(vec)> df, double eps, double alpha_0 = 1, double lambda = 0.5, double delta = 0.9)
{
  vec xk = x0;
  vec grad_xk = df(xk);
  double alpha_k = alpha_0;
  while (grad_xk.lenSquared() >= eps * eps)
  {
    double newF = f(xk - grad_xk * alpha_k);
    double oldF = f(xk);
    bool was_nan = false;
    double alpha_prev = alpha_k;

    if (isnan(newF) || newF == INFINITY)
      was_nan = true;

    while (isnan(newF) || newF == INFINITY || newF - oldF > -delta * alpha_k * grad_xk.lenSquared()) {
      alpha_k = alpha_k * lambda;
      newF = f(xk - grad_xk * alpha_k);
    }
    
    xk = xk - grad_xk * alpha_k;

    if (was_nan)
      alpha_k = alpha_prev;

    grad_xk = df(xk);
  }
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
vec BFGSMethod(vec x0, std::function<double(vec)> f, std::function<vec(vec)> df, double eps, int n)
{
  vec xk = x0;
  vec wk = df(xk) * -1;
  matr Ak(xk.size(), xk.size());
  
  int k = 1;

  for (int i = 0; i < xk.size(); i++)
    Ak[i][i] = 1;
  
  while (wk.lenSquared() >= eps * eps)
  {
    vec pk = Ak * wk;
    
    // trying to find optimal alpha_k. Ughhhhhh
    double alpha_k = FindAlpha(xk, pk, f, df, eps);

    vec delta_x_k = pk * alpha_k;
    vec delta_w_k = wk;

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
      double k = 1.0 / (delta_x_k * delta_w_k);
      for (int i = 0; i < xk.size(); i++)
        I[i][i] = 1;

      matr m1(Ak.sizeH(), Ak.sizeW());
      matr m2(Ak.sizeH(), Ak.sizeW());
      matr m3(Ak.sizeH(), Ak.sizeW());
      m1.fromVectors(delta_x_k, delta_w_k);
      m2.fromVectors(delta_w_k, delta_x_k);
      m3.fromVectors(delta_x_k, delta_x_k);
      m1 = m1 * k;
      m2 = m2 * k;
      m3 = m3 * k;

      Ak = (I - m1) * Ak * (I - m2) + m3;
    }
  }
  return xk;
}


#endif