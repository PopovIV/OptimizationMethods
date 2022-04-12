#ifndef GRADIENTSOLVER_H_INCLUDED__
#define GRADIENTSOLVER_H_INCLUDED__

#include <limits>
#include <functional>
#include <algorithm>


#define GOLDEN_RATIO 1.6180339887498948482045868343656381177203091798057628621

// (a, b) - is where we trying to find solution
// f - our function
// eps - our error
double GoldenRatioMethod(double a, double b, std::function<double(double)> f, double eps)
{
  int flag = 2;
  double y1 = 0;
  double y2 = 0;

  while (fabs(b - a) >= eps)
  {
    double x1 = b - (b-a) / GOLDEN_RATIO;
    double x2 = a + (b-a) / GOLDEN_RATIO;
    
    if (flag == 2)
    {
        y1 = f(x1);
        y2 = f(x2);
    }
    if (flag == 0)
    {
        y2 = y1;
        y1 = f(x1);
    }
    if (flag == 1)
    {
        y1 = y2;
        y2 = f(x2);
    }
    if (y1 >= y2)
      a = x1, flag = 1;
    else
      b = x2, flag = 0;
  }
  return (b+a) / 2.0;
}

// (a, b) - is where we trying to find solution
// f - our function
// eps - our error
double DichotomyMethod(double a, double b, std::function<double(double)> f, double eps)
{
  while (fabs(b - a) >= eps)
  {
    double m = (b + a) / 2.0;
    double delta = (b-a) / 100.0;
    double y1 = f(m - delta);
    double y2 = f(m + delta);
    if (y1 >= y2)
      a = m;
    else
      b = m;
  }
  return (b + a) / 2.0;
}


#endif