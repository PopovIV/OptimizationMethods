#include <iostream>
#include "GradientSolver.h"
#include <iomanip>
#include <cmath>
#include <limits>
int main(void) { 
  int counter = 0;

  // plot can be found here https://www.desmos.com/calculator/v3jqvsfynk
  auto f = [&](double x) -> double {
    counter++;
    return pow(x, 6) + 3 * x*x + 6 * x - 1;
  };
  
  double EPSes[] = {0.1, 0.01, 0.001, 0.0001};
  double q_true = -0.754877666246693;
  std::cout << "Golden Ratio Method" << std::endl;
  for (int i = 0; i < sizeof(EPSes) / sizeof(EPSes[0]); i++)
  {
    counter = 0;
    double q = GoldenRatioMethod(-1, 0, f, EPSes[i]);

    std::cout << "Eps=" << EPSes[i] 
        << " found min=" <<  f(q) 
        << " at x=" << q << " error=" << std::fixed << std::setprecision(6) << fabs(q - q_true)
        << " func error=" << fabs(f(q) - f(q_true))
        << " counter=" << counter << std::endl;
  }

  std::cout << "Dichotomy Method" << std::endl;
  for (int i = 0; i < sizeof(EPSes) / sizeof(EPSes[0]); i++)
  {
    counter = 0;
    double q = DichotomyMethod(-1, 0, f, EPSes[i]);

    std::cout << "Eps=" << EPSes[i]
      << " found min=" << f(q)
      << " at x=" << q << " error=" << fabs(q - q_true)
      << " func error=" << fabs(f(q) - f(q_true))
      << " counter=" << counter << std::endl;
  }

  return 0;

}