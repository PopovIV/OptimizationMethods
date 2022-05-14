#include <iostream>
#include <iomanip>
#include "MultiGradientSolver.h"

int main(void) { 
  int counter = 0;
  int iterations = 0;
  auto f = [&](vec x) -> double {
    counter++;
    return x[0] * x[0] + 2 * x[1] * x[1] + exp(x[0] * x[0] + x[1] * x[1]);
  };

  auto df = [&](vec x) -> vec {
    vec res(2);
    res[0] = 2 * x[0] * (exp(x[0] * x[0] + x[1] * x[1]) + 1);
    res[1] = 2 * x[1] * (exp(x[0] * x[0] + x[1] * x[1]) + 2);
    return res;
  };

  // 0.01
  // числодроблений
  // шаг
  // коэффициент сжатия. Qk
  // (1-eps) / r


  // Xk, Wk, Qk

  double EPSes[] = {0.01};
  vec x_first({ 1, 1 });
  vec q_true({ 0, 0 });
  std::cout << std::fixed << std::setprecision(15);
  std::cout << "*********************Gradient SplitStep Method" << std::endl;
  for (int i = 0; i < sizeof(EPSes) / sizeof(EPSes[0]); i++)
  {
    counter = 0;
    iterations = 0;
    auto q = GradMethodSplitStep(x_first, f, df, EPSes[i], iterations);

    std::cout << "Eps=" << EPSes[i] 
        << " number of iter=" << iterations
        << " found min=" << f(q)
        << " error=" << (q - q_true).len()
        << " counter=" << counter << " at x = ";
    q.print();
  }

  std::cout << "R = " << fabs(df(vec{ 1, 1 })[1]) << std::endl;
  std::cout << "(1-eps)/R=" << (1 - 0.01) / fabs(df(vec{ 1, 1 })[1]) << std::endl;

  
  std::cout << "************BFGS Method" << std::endl;
  vec ansr(0);
  for (int i = 0; i < sizeof(EPSes) / sizeof(EPSes[0]); i++)
  {
    counter = 0;
    iterations = 0;
    auto q = BFGSMethod(x_first, f, df, EPSes[i], 10, iterations);

    std::cout << "Eps=" << EPSes[i]
      << " number of iter=" << iterations
      << " found min=" << f(q)
      << " error=" << (q - q_true).len()
      << " counter=" << counter << " at x = ";
    q.print();
    ansr = q;
  }
  df(ansr).print();

  return 0;

}