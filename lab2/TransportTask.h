#ifndef TRANSPORTTASK_H_INCLUDED__
#define TRANSPORTTASK_H_INCLUDED__

#include "math.h"
#include <limits>
#include <algorithm>

class TransportTask {
private:

  // from lection
  vec a; // size is N 
  vec b; // size is M
  matr C; // size is [N, M]

public:
  // constructor from Common Form
  TransportTask(
    matr c,    // size is [N, M]  and sum(sum(c[i][j] * x[i][j])) -> min
    vec a,  // size is N
    vec b   // size is M
  ) : C(c.sizeH(), c.sizeW()), a(a), b(b)
  {
    C = c;
  }
  #if 0 // I think you may need it
  vec extremePointMethod(void) {
  
    if(M > N)
      throw std::exception("ERROR: incorrect size of task(M > N) in extreme point method");

    std::vector<int> setOfIndexies(A.sizeW());
    for (int i = 0; i < setOfIndexies.size(); i++)
        setOfIndexies[i] = i;
    vec solution(0);
    double minOfFunction = std::numeric_limits<double>::max();
    std::vector<std::vector<int>> vectorOfIndexies = combine(N, M);
    for(auto& indexies : vectorOfIndexies){//checks every possible combination of columns
      //get matrix for linear system
      matr subMatr(0,0);
      try {
        subMatr = A.getSubMatrix(indexies);
      }
      catch (std::exception& ex) {
        throw ex;
      }
      //check det
      try {
        if (subMatr.determinant() == 0)
          continue;
      }
      catch (std::exception& ex) {
        throw ex;
      }
      //solve linear system
      vec systemSolution = rotationMethod(subMatr, b);

      //check that solution is >= 0
      bool flag = true;
      for (int i = 0; i < systemSolution.size(); i++)
        if (systemSolution[i] < 0) {
          flag = false;
          break;
        }//bad practice, ask Vova how to do better
      if (!flag)
        continue;
      //add zeros to solution 
      vec newSolution(c.size());
      for (int i = 0; i < indexies.size(); i++)
        newSolution[indexies[i]] = systemSolution[i];

      //check if that vector is maybe better solution
      double functionVal = 0;
      for (int i = 0; i < newSolution.size(); i++)
        functionVal += newSolution[i] * c[i];
      if (functionVal < minOfFunction) {
        minOfFunction = functionVal;
        solution = newSolution;
      }
    }

    if (solution.size() == 0)
      throw std::exception("No solutions found in extreme point");
    return solution;

  }
  #endif

  //function to get correct answer
  #if 0
  vec retrieveCorrectAnswer(vec answer){
    vec c_true(cOriginalSize);

    //Fill indexies that we didn't change
    int i = 0;
    for (i = 0; i < (cOriginalSize - numOfUnsigned); i++)
      c_true[i] = answer[i];

    //Fill indexies that we did change
    for (int j = answer.size() - numOfUnsigned; j < answer.size(); j++, i++)
      c_true[i] = answer[i] - answer[j];

    //return
    return c_true;
  }
  #endif
};

#endif