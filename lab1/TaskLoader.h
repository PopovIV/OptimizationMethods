#ifndef TASKLOADER_H_INCLUDED__
#define TASKLOADER_H_INCLUDED__

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "LinearProgTask.h"

class TaskLoader
{
  /* I believe file are written like
  * 1 2 0 0 0 0 == 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 1 2 0 0 0 0 == 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 1 2 0 0 0 0 == 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 1 2 0 0 0 0 >= 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 1 2 0 0 0 0 >= 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 1 2 0 0 0 0 <= 0  - coeffs from line in matr A, symbol of sign and element from vec b
  * 0 3               - index of params which are >= 0
  * 5                 - index of params which are <= 0
  * 1 1 1 -1 1 -1     - vector c
  */

private:
  std::vector<std::pair<int, int>> swapped;
  std::vector<int> negated;

public:
  // very stupid one.
  std::pair<LinearProgTask, LinearProgTask> load(std::string filename) {
    std::ifstream file(filename);
    std::string   line;
    std::vector<std::vector<std::string>> linearSystem;//
    std::vector<int> greaterIndices;
    std::vector<int> lessIndices;
    std::vector<double> c;

    int state = 0; // 0 - system, 1 - greater than 0, 2 - less than 0, 3 - vector c

    int inequalities = 0;
    int equalities = 0;
    bool nextIter = false;
    while (std::getline(file, line)) {
      std::vector<std::string> result;
      std::istringstream iss(line);

      for (std::string s; iss >> s; )
        result.push_back(s);

      if (state == 0) {
        for (int i = 0; i < result.size(); i++)
          if (result[i] == ">=" || result[i] == "<=" || result[i] == "==") {
            if (result[i] == "==")
              equalities++;
            else
              inequalities++;
            linearSystem.push_back(result);
            nextIter = true;
            break;
          }

         if (nextIter) {
           nextIter = false;
           continue;
         }

         state++;
      }

      if (state == 1) {
        for (int i = 0; i < result.size(); i++)
          greaterIndices.push_back(std::stoi(result[i]));
      }
      else if (state == 2) {
        for (int i = 0; i < result.size(); i++)
          lessIndices.push_back(std::stoi(result[i]));
      }
      else if (state == 3) {
        for (int i = 0; i < result.size(); i++)
          c.push_back(std::stod(result[i]));
      }
      state++;
    }

    // now lets fix everything and remember what we fixed/ 
    // FIRST - all elements must be in matrix and vectors
    vec c_to_build(c.size());
    matr A1(inequalities, c.size()), A2(equalities, c.size());
    vec b1(inequalities), b2(equalities);

    // make c
    for (int i = 0; i < c.size(); i++)
      c_to_build[i] = c[i];

    int A1Index = 0, A2Index = 0;

    for (auto& line : linearSystem) {
      if (line[line.size() - 2] == "==") {
        for (int i = 0; i < line.size() - 2; i++) {
          A2[A2Index][i] = std::stod(line[i]);
        }
        b2[A2Index++] = std::stod(line[line.size() - 1]);
      }
      if (line[line.size() - 2] == ">=") {
        // we must change sign, so multiply everything by -1
        for (int i = 0; i < line.size() - 2; i++) {
          A1[A1Index][i] = -1 * std::stod(line[i]);
        }
        b1[A1Index++] = -1 * std::stod(line[line.size() - 1]);
      }
      if (line[line.size() - 2] == "<=") {
        for (int i = 0; i < line.size() - 2; i++) {
          A1[A1Index][i] = std::stod(line[i]);
        }
        b1[A1Index++] = std::stod(line[line.size() - 1]);
      }
  }

    // SECOND - change all elements which are <= 0 to be >= 0. 
    for (int i = 0; i < lessIndices.size(); i++) {
    negated.push_back(lessIndices[i]);

    c[lessIndices[i]] *= -1;

    for (int j = 0; j < A1.sizeH(); j++)
      A1[j][lessIndices[i]] *= -1;

    for (int j = 0; j < A2.sizeH(); j++)
      A2[j][lessIndices[i]] *= -1;
    }

    // THIRD - make every element with sign restrictions a first element
    std::vector<int> allSignedIndices;
    allSignedIndices.insert(allSignedIndices.end(), greaterIndices.begin(), greaterIndices.end());
    allSignedIndices.insert(allSignedIndices.end(), lessIndices.begin(), lessIndices.end());
    std::sort(allSignedIndices.begin(), allSignedIndices.end());
    for (int i = 0; i < allSignedIndices.size(); i++) {
      if (allSignedIndices[i] != i) {
        swapped.push_back(std::make_pair(allSignedIndices[i], i));
        // when we swap elements, we must swap them in vector 'c' and matrices A1 and A2

        std::swap(c[allSignedIndices[i]], c[i]);

        for (int j = 0; j < A1.sizeH(); j++)
          std::swap(A1[j][allSignedIndices[i]], A1[j][i]);

        for (int j = 0; j < A2.sizeH(); j++)
          std::swap(A2[j][allSignedIndices[i]], A2[j][i]);

        allSignedIndices[i] = i;
      }
    }

    //FORTH - make pair of LinearProgTask and dual LinearProgTask
    
    //Create dual task from A1, b1, A2, b2,...
    matr A_dual(A1);
    vec c_dual(b1);
    int N1 = allSignedIndices.size();
    int numOfSigned = A_dual.sizeH();

    //In our LPP we want to find min, so in dual we want to find max
    //all <= to >=
    for (int i = 0; i < A1.sizeH(); i++) {
      for (int j = 0; j < A1.sizeW(); j++)
        A_dual[i][j] *= -1;
      c_dual[i] *= -1;
    }

    //unite and transpose
    try {
      A_dual.concatinateDown(A2);
    }
    catch (std::exception& ex) {
      throw ex;
    }
    try {
      A_dual.transpose();
    }
    catch(std::exception& ex){
      throw ex;
    }
    try {
      c_dual.concatinate(b2);
    }
    catch (std::exception& ex) {
      throw ex;
    }
    //prepare
    matr A1_dual(N1, A_dual.sizeW()), A2_dual(A_dual.sizeH() - N1, A_dual.sizeW());
    vec b1_dual(N1), b2_dual(c.size() - N1);

    for (int i = 0; i < N1; i++) {
      b1_dual[i] = c[i];
      for (int j = 0; j < A_dual.sizeW(); j++)
        A1_dual[i][j] = A_dual[i][j];
    }

    for (int i = N1; i < A_dual.sizeH(); i++) {
      b2_dual[i - N1] = c[i];
      for (int j = 0; j < A_dual.sizeW(); j++)
        A2_dual[i - N1][j] = A_dual[i][j];
    }

    for (int i = 0; i < c_dual.size(); i++)
      c_dual[i] *= -1;

    try {
      return std::make_pair(
          LinearProgTask(c_to_build, A1, b1, A2, b2, N1),// common
          LinearProgTask(c_dual, A1_dual, b1_dual, A2_dual, b2_dual, numOfSigned)//dual
      );
    }
    catch (std::exception& ex) {
      throw ex;
    }

  }

  vec retrieveCorrectAnswer(vec x){
    vec x_true = x;

    //we must do every thing but backwards

    // THIRD - swap
    for (int i = swapped.size() - 1; i >= 0; i--)
    std::swap(x_true[swapped[i].first], x_true[swapped[i].second]);

    // SECOND - negate
    for (int i = 0; i < negated.size(); i++)
    x_true[negated[i]] *= -1;

    //return
    return x_true;
  }
};

#endif