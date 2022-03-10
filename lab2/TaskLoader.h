#ifndef TASKLOADER_H_INCLUDED__
#define TASKLOADER_H_INCLUDED__

#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "TransportTask.h"

class TaskLoader
{
  /* I believe file are written like
  * 2 3              - number of storages and stores
  * 1 2 3            - matrix of cost 
  * 4 5 6            - matrix of cost
  * 10 10            - amount we have in storages
  * 3 12             - amount we need in stores 
  */

private:
  bool IsPhantomStore = false;
  bool IsPhantomStorage = false;

public:
  // very stupid one.
  TransportTask load(std::string filename) {
    std::ifstream file(filename);
    std::string   line;
    int N, M;
    std::vector<std::vector<std::string>> matrixC;//
    std::vector<double> A;
    std::vector<double> B;

    int state = 0; // 0 - number of storages and stores, 1 - matrix, 2 - amount in storages, 3 amount in stores

    bool nextIter = false;
    while (std::getline(file, line)) {
      std::vector<std::string> result;
      std::istringstream iss(line);

      for (std::string s; iss >> s; )
        result.push_back(s);

      if (state == 0) {
          N = std::stoi(result[0]);
          M = std::stoi(result[1]);
      }
      else if (state == 1) {
        matrixC.push_back(result);
        if (matrixC.size() != N)
          continue;
      }
      else if (state == 2) {
        for (int i = 0; i < result.size(); i++)
          A.push_back(std::stod(result[i]));
      }
      else if (state == 3) {
        for (int i = 0; i < result.size(); i++)
          B.push_back(std::stod(result[i]));
      }
      state++;
    }

    // now lets fix everything and remember what we fixed/ 
    // FIRST - all elements must be in matrix and vectors
    matr C(N, M);
    vec a(N), b(M);
    double sumA = 0;
    double sumB = 0;

    // make a
    for (int i = 0; i < a.size(); i++)
    {
      a[i] = A[i];
      sumA += A[i];
    }
    
    // make b
    for (int i = 0; i < b.size(); i++)
    {
      b[i] = B[i];
      sumB += B[i];
    }

    // make C
    for (int i = 0; i < C.sizeH(); i++)
      for (int j = 0; j < C.sizeW(); j++)
        C[i][j] = std::stod(matrixC[i][j]);


    // check inequality
    if (sumA != sumB)
    {
      // Oh no! Open Transport task!!!.
      // Close it ASAP
      if (sumA < sumB)
      {
        // add Phantom store
        IsPhantomStore = true;
        a.append(sumB - sumA);

        // We are Phantom Maticies and we will steal your Vector
        matr phantom(1, b.size());
        C.concatinateDown(phantom);
      }
      else
      {
        // add Phantom store
        // add Phantom storage
        IsPhantomStorage = true;
        b.append(sumA - sumB);

        // We are Phantom Maticies and we will steal your Vector
        matr phantom(a.size(), 1);
        C.concatinateRight(phantom);
      }
    }

    try {
      return TransportTask(C, a, b);
    }
    catch (std::exception& ex) {
      throw ex;
    }

  }

  matr retrieveCorrectAnswer(matr x){
    
    if (IsPhantomStore)
    {
      matr new_x(x.sizeH() - 1, x.sizeW());
      for (int i = 0; i < new_x.sizeH(); i++)
        for (int j = 0; j < new_x.sizeW(); j++)
          new_x[i][j]= x[i][j];
      x = new_x;
    }
    
    if (IsPhantomStorage)
    {
      matr new_x(x.sizeH(), x.sizeW() - 1);
      for (int i = 0; i < new_x.sizeH(); i++)
        for (int j = 0; j < new_x.sizeW(); j++)
          new_x[i][j] = x[i][j];
      x = new_x;
    }

    return x;
  }
};

#endif