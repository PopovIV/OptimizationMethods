#ifndef LINEARPROGTASK_H_INCLUDED__
#define LINEARPROGTASK_H_INCLUDED__

#include "math.h"
#include <limits>
#include <algorithm>
#include "OtherSimplex.h"

class LinearProgTask {
private:

  int N; // number of parameters
  int M; // number of equations

  // !!!! data stored in canonical form !!!

  // that means that we trying to find 
  // vec x; with size N, and all elements are >= 0

  // from lection
  vec c; // size is N 
  matr A; // size is [N, M]
  vec b; // size is M

  //additional data to get answer for common form
  int cOriginalSize;
  int numOfUnsigned;
  std::vector<int> startingBasis;

  void genCombinations(int n, int k, int i, int nc, std::vector<int>& temp, std::vector<std::vector<int>>& ret) {
    if (i == n) {
      if (nc == k)
        ret.emplace_back(temp);
      return;
    }
    if (nc < k)
      temp[nc] = i;
    genCombinations(n, k, i + 1, nc + 1, temp, ret);
    genCombinations(n, k, i + 1, nc, temp, ret);
  }
  std::vector<std::vector<int>> combine(int n, int k) {
    std::vector<int> temp(k, 0);
    std::vector<std::vector<int>> ret;
    genCombinations(n, k, 0, 0, temp, ret);
    return ret;
  }
  //function to solve linear system of eq...
  vec rotationMethod(matr A, vec B) {

    double C, S;
    vec X(B.size());
    int j = 0;

    for (int i = 0; i < A.sizeW()-1; i++) {
      for (int m = i + 1; m < A.sizeW(); m++) {
        C = A[i][j] / sqrt(A[i][j] * A[i][j] + A[m][j] * A[m][j]);
        S = A[m][j] / sqrt(A[i][j] * A[i][j] + A[m][j] * A[m][j]);
        for (int k = j; k < A.sizeW(); k++) {
          double A1 = C * A[i][k] + S * A[m][k];
          double A2 = C * A[m][k] - S * A[i][k];
          A[i][k] = A1;
          A[m][k] = A2;
        }
        double B1 = C * B[i] + S * B[m];
        double B2 = C * B[m] - S * B[i];
        B[i] = B1;
        B[m] = B2;
      }
      j++;
    }

    //get X
    for (int i = A.sizeW() - 1; i >= 0; i--) {
      double alpha = B[i];
      for (int j = A.sizeW() - 1; j > i; j--)
        alpha -= A[i][j] * X[j];
      X[i] = alpha / A[i][i];
    }

    return X;

  }

public:
  vec getC()
  {
    return c;
  }
  // constructor from Common Form
  LinearProgTask(
    vec c,    // size is N and c[N] * x[N] -> min
    matr A1,  // size is [N, M1] and A1 * x <= b1
    vec b1,   // size is M1
    matr A2,  // size is [N, M2] and A2 * x == b2
    vec b2,   // size is M2
    int N1    // x[N1] >= 0
  ): c(c), A(A1), b(b1), cOriginalSize(c.size())
  {
    numOfUnsigned = A1.sizeW() - N1 > 0 ? A.sizeW() - N1 : 0;
    //Unite A1 and A2 to one matrix
    A.concatinateDown(A2);
    //Unite b1 and b2 to one vector
    b.concatinate(b2);
    //FIRST STEP 
    //Change <= to =
    matr additionalMatrix(A.sizeH(), A1.sizeH());
    for (int i = 0, j = 0; i < additionalMatrix.sizeW(); i++) {
      additionalMatrix[i][j++] = 1;
      this->c.append(0);//!!

      this->startingBasis.push_back(this->c.size() - 1);
    }
    A.concatinateRight(additionalMatrix);

    //SECOND STEP
    //Check b for negative elements, if element is negetive => a[i] *= -1, b[i] *= -1
    for (int i = 0; i < b.size(); i++) {
      if (b[i] < 0) {
        b[i] *= -1;
        for(int j = 0; j < A.sizeW(); j++)
          A[i][j] *= -1;
      }
    }

    //THIRD STEP 
    //Check remaining not >= sign variables 
    for (int i = N1; i < A1.sizeW(); i++) {
      matr additionalVector(A.sizeH(), 1);
      for (int j = 0; j < additionalVector.sizeH(); j++)
        additionalVector[j][0] = -1 * A[j][i];
      A.concatinateRight(additionalVector);
      this->c.append(-1 * c[i]);//!!
    }

    N = A.sizeW();
    M = A.sizeH();

  }

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

  std::vector<std::pair<int, int>> fixSimplexTable(std::vector<std::vector<double>> &simplexTable)
  {
    std::vector<std::pair<int, int>> swappers;
    // forward
    for (int line = 0; line < M; line++)
    {
      // if something goes wrong
      if (simplexTable[line][line] == 0)
      {
        // find in vertical line
        int j = line + 1;
        for (; j < M; j++)
          if (simplexTable[j][line] != 0)
          {
            auto tmp = simplexTable[j];
            simplexTable[j] = simplexTable[line];
            simplexTable[line] = tmp;
            break;
          }
        if (j == M)
        {
          printf("TODO");
          // try to find horisontaly
        }
      }

      // divide line
      for (int j = 0; j < N + 1; j++)
        simplexTable[line][j] /= simplexTable[line][line];

      // substract and divide other lines
      for (int l2 = line + 1; l2 < M; l2 ++)
        for (int j = 0; j < N + 1; j++)
          simplexTable[l2][j] -= simplexTable[line][j] * simplexTable[l2][line];
    }

    // backward

    
    return swappers;
  }


  bool isIn(int A, std::vector<int> B)
  {
    for (auto &i : B)
      if (A == i)
        return true;
    return false;
  }

  bool checkSolution(vec answer)
  {
    vec check = A * answer;
    
    for (int i = 0; i < b.size(); i++)
      if (check[i] - b[i] > 1e-6)
        return false;
    return true;
  }

  std::vector<int> doSimplex(std::vector<std::vector<double>> &simplexTable, std::vector<int> basis)
  {
    int b_col = simplexTable[0].size() - 1;
    int c_row = simplexTable.size() - 1;

    while (true)
    {
      int main_column_index = 0;
      int main_row_index = -1;

      // step 1 - check if optimal
      bool is_optimal = true;
      for (int i = 0; i < b_col; i++)
        if (simplexTable[c_row][i] < 0)
        {
          is_optimal = false;
          main_column_index = i;
          break;
        }
      if (is_optimal)
        break;

      // step 2 determine main column
      // Blend rule
      for (int i = 0; i < b_col; i++)
      {
        if ((simplexTable[c_row][i] < 0) && !isIn(i, basis))
        {
          main_column_index = i;
          break;
        }
      }
      
      // step 3 determine main row
      // blends rule
      std::vector<double> delta;
      delta.resize(c_row);
      for (int i = 0; i < delta.size(); i++)
        if (simplexTable[i][main_column_index] > 0)
        {
          delta[i] = simplexTable[i][b_col] / simplexTable[i][main_column_index];
        }
        else
        {
          delta[i] = INFINITY;
        }


      for (int i = 0; i < delta.size(); i++)
      {
        if (delta[i] != INFINITY && (main_row_index == -1 || delta[i] < delta[main_row_index]))
          main_row_index = i;
      }

      if (main_row_index == -1)
      {
        printf("Function is not limited\n");
        return std::vector<int>();
      }
      

      // step 4 update table
      basis[main_row_index] = main_column_index;

      std::vector<std::vector<double>> newSimplexTable = simplexTable;

      // step 4.1 update main row
      for (int i = 0; i < b_col + 1; i++)
        newSimplexTable[main_row_index][i] = simplexTable[main_row_index][i] / simplexTable[main_row_index][main_column_index];

      // step 4.2 update other rows
      for (int i = 0; i < c_row + 1; i++)
      {
        if (i == main_row_index)
          continue;

        for (int j = 0; j < b_col + 1; j++)
          newSimplexTable[i][j] = simplexTable[i][j] - simplexTable[i][main_column_index] * newSimplexTable[main_row_index][j];
      }
      simplexTable = newSimplexTable;
    }
    return basis;
  }

  std::vector<int> createBasis(std::vector<std::vector<double>> &simplexTable)
  {
    std::vector<int> basis;
    int A_width = simplexTable[0].size() - 1;
    int A_height = simplexTable.size() - 1;

    while (basis.size() != A_height)
    {
      int chosen_line = -1;
      int chosen_col = -1;
      
      for (int col = 0; col < A_width; col++)
      {
        if (isIn(col, basis))
          continue;
        int chosen = -1;

        if (simplexTable[basis.size()][col] != 0)
          chosen = basis.size();

        if (chosen == -1)
          continue;

        chosen_line = chosen;
        chosen_col = col;
        break;
      }

      if (chosen_line == -1 || chosen_col == -1)
      {
        printf("VOVAS OBOSRALSYA AGAIN");
      }

      basis.push_back(chosen_col);

      double div = simplexTable[chosen_line][chosen_col];
      // this line
      for (int i = 0; i < simplexTable[chosen_line].size(); i++)
        simplexTable[chosen_line][i] /= div;

      // other lines
      for (int i = 0; i < A_height; i++)
      {
        if (i == chosen_line)
          continue;

        double coeff = simplexTable[i][chosen_col];
        for (int j = 0; j < simplexTable[i].size(); j++)
        {
          simplexTable[i][j] -= coeff * simplexTable[chosen_line][j];
        }
      }
    }

    //fix b
    bool must_fix = false;
    for (int i = 0; i < A_height; i++)
    {
      if (simplexTable[i][simplexTable[i].size() - 1] < 0) // if b[i] < 0 we must fix it
      {
        must_fix = true;
        break;
      }
    }

    while (must_fix)
    {
      int max_b_index = 0;
      int max_el_index = 0;

      // find max abs(b)
      for (int i = 0; i < A_height; i++)
      {
        if (simplexTable[i][simplexTable[i].size() - 1] < simplexTable[max_b_index][simplexTable[i].size() - 1])
        {
          max_b_index = i;
        }
      }

      for (int i = 0; i < A_width; i++)
      {
        if (simplexTable[max_b_index][i] < simplexTable[max_b_index][max_el_index])
        {
          max_el_index = i;
        }
      }

      // update matr
      double div = simplexTable[max_b_index][max_el_index];

      for (int i = 0; i < simplexTable[max_b_index].size(); i++)
        simplexTable[max_b_index][i] /= div;

      for (int k = 0; k < A_height; k++)
      {
        if (k == max_b_index)
          continue;

        double coeff = simplexTable[k][max_el_index];
        for (int i = 0; i < simplexTable[k].size(); i++)
          simplexTable[k][i] -= coeff * simplexTable[max_b_index][i];
      }

      //change basis
      basis[max_b_index] = max_el_index;

      // update state
      must_fix = false;
      for (int i = 0; i < A_height; i++)
      {
        if (simplexTable[i][simplexTable[i].size() - 1] < 0) // if b[i] < 0 we must fix it
        {
          must_fix = true;
        }
      }
    }
    
    // fix C
    for (int i = 0; i < basis.size(); i++)
      if (simplexTable[A_height][basis[i]] != 0)
      {
        double val = simplexTable[A_height][basis[i]];
        for (int j = 0; j < A_width; j++)
          simplexTable[A_height][j] -= simplexTable[i][j] * val;
      }

    return basis;
  }

  vec tableSimplexMethod(void) {
    ld _A[100][100] = {0};
    ld _B[100] = {0};
    ld _C[100] = {0};
    
    // fill A
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        _A[i][j] = A[i][j];

    // fill b
    for (int i = 0; i < M; i++)
      _B[i] = b[i];

    // fill c
    for (int i = 0; i < N; i++)
      _C[i] = -c[i];
    
    vvd newA(M);
    vd newb(_B, _B + M);
    vd newc(_C, _C + N);
    for (int i = 0; i < M; i++) newA[i] = vd(_A[i], _A[i] + N);

    // use someone simplex
    LPSolver solver(newA, newb, newc);
    vd x;
    ld value = solver.Solve(x);

    vec result(N);

    for (int j = 0; j < N; j++)
      result[j] = x[j];

    return result;

#if 0 // HA-HA our simplex have a bug!!!
    std::vector<std::vector<double>> simplexTable;
    std::vector<int> basis;
    vec result(N);

    int b_col = N;
    int c_row = M;

    // FILL SIMPLEX TABLE
    simplexTable.resize(M + 1);
    for (int i = 0; i < M + 1; i++)
      simplexTable[i].resize(N + 1);

    // fill A
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        simplexTable[i][j] = A[i][j];

    // fill b
    for (int i = 0; i < M; i++)
      simplexTable[i][b_col] = b[i];

    // fill c
    for (int i = 0; i < N; i++)
      simplexTable[c_row][i] = c[i];

    // fill f
    simplexTable[M][N] = 0;

    basis = createBasis(simplexTable);

    auto new_basis = doSimplex(simplexTable, basis);
    if (new_basis.size() == 0)
      return vec(N);

    basis = new_basis;

    // send 
    for (int i = 0; i < result.size(); i++)
    {
      int k = -1;

      for (int j = 0; j < basis.size(); j++)
        if (basis[j] == i)
        {
          k = j;
          break;
        }
      if (k != -1)
        result[i] = simplexTable[k][b_col];
      else
        result[i] = 0;
    }

    return result;
#endif
  }


  //function to get correct answer
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

  vec retrieveCorrectAnswerFromDual(vec answer) {
    vec c_true(cOriginalSize);

    //Fill indexies that we didn't change
    int i = 0;
    for (i = 0; i < (cOriginalSize - numOfUnsigned); i++)
        c_true[i] = answer[i];

    //Fill indexies that we did change
    for (int j = answer.size() - numOfUnsigned; j < answer.size(); j++, i++)
      c_true[i] = answer[i] - answer[j];

    //sdvig
    int k = A.sizeW() - cOriginalSize - numOfUnsigned;
    for(int i = 0; i < k; i++){
      double a = c_true[c_true.size() - 1];
      for (int j = c_true.size() - 1; j > 0; j--)
        c_true[j] = c_true[j - 1];
      c_true[0] = a;
    }

    //return
    return c_true;
  }

};

#endif