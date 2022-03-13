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

  //utility functions for extreme point method
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
  vec Gauss(matr A, vec B) {

    //vec B to matrix
    matr b(B.size(), 1);
    for (int i = 0; i < B.size(); i++)
      b[i][0] = B[i];
    A.concatinateRight(b);
    int N = A.sizeH();
    // Triangularization
    for (int k = 0; k < N; k++) {
      int i_max = k;
      int v_max = fabs(A[i_max][k]);

      for (int i = k + 1; i < N; i++)
        if (fabs(A[i][k]) > v_max)
          v_max = A[i][k], i_max = i;
      if (i_max != k)
        for (int t = 0; t <= N; t++) {
          double tmp = A[k][t];
          A[k][t] = A[i_max][t];
          A[i_max][t] = tmp;
        }
      for (int i = k + 1; i < N; i++) {
        double f = A[i][k] / A[k][k];
        for (int j = k + 1; j <= N; j++)
          A[i][j] -= A[k][j] * f;
        A[i][k] = 0;
      }
    }
    // Resolution
    vec x(N);
    for (int i = N - 1; i >= 0; i--) {
      x[i] = A[i][N];
      for (int j = i + 1; j < N; j++)
        x[i] -= A[i][j] * x[j];
      x[i] = x[i] / A[i][i];
    }

    return x;

  }

public:
  // constructor from Common Form
  TransportTask(
    matr c,    // size is [N, M]  and sum(sum(c[i][j] * x[i][j])) -> min
    vec a,  // size is N
    vec b   // size is M
  ) : C(c.sizeH(), c.sizeW()), a(a), b(b)
  {
    C = c;

    NorthWestMethod(a, b);
  }

  // for first solution
  std::pair<matr, std::vector<std::pair<int, int>>> NorthWestMethod(vec A, vec B)
  {
    matr X(A.size(), B.size());
    int nw_cell_x = 0;
    int nw_cell_y = 0;
    std::vector<std::pair<int, int>> basis;
    while (nw_cell_y < A.size() && nw_cell_x < B.size())
    {
      basis.push_back(std::make_pair(nw_cell_y, nw_cell_x));
      double col_amount = B[nw_cell_x];
      double row_amount = A[nw_cell_y];

      for (int i = 0; i < nw_cell_x; i++)
        row_amount -= X[nw_cell_y][i];

      for (int i = 0; i < nw_cell_y; i++)
        col_amount -= X[i][nw_cell_x];

      if (col_amount < row_amount)
      {
        X[nw_cell_y][nw_cell_x] = col_amount;
        nw_cell_x++;
      }
      else
      {
        X[nw_cell_y][nw_cell_x] = row_amount;
        nw_cell_y++;
      }
    }
    if (nw_cell_y < A.size() - 1 || nw_cell_x < B.size() - 1)
      printf("Vovas made a huge mistake");
    return std::make_pair(X, basis);
  }

  std::pair<vec, vec> BuildPotentials(std::vector<std::pair<int, int>> basis)
  {
    vec U(a.size());
    vec V(b.size());
    std::vector<int> filled_U;
    std::vector<int> filled_V;

    filled_U.push_back(basis[0].second);
    // 3.1) let some value became 0
    U[basis[0].second] = 0;
    bool is_fill_U = false;
    bool is_done = false;
    // 3.1) repeat until we done
    while (!is_done)
    {
      is_done = true;
      bool is_changed = false;
      int remember_x = -1;

      // 3.1) let some value became 0
      for (int i = 0; i < basis.size(); i++)
      {
        int cell_y = basis[i].first;
        int cell_x = basis[i].second;

        bool is_x = false;
        bool is_y = false;

        for (int j = 0; j < filled_U.size(); j++)
          if (cell_y == filled_U[j])
          {
            is_y = true;
            break;
          }

        for (int j = 0; j < filled_V.size(); j++)
          if (cell_x == filled_V[j])
          {
            is_x = true;
            break;
          }

        if (is_x && !is_y)
        {
          U[cell_y] = C[cell_y][cell_x] - V[cell_x];
          filled_U.push_back(cell_y);
          is_changed = true;
        }
        else if (!is_x && is_y)
        {
          V[cell_x] = C[cell_y][cell_x] - U[cell_y];
          filled_V.push_back(cell_x);
          is_changed = true;
        }
        else if (!is_x && !is_y)
        {
          remember_x = cell_x;
        }

        if (!(is_x && is_y))
        {
          is_done = false;
        }
      }
      // if we iterate and something do not changed
      if (is_changed == false && is_done == false)
      {
        V[remember_x] = 0;
        filled_V.push_back(remember_x);
      }
    }
    return std::make_pair(U, V);
  }

  // check all possible cycles and find shortest one
  // Using DFS. Aghhhh
  unsigned int FindCycle(
    std::vector<std::vector<int>> vertices_for_x,
    std::vector<std::vector<int>> vertices_for_y,
    std::pair<int, int> start,
    std::vector<std::pair<int, int>>& cur_path,
    bool go_vertical = false)
  {
    if (cur_path.size() == 0)
    {
      cur_path.push_back(start);
    }

    // if we added some shit
    for (int i = 0; i < cur_path.size() - 1; i++)
      if (cur_path[i] == cur_path[cur_path.size() - 1])
        return 0;

    // end of recursion
    if ((cur_path.size() > 2) &&
      (cur_path[cur_path.size() - 1].first == start.first ||
        cur_path[cur_path.size() - 1].second == start.second))
    {
      return cur_path.size();
    }

    int y_to_check = cur_path[cur_path.size() - 1].first;
    int x_to_check = cur_path[cur_path.size() - 1].second;

    unsigned int best = -1;// mean largest number possible;
    std::vector<std::pair<int, int>> best_res;

    // check all horisontal
    if (!go_vertical)
    {
      for (int i = 0; i < vertices_for_y[y_to_check].size(); i++)
      {
        cur_path.push_back(std::make_pair(y_to_check, vertices_for_y[y_to_check][i]));

        auto res = FindCycle(vertices_for_x, vertices_for_y, start, cur_path, true);

        if (res != 0 && cur_path.size() < best)
        {
          best_res = cur_path;
          best = res;
        }

        cur_path.pop_back();
      }
    }
    else // check all vertical
    {
      for (int i = 0; i < vertices_for_x[x_to_check].size(); i++)
      {
        cur_path.push_back(std::make_pair(vertices_for_x[x_to_check][i], x_to_check));

        auto res = FindCycle(vertices_for_x, vertices_for_y, start, cur_path, false);

        if (res != 0 && cur_path.size() < best)
        {
          best_res = cur_path;
          best = res;
        }

        cur_path.pop_back();
      }
    }
    if (best != -1)
      cur_path = best_res;
    return best_res.size();
  }

  matr PotentialMethod(void) {

    // 1) Find first solution
    auto Xb = NorthWestMethod(a, b);

    matr X = Xb.first;

    // 2) find basis
    std::vector<std::pair<int, int>> basis = Xb.second;
    std::vector<std::vector<int>> basis_y_for_x;
    std::vector<std::vector<int>> basis_x_for_y;

    basis_y_for_x.resize(b.size());
    basis_x_for_y.resize(a.size());

    for (int k = 0; k < basis.size(); k++)
    {
      int i = basis[k].first;
      int j = basis[k].second;
      basis_x_for_y[i].push_back(j);
      basis_y_for_x[j].push_back(i);
    }

    // 4) cycle of DEATH begins
    while (true)
    {
      // 3) build potentials
      auto uv = BuildPotentials(basis);
      vec U = uv.first;
      vec V = uv.second;

      // check optimality
      int cur_min_x = 0;
      int cur_min_y = 0;
      double value = C[0][0] - U[0] - V[0];

      for (int i = 0; i < C.sizeH(); i++)
        for (int j = 0; j < C.sizeW(); j++)
          if (C[i][j] - U[i] - V[j] < value)
          {
            cur_min_x = j;
            cur_min_y = i;
            value = C[i][j] - U[i] - V[j];
          }
      if (value >= 0)
        break;

      // if not optimal, we need to rotate something

      // 5) find cycle
      std::vector<std::pair<int, int>> cycle;
      int res = FindCycle(basis_y_for_x, basis_x_for_y, std::make_pair(cur_min_y, cur_min_x), cycle);

      if (res < 3)
        throw std::exception("Cannot find cycle");

      std::pair<int, int> min_cell = cycle[1];
      // 6) find in cycle minimum on odd (because we start with 0) position
      for (int i = 1; i < cycle.size(); i += 2)
        if (X[cycle[i].first][cycle[i].second] < X[min_cell.first][min_cell.second])
          min_cell = cycle[i];
      double val = X[min_cell.first][min_cell.second];
      // 7) change numbers in solution
      for (int i = 0; i < cycle.size(); i++)
        X[cycle[i].first][cycle[i].second] += (((i + 1) % 2) * 2 - 1) * val; // cycle[0] += val, cycle[1] -= val and so on

      // 8) remove zero cell from basis
      for (int i = 0; i < basis.size(); i++)
        if (basis[i] == min_cell)
        {
          basis.erase(basis.begin() + i);
          break;
        }

      for (int i = 0; i < basis_x_for_y[min_cell.first].size(); i++)
        if (basis_x_for_y[min_cell.first][i] == min_cell.second)
        {
          basis_x_for_y[min_cell.first].erase(basis_x_for_y[min_cell.first].begin() + i);
          break;
        }

      for (int i = 0; i < basis_y_for_x[min_cell.second].size(); i++)
        if (basis_y_for_x[min_cell.second][i] == min_cell.first)
        {
          basis_y_for_x[min_cell.second].erase(basis_y_for_x[min_cell.second].begin() + i);
          break;
        }
      // 9) add new cell in basis
      basis.push_back(std::make_pair(cur_min_y, cur_min_x));
      basis_x_for_y[cur_min_y].push_back(cur_min_x);
      basis_y_for_x[cur_min_x].push_back(cur_min_y);
    }

    double sum = 0;
    for (int i = 0; i < C.sizeH(); i++)
      for (int j = 0; j < C.sizeW(); j++)
        sum += X[i][j] * C[i][j];
    std::cout << "Total sum: " << sum << std::endl;

    return X;
  }

  //Extreme point method
  matr extremePointMethod(void) {

    //fistly transform transport task to linear program task
    //free column
    vec f(a);
    f.concatinate(b);
    //matrix and vector function to minize
    matr A(f.size(), C.sizeH() * C.sizeW());
    vec z(0);
    for (int i = 0, A_j = 0; i < C.sizeH(); i++, A_j += C.sizeW())
      for (int j = 0; j < C.sizeW(); j++) {
        z.append(C[i][j]);
        A[i][j + A_j] = 1;
        A[C.sizeH() + j][j + A_j] = 1;
      }
    A.deleteLastRow();
    f.pop();

    //method itself
    std::vector<int> setOfIndexies(A.sizeW());
    for (int i = 0; i < setOfIndexies.size(); i++)
      setOfIndexies[i] = i;
    vec solution(0);
    double minOfFunction = std::numeric_limits<double>::max();
    std::vector<std::vector<int>> vectorOfIndexies = combine(A.sizeW(), A.sizeH());
    for (auto& indexies : vectorOfIndexies) {//checks every possible combination of columns
      //get matrix for linear system
      matr subMatr(0, 0);
      try {
        subMatr = A.getSubMatrix(indexies);
      }
      catch (std::exception& ex) {
        throw ex;
      }
      //check det
      //try {
      //  if (subMatr.determinant() == 0)
      //    continue;
      //}
      //catch (std::exception& ex) {
      //  throw ex;
      //}
      //solve linear system

      vec systemSolution = Gauss(subMatr, f);

      //check that solution is >= 0
      bool flag = true;
      for (int i = 0; i < systemSolution.size(); i++)
        if (systemSolution[i] < 0 || std::isnan(systemSolution[i])) {
          if (fabs(systemSolution[i]) == 0) {
            systemSolution[i] = 0;
            break;
          }
          flag = false;
          break;
        }//bad practice, ask Vova how to do better
      if (!flag)
        continue;
      //add zeros to solution 
      vec newSolution(z.size());
      for (int i = 0; i < indexies.size(); i++)
        newSolution[indexies[i]] = systemSolution[i];

      //check if that vector is maybe better solution
      double functionVal = 0;
      for (int i = 0; i < newSolution.size(); i++)
        functionVal += newSolution[i] * z[i];
      if (functionVal < minOfFunction) {
        minOfFunction = functionVal;
        solution = newSolution;
      }
    }

    if (solution.size() == 0)
      throw std::exception("No solutions found in extreme point");

    std::cout << "Total sum: " << solution * z << std::endl;

    matr sol(C.sizeH(), C.sizeW());
    for (int i = 0, k = 0; i < C.sizeH(); i++)
      for (int j = 0; j < C.sizeW(); j++)
        sol[i][j] = solution[k++];
    return sol;

  }

};

#endif