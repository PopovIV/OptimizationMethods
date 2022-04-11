#ifndef MATH_H_INCLUDED__
#define MATH_H_INCLUDED__

#include <vector>
#include <iostream>//for debug

class vec
{
private:
  std::vector<double> A;
public:
  vec(std::initializer_list<double> data)
  {
    A.resize(data.size());
    for (int i = 0; i < A.size(); i++)
      A[i] = *(data.begin() + i);
  }

  vec(int size)
  {
    A.resize(size);
    for (int i = 0; i < A.size(); i++)
      A[i] = 0;
  }

  double & operator[](int i)
  {
    return A[i];
  }

  int size()
  {
    return A.size();
  }

  double operator*(vec other)
  {
    double res = 0;
    if (A.size() != other.A.size())
      throw std::exception("ERROR: incorrect size of vector in mul operator");
    for (int i = 0; i < A.size(); i++)
      res += A[i] * other[i];

    return res;
  }

  vec operator-(vec other)
  {
    vec res(A.size());
    if (A.size() != other.A.size())
      throw std::exception("ERROR: incorrect size of vector in mul operator");
    for (int i = 0; i < A.size(); i++)
      res[i] = A[i] - other[i];

    return res;
  }

  vec operator+(vec other)
  {
    vec res(A.size());
    if (A.size() != other.A.size())
      throw std::exception("ERROR: incorrect size of vector in mul operator");
    for (int i = 0; i < A.size(); i++)
      res[i] = A[i] + other[i];

    return res;
  }

  vec operator*(double other)
  {
    vec res(A.size());
    for (int i = 0; i < A.size(); i++)
      res[i] = A[i] * other;

    return res;
  }

  void append(double q)
  {
    A.push_back(q);
  }

  void concatinate(vec q)
  {
    A.insert(A.end(), q.A.begin(), q.A.end());
  }

  double lenSquared(void)
  {
    double sum_of_elems = 0;
    for (auto& n : A)
      sum_of_elems += n * n;
    return sum_of_elems;
  }

  double len(void)
  {
    return sqrt(lenSquared());
  }

  //for debug
  void print(void) {
    for (auto f : A)
      std::cout << f << " ";
    std::cout << std::endl;
  }

};

class matr
{
private:
  std::vector<std::vector<double>> A;

  // Получение матрицы без i-й строки и j-го столбца
  matr getMatr(int i, int j) {

    int di = 0, dj = 0;
    matr res(sizeH() - 1, sizeW() - 1);
    for (int ki = 0; ki < sizeW() - 1; ki++) {
      if (ki == i)
        di = 1;
      dj = 0;
      for (int kj = 0; kj < sizeW() - 1; kj++) {
        if (kj == j)
          dj = 1;
        res[ki][kj] = A[ki + di][kj + dj];
      }
    }

    return res;

  };

public:
  matr(std::initializer_list<std::initializer_list<double>> data)
  {
    A.resize(data.size());
    for (int i = 0; i < A.size(); i++)
      A[i].resize((data.begin() + i)->size());

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A[i].size(); j++)
        A[i][j] = *((*(data.begin() + i)).begin() + j);
  }

  matr(int sizeH, int sizeW)
  {
    if (sizeH < 0)
      sizeH = 0;
    A.resize(sizeH);
    for (int i = 0; i < A.size(); i++)
      A[i].resize(sizeW);

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A[i].size(); j++)
        A[i][j] = 0;
  }

  //Is it bad? Don't think so
  matr(std::vector<std::vector<double>> matrix) {
    A = matrix;
  }

  int sizeH()
  {
    return A.size();
  }

  int sizeW()
  {
    return A.size() > 0 ? A[0].size() : 0;
  }

  std::vector<double>& operator[](int i)
  {
    return A[i];
  }

  vec operator*(vec other)
  {
    vec res(this->sizeH());
    if (this->sizeW() != other.size())
      throw std::exception("ERROR: incorrect size of vector in mul operator");
    for (int i = 0; i < this->sizeH(); i++)
      for (int j = 0; j < this->sizeW(); j++)
        res[i] += A[i][j] * other[j];
    return res;
  }

  matr operator*(matr other)
  {
    matr res(this->sizeH(), other.sizeW());
    if (this->sizeW() != other.sizeH())
      throw std::exception("ERROR: incorrect size of matrix in mul operator");
    for (int i = 0; i < this->sizeH(); i++)
      for (int j = 0; j < other.sizeW(); j++)
        for (int k = 0; k < this->sizeW(); k++)
          res[i][j] += A[i][k] * other[k][j];
    return res;
  }

  matr operator-(matr other)
  {
    matr res(this->sizeH(), this->sizeW());
    if (this->sizeW() != other.sizeW() || this->sizeH() != other.sizeH())
      throw std::exception("ERROR: incorrect size of matrix in mul operator");
    for (int i = 0; i < this->sizeH(); i++)
      for (int j = 0; j < this->sizeW(); j++)
         res[i][j] = A[i][j] - other[i][j];
    return res;
  }

  matr operator+(matr other)
  {
    matr res(this->sizeH(), this->sizeW());
    if (this->sizeW() != other.sizeW() || this->sizeH() != other.sizeH())
      throw std::exception("ERROR: incorrect size of matrix in mul operator");
    for (int i = 0; i < this->sizeH(); i++)
      for (int j = 0; j < this->sizeW(); j++)
        res[i][j] = A[i][j] + other[i][j];
    return res;
  }
  matr operator*(double other)
  {
    matr res(this->sizeH(), this->sizeW());
    for (int i = 0; i < this->sizeH(); i++)
      for (int j = 0; j < this->sizeW(); j++)
          res[i][j] = A[i][j] * other;
    return res;
  }

  void concatinateRight(matr other)
  {
    if (this->sizeH() != other.sizeH())
      throw std::exception("ERROR: incorrect size of matrix in concatinateRight function");
    for (int i = 0; i < this->sizeH(); i++)
      A[i].insert(A[i].end(), other.A[i].begin(), other.A[i].end());
  }

  void concatinateDown(matr other)
  {
    if (other.sizeW() == 0)
      return;
    if (this->sizeH() == 0) {
      this->A = other.A;
      return;
    }
    if (this->sizeW() != other.sizeW())
      throw std::exception("ERROR: incorrect size of matrix in concatinateDown function");
    A.insert(A.end(), other.A.begin(), other.A.end());
  }
  //my additions for math.h
  void transpose(void) {

    matr tmp(sizeW(), sizeH());
    for (int i = 0; i < tmp.sizeH(); i++)
      for (int j = 0; j < tmp.sizeW(); j++)
        tmp[i][j] = A[j][i];
    A = tmp.A;
  }

  //method to get submatrix by vector of indexies
  matr getSubMatrix(std::vector<int> indexies) {

    if (sizeW() < indexies.size())
      throw std::exception("ERROR: incorrect size of matrix in getSubMatrix function");

    matr subMatr(sizeH(), indexies.size());
    for (int i = 0; i < subMatr.sizeH(); i++)
      for (int j = 0; j < subMatr.sizeW(); j++)
        subMatr[i][j] = A[i][indexies[j]];

    return subMatr;

  }

  //find determinant of square matrix only
  double determinant(void) {

    if (sizeH() != sizeW() || sizeH() == 0)
      throw std::exception("ERROR: incorrect size of matrix in determinant function");

    switch (sizeH()) {
    case 1:
      return A[0][0];
    case 2:
      return A[0][0] * A[1][1] - (A[1][0] * A[0][1]);
    default:
      int k = 1;//(-1)^1
      double result = 0;
      for (int i = 0; i < sizeH(); i++) {
        result += k * A[i][0] * getMatr(i, 0).determinant();
        k *= -1;
      }
      return result;
    }

  }

  void deleteLastRow(void) {
    A.pop_back();
  }

  void fromVectors(vec v1, vec v2) {
    A.clear();
    
    A.resize(v1.size());
    for (auto &el : A)
      el.resize(v2.size());

    for (int i = 0; i < sizeH(); i++)
      for (int j = 0; j < sizeW(); j++)
        A[i][j] = v1[i] * v2[j];
  }

  //for debug
  void print(void) {
    for (int i = 0; i < A.size(); i++) {
      for (int j = 0; j < A[0].size(); j++)
        std::cout << A[i][j] << " ";
      std::cout << std::endl;
    }
  }

};
#endif