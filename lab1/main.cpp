#include "TaskLoader.h"

int main(void) {
 
  //parse to canon from file
  TaskLoader t;
  std::pair<LinearProgTask, LinearProgTask> task = t.load("task.txt");
  //solve LPP using two methods
  vec sol(0), sol2(0);
  
  try {
     sol = task.first.tableSimplexMethod();
  }
  catch (std::exception& ex) {
    std::cout << ex.what() << std::endl;
  }
  try{
    sol2 = task.first.extremePointMethod();
  }
  catch (std::exception& ex) {
    std::cout << ex.what() << std::endl;
  }
  //print answers
  std::cout << "Solution of canon form "<<std::endl; 
  std::cout << "Simplex:" << std::endl;
  sol.print();
  std::cout << "Extreme point:" << std::endl;
  sol2.print();
  std::cout << "Solution of common form " << std::endl;
  std::cout << "Simplex:" << std::endl;
  task.first.retrieveCorrectAnswer(sol).print();
  std::cout << "x * C :" <<  sol * task.first.getC() << std::endl;
  std::cout << "Extreme point:" << std::endl;
  if (sol2.size() != 0)
    task.first.retrieveCorrectAnswer(sol2).print();
  std::cout << "x * C :" << sol2 * task.first.getC() << std::endl;
  std::cout << std::endl;
  
  //solve dual LPP using two methods
  try {
     sol = task.second.tableSimplexMethod();
  }
  catch (std::exception& ex) {
    std::cout << ex.what() << std::endl;
  }
  try{
    sol2 = task.second.extremePointMethod();
  }
  catch (std::exception& ex) {
    sol2 = vec(0);
    std::cout << ex.what() << std::endl;
  }
  //print answers
  std::cout << "Dual task" << std::endl;
  std::cout << "Solution of canon form " << std::endl;
  std::cout << "Simplex:" << std::endl;
  sol.print();
  std::cout << "Extreme point:" << std::endl;
  sol2.print();
  std::cout << "Solution of common form " << std::endl;
  std::cout << "Simplex:" << std::endl;
  task.second.retrieveCorrectAnswer(sol).print();
  std::cout << "y * C :" << sol * task.second.getC() << std::endl;
  std::cout << "Extreme point:" << std::endl;
  if (sol2.size() != 0)
    task.second.retrieveCorrectAnswer(sol2).print();
  std::cout << "y * C :" << sol2 * task.second.getC() << std::endl;

  return 0;

}