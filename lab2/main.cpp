#include "TaskLoader.h"

int main(void) {
 
  //parse to task from file
  TaskLoader t;
  TransportTask task = t.load("14.txt");

  //indexies that the same for potentials fo extreme point
  //HOOOORAY
  //auto v = task.BuildPotentials({ {0,0}, {1,0}, {2,1}, {2,2}, {3,3}, {1,4}, {3,4}, {0,1} });

  ////solve with extreme points method
  //std::cout << "Extreme point method:" << std::endl;
  //try {
  //  t.retrieveCorrectAnswer(task.extremePointMethod()).print();
  //}
  //catch(std::exception &ex){
  //  std::cout << ex.what() << std::endl;
  //}
  std::cout << "\nPotentials method:" << std::endl;
  try {
    t.retrieveCorrectAnswer(task.PotentialMethod()).print();
  }
  catch(std::exception &ex){
    std::cout << ex.what() << std::endl;
  }

  return 0;

}