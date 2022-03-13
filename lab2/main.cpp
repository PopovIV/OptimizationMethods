#include "TaskLoader.h"

int main(void) {
 
  //parse to task from file
  TaskLoader t;
  TransportTask task = t.load("14.txt");

  //solve with extreme points method
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