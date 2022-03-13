#include "TaskLoader.h"

int main(void) {
 
  //parse to task from file
  TaskLoader t;
  // taken from https://www.matburo.ru/Examples/Files/Transport1.pdf
  TransportTask task = t.load("task.txt");
  

  auto q = task.PotentialMethod();
  
  // Answer is here
  //t.retrieveCorrectAnswer(task.x);
  return 0;

}