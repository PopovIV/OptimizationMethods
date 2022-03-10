#include "TaskLoader.h"

int main(void) {
 
  //parse to task from file
  TaskLoader t;
  TransportTask task = t.load("task.txt");
  
  // Answer is here
  //t.retrieveCorrectAnswer(task.x);
  return 0;

}