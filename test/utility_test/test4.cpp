#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "toofast.h"

using namespace std;

// test progress

int main(){
  string title;

  title = "Testing Progress String";

  progress prog;
  prog.start(title);
  prog.start_process("Nested 1", 1);
  prog.start_process("Nested 2", 2);
  prog.start_process("Nested 3", 3);
  prog.start_process("Nested 4", 4);
  prog.end_process(4);
  prog.end_process(3);
  prog.end_process(2);
  prog.end_process(1);

  prog.start_process("Subtitle 1", 1);
  prog.end_process(1);
  prog.start_process("Subtitle 2", 2);
  prog.end_process(2);
  prog.start_process("Subtitle 3", 3);
  prog.end_process(3);
  prog.start_process("Subtitle 4", 4);
  prog.end_process(4);

  prog.end();

  return 0;

}
