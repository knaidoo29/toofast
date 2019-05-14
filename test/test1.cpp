#include <iostream>
#include "toofast.h"

using namespace std;

// Progress bar test

int main(){

  long int size = 1000000000;

  for(long int i = 0; i < size; i++){
    progress_bar(i, size);
  }

  for(long int i = 0; i < size; i++){
    progress_bar(i, size, true);
  }

}
