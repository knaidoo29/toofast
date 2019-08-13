#include "../include/twofast.h"
#include <iostream>
#include <string>

using namespace std;

void progress::print_line(char symbol){
  /* Prints line of a symbol. 
  
  Parameters
  ----------
  char symbol:
    A symbol to be printed on a line.
  */
  for(int i = 0; i < 2*length; i++){
    if(i+1 == 2*length){
      cout << symbol << endl;
    }
    else{
      cout << symbol;
    }
  }
}

void progress::print_title(string title, char symbol){
  /* Prints a title with a line of 'symbol' on either side. 

  Parameters
  ----------
  string title:
    A title to write to a line.
  char symbol:
    A symbol to be printed on a line.
  */
  int title_length = title.length();
  if(title_length%2 != 0){
    title += " ";
  }
  for(int i = 0; i < length - title.length()/2 - 1; i++){
    cout << symbol;
  }
  cout << " " << title << " ";
  for(int i = 0; i < length - title.length()/2 - 1; i++){
    if(i+1 == length - title.length()/2 - 1){
      cout << symbol << endl;
    }
    else{
      if(i+1 < length - title.length()/2 - 1){
        cout << symbol;
      }
      else{
        cout << " ";
      }
    }
  }
}

void progress::start(string title){
  /* Prints the start of a process.

  Parameters
  ----------
  string title:
    Title of process being started.
  */
  cout << "\n";
  progress::print_line('=');
  progress::print_title(title, '#');
  progress::print_line('-');
  cout << "\n";
}

void progress::end(){
  /* Prints the end of a process.

  Parameters
  ----------
  string title:
    Title of process that is finished.
  */
  string end_title = "DONE";
  cout << "\n";
  progress::print_line('-');
  progress::print_title(end_title, '#');
  progress::print_line('=');
  cout << "\n";
}

void progress::start_process(string process_name, int level){
  /* Prints the start of a process.

  Parameters
  ----------
  string process_name:
    Name of process beginning.
  int level:
    Determines what symbol to use around the title.
  */
  char symbol;
  switch (level) {
    case 1:
      symbol = s1;
      break;
    case 2:
      symbol = s2;
      break;
    case 3:
      symbol = s3;
      break;
    case 4:
      symbol = s4;
      break;
    default:
      symbol = s1;
      break;
  }
  cout << "\n";
  progress::print_title(process_name, symbol);
  cout << "\n";
}

void progress::end_process(int level){
  /* Prints the end of a process.

  Parameters
  ----------
  string process_name:
    Name of process finishing.
  int level:
    Determines what symbol to use around the title.
  */
  string end_process = "DONE";
  char symbol;
  switch (level) {
    case 1:
      symbol = s1;
      break;
    case 2:
      symbol = s2;
      break;
    case 3:
      symbol = s3;
      break;
    case 4:
      symbol = s4;
      break;
    default:
      symbol = s1;
      break;
  }
  cout << "\n";
  progress::print_title(end_process, symbol);
  cout << "\n";
}
