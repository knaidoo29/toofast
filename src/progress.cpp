#include "../include/progress.h"

#include <iostream>
#include <string>

void progress::print_line(char symbol){
  for(int i = 0; i < 2*length; i++){
    if(i+1 == 2*length){
      cout << symbol << endl;
    }
    else{
      cout << symbol;
    }
  }
}

void progress::print_title(string title, char symbol, int level){
  int title_length = title.length();
  if(title_length%2 != 0){
    title += " ";
  }
  for(int i = 0; i < length - title.length()/2 - 1; i++){
    if(i < 4*(level - 1)){
      cout << " ";
    }
    else{
      cout << symbol;
    }
  }
  cout << " " << title << " ";
  for(int i = 0; i < length - title.length()/2 - 1; i++){
    if(i+1 == length - title.length()/2 - 1){
      if(i+1 < length - title.length()/2 - 1 - 4*(level - 1)){
        cout << symbol << endl;
      }
      else{
        cout << " " << endl;
      }
    }
    else{
      if(i+1 < length - title.length()/2 - 1 - 4*(level - 1)){
        cout << symbol;
      }
      else{
        cout << " ";
      }
    }
  }
}

void progress::start(string title){
  cout << "\n";
  progress::print_line('=');
  progress::print_title(title, '#', 1);
  progress::print_line('-');
  cout << "\n";
}

void progress::end(){
  string end_title = "DONE";
  cout << "\n";
  progress::print_line('-');
  progress::print_title(end_title, '#', 1);
  progress::print_line('=');
  cout << "\n";
}

void progress::start_process(string process_name, int level){
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
      symbol = s4;
      break;
  }
  cout << "\n";
  progress::print_title(process_name, symbol, level);
  cout << "\n";
}

void progress::end_process(int level){
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
      symbol = s4;
      break;
  }
  cout << "\n";
  progress::print_title(end_process, symbol, level);
  cout << "\n";
}
