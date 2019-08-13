//#ifndef UTIL_PROGRESS_H
//#define UTIL_PROGRESS_H

#include <iostream>

using namespace std;

template<typename TYPE>
void progress_bar(TYPE i, TYPE size, bool include_index = false){
  /* Prints a dynamic progress bar in the terminal.

  Parameters
  ----------

  type i:
    Index of for-loop to indicate how much progress has been made.
  type size:
    Full length of for-loop.
  bool include_index:
    Outputs the progression of the index itself as well as the progress bar.
  */
  TYPE bar_length = 50;
  double j, sizej;
  if(i != 0){
    if((bar_length*i)/size != (bar_length*(i+1))/size or i==size){
      cout.flush();
      cout << "\r\t|";
      for(int j = 0; j < (bar_length*(i+1))/size; j++){
        cout << "#";
      }
      for(int j = (bar_length*(i+1))/size; j < bar_length; j++){
        cout << "-";
      }
      cout << "|  ";
      j = (double) i+1;
      sizej = (double) size;
      cout << 100.*(j / sizej) << "%";
      if(include_index == true){
        cout << "  [" << i+1 << "/" << size << "]";
      }
      if(i+1 == size){
        cout << "  " << endl;
      }
    }
  }
}

class progress{
  /* A class for organising and printing c++ script progression. */
  public:
    progress(){
      length = 50;
      s1 = '-';
      s2 = '*';
      s3 = '+';
      s4 = '~';
    }
    void start(string title);
    void end();
    void start_process(string process_name, int level);
    void end_process(int level);
  private:
    int length;
    char s1, s2, s3, s4;
    void print_line(char symbol);
    void print_title(string title, char symbol);
};

//#endif