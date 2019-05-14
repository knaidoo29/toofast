#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "toofast.h"

using namespace std;

// test read ascii table

int main(){

  string filename = "levy_flight_100.txt";

  vector<float>position;

  int row_length, column_length;

  row_length = read_ascii_table(filename, position);
  column_length = position.size()/row_length;

  cout << row_length << endl;
  cout << column_length << endl;

  vector<float>x;
  vector<float>y;
  vector<float>z;

  extract_from_table(position, 0, row_length, x);
  extract_from_table(position, 1, row_length, y);
  extract_from_table(position, 2, row_length, z);

}
