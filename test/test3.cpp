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

  print_vector(x);

  cout << "- Write to file:" << endl;

  string filename2;
  filename2 = "levy_flight_100_out.txt";
  Writer wr;
  wr.add2header("levy flight positions");
  wr.store(x, column_length, "x");
  wr.store(y, column_length, "y");
  wr.store(z, column_length, "z");
  wr.write2file(filename2);

}
