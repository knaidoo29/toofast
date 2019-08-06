#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "toofast.h"

using namespace std;

int main(){

  // ------------------------------------------------------------------------ //

  string filename_data, filename_rand, filename_output;
  int numbins;
  float minimum, maximum;
  bool uselog;

  // Paramfile -- Proxy

  filename_data = "data.txt";
  filename_rand = "randoms.txt";
  filename_output = "xi.txt";

  minimum = 0.01;
  maximum = 75.;
  numbins = 40;

  uselog = true;

  // ------------------------------------------------------------------------ //

  string title = "Testing Two Point Correlation Function";
  progress prog;
  prog.start(title);

  // ------------------------------------------------------------------------ //

  prog.start_process("Read Data", 2);

  vector<float>pos_data;

  int row_length, column_length;

  row_length = read_ascii_table(filename_data, pos_data);
  column_length = pos_data.size()/row_length;

  vector<float>x_data;
  vector<float>y_data;
  vector<float>z_data;

  extract_from_table(pos_data, 0, row_length, x_data);
  extract_from_table(pos_data, 1, row_length, y_data);
  extract_from_table(pos_data, 2, row_length, z_data);

  // ------------------------------------------------------------------------ //

  prog.start_process("Read Randoms", 2);

  vector<float>pos_rand;

  row_length = read_ascii_table(filename_rand, pos_rand);
  column_length = pos_rand.size()/row_length;

  vector<float>x_rand;
  vector<float>y_rand;
  vector<float>z_rand;

  extract_from_table(pos_rand, 0, row_length, x_rand);
  extract_from_table(pos_rand, 1, row_length, y_rand);
  extract_from_table(pos_rand, 2, row_length, z_rand);

  // ------------------------------------------------------------------------ //

  prog.start_process("Calculate Xi(r)", 2);

  prog.start_process("Get DD", 3);

  float dd[numbins];

  get_dd(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data);

  print_array(dd, numbins);
  /*
  // ------------------------------------------------------------------------ //

  prog.start_process("Get RR", 3);

  float rr[numbins];

  get_dd(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand);

  // ------------------------------------------------------------------------ //

  prog.start_process("Get DR", 3);

  float dr[numbins];

  get_dr(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
    x_rand, y_rand, z_rand);
  */
  // ------------------------------------------------------------------------ //

  prog.start_process("Print xi(r)", 3);

  float r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);
  /*
  float xi[numbins];

  get_xi(dd, dr, rr, numbins, x_data.size(), x_rand.size(), xi);

  cout << '\t' << '\t' << "r" << '\t' << "Xi(r)" << endl;

  for(int i = 0; i < numbins; i++){
    cout << '\t' << '\t' <<  r[i] << '\t' << xi[i] << endl;
  }
  */
  // ------------------------------------------------------------------------ //

  prog.start_process("Save data", 2);

  Writer wr;
  wr.add2header("Levy flight Xi estimation");
  wr.store(r, numbins, "r");
  /*wr.store(dd, numbins, "DD");
  wr.store(rr, numbins, "RR");
  wr.store(dr, numbins, "DR");
  wr.store(xi, numbins, "xi");
  wr.write2file(filename_output);
  */
  prog.end();

  return 0;
}
