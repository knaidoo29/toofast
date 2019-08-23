#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "twofast.h"

using namespace std;

int main(int argc, char** argv){

  progress prog;
  prog.start("TwoFast: 2-Point");

  // READING PARAMFILE --------------------------------------------------------------- //

  prog.start_process("Reading Parameter File", 1);

  string prefix, paramfile, whichparam, data_file, rand_file, out_file;
  bool data_file_found = false, rand_file_found = false, out_file_found = false;

  string mode, calc_DD, calc_DR, calc_RR, calc_xi, uselog_str, minimum_str, maximum_str, numbins_str;
  bool mode_found = false, calc_DD_found = false, calc_DR_found = false, calc_RR_found = false, calc_xi_found = false;
  bool uselog_found = false, minimum_found = false, maximum_found = false, numbins_found = false;

  paramfile = argv[1];

  cout << "Paramfile = " << paramfile << endl;

  whichparam = "data_file";
  extract_param(paramfile, whichparam, &data_file, &data_file_found);
  summarise_param(whichparam, data_file, data_file_found);

  whichparam = "rand_file";
  extract_param(paramfile, whichparam, &rand_file, &rand_file_found);
  summarise_param(whichparam, rand_file, rand_file_found);

  whichparam = "out_file";
  extract_param(paramfile, whichparam, &out_file, &out_file_found);
  summarise_param(whichparam, out_file, out_file_found);

  whichparam = "mode";
  extract_param(paramfile, whichparam, &mode, &mode_found);
  for_each(mode.begin(), mode.end(), to_lowercase);
  summarise_param(whichparam, mode, mode_found);

  whichparam = "calc_DD";
  extract_param(paramfile, whichparam, &calc_DD, &calc_DD_found);
  for_each(calc_DD.begin(), calc_DD.end(), to_lowercase);
  summarise_param(whichparam, calc_DD, calc_DD_found);

  whichparam = "calc_DR";
  extract_param(paramfile, whichparam, &calc_DR, &calc_DR_found);
  for_each(calc_DR.begin(), calc_DR.end(), to_lowercase);
  summarise_param(whichparam, calc_DR, calc_DR_found);

  whichparam = "calc_RR";
  extract_param(paramfile, whichparam, &calc_RR, &calc_RR_found);
  for_each(calc_RR.begin(), calc_RR.end(), to_lowercase);
  summarise_param(whichparam, calc_RR, calc_RR_found);

  whichparam = "calc_xi";
  extract_param(paramfile, whichparam, &calc_xi, &calc_xi_found);
  for_each(calc_xi.begin(), calc_xi.end(), to_lowercase);
  summarise_param(whichparam, calc_xi, calc_xi_found);

  whichparam = "uselog";
  extract_param(paramfile, whichparam, &uselog_str, &uselog_found);
  for_each(uselog_str.begin(), uselog_str.end(), to_lowercase);
  summarise_param(whichparam, uselog_str, uselog_found);

  whichparam = "minimum";
  extract_param(paramfile, whichparam, &minimum_str, &minimum_found);
  summarise_param(whichparam, minimum_str, minimum_found);

  whichparam = "maximum";
  extract_param(paramfile, whichparam, &maximum_str, &maximum_found);
  summarise_param(whichparam, maximum_str, maximum_found);

  whichparam = "numbins";
  extract_param(paramfile, whichparam, &numbins_str, &numbins_found);
  summarise_param(whichparam, numbins_str, numbins_found);

  int numbins;
  float minimum, maximum;
  bool uselog;

  if(uselog_str == "yes"){
    uselog = true;
  }
  else if(uselog_str == "no"){
    uselog = false;
  }
  else{
    cout << "Error: uselog is ill defined -- found in paramfile = " << uselog_str << endl;
    cout << "Should be either yes/no" << endl;
    cout << "Will assume 'no'" << endl;
    uselog = false;
  }

  minimum = stof(minimum_str);
  maximum = stof(maximum_str);
  numbins = stof(numbins_str);

  // --------------------------------------------------------------------------------- //

  prog.start_process("Reading data", 1);

  vector<float>pos_data;

  int row_length, column_length;

  if(data_file_found == true){
    row_length = read_ascii_table(data_file, pos_data);
    column_length = pos_data.size()/row_length;
  }

  vector<float>x_data;
  vector<float>y_data;
  vector<float>z_data;
  vector<float>phi_data;
  vector<float>theta_data;

  if(data_file_found == true){
    if(mode == "2d" || mode == "3d"){
      extract_from_table(pos_data, 0, row_length, x_data);
      extract_from_table(pos_data, 1, row_length, y_data);
      if(mode == "3d"){
        extract_from_table(pos_data, 2, row_length, z_data);
      }
    }
    else if(mode == "tomo"){
      extract_from_table(pos_data, 0, row_length, phi_data);
      extract_from_table(pos_data, 1, row_length, theta_data);
    }
  }

  if((data_file_found == true) && (calc_DD == "yes")){
    cout << "Loaded data file." << endl;
  }

  // --------------------------------------------------------------------------------- //

  vector<float>pos_rand;

  if(rand_file_found == true){
    row_length = read_ascii_table(rand_file, pos_rand);
    column_length = pos_rand.size()/row_length;
  }

  vector<float>x_rand;
  vector<float>y_rand;
  vector<float>z_rand;
  vector<float>phi_rand;
  vector<float>theta_rand;

  if(data_file_found == true){
    if(mode == "2d" || mode == "3d"){
      extract_from_table(pos_rand, 0, row_length, x_rand);
      extract_from_table(pos_rand, 1, row_length, y_rand);
      if(mode == "3d"){
        extract_from_table(pos_rand, 2, row_length, z_rand);
      }
    }
    else if(mode == "tomo"){
      extract_from_table(pos_rand, 0, row_length, phi_rand);
      extract_from_table(pos_rand, 1, row_length, theta_rand);
    }
  }

  if((rand_file_found == true) && ((calc_DR == "yes") || (calc_RR == "yes"))){
    cout << "Loaded randoms file." << endl;
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating DD", 1);

  vector<float>dd;

  for(int i = 0; i < numbins; i++){
    dd.push_back(0.);
  }

  if(calc_DD == "yes"){
    if(mode == "2d"){
      get_dd_2d(dd, minimum, maximum, numbins, uselog, x_data, y_data);
    }
    else if(mode == "3d"){
      get_dd_3d(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data);
    }
    else if(mode == "tomo"){
      get_dd_tomo(dd, minimum, maximum, numbins, uselog, phi_data, theta_data);
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating DR", 1);

  vector<float>dr;

  for(int i = 0; i < numbins; i++){
    dr.push_back(0.);
  }

  if(calc_DR == "yes"){
    if(mode == "2d"){
      get_dr_2d(dr, minimum, maximum, numbins, uselog, x_data, y_data, x_rand, y_rand);
    }
    else if(mode == "3d"){
      get_dr_3d(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
        x_rand, y_rand, z_rand);
    }
    else if(mode == "tomo"){
      get_dr_tomo(dr, minimum, maximum, numbins, uselog, phi_data, theta_data, phi_rand, theta_rand);
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating RR", 1);

  vector<float>rr;

  for(int i = 0; i < numbins; i++){
    rr.push_back(0.);
  }

  if(calc_RR == "yes"){
    if(mode == "2d"){
      get_dd_2d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand);
    }
    else if(mode == "3d"){
      get_dd_3d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand);
    }
    else if(mode == "tomo"){
      get_dd_tomo(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand);
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Save data", 1);

  float r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);

  Writer wr;
  wr.add2header("Levy flight Xi estimation");
  wr.store(r, numbins, "r");
  if(calc_DD == "yes"){
    wr.store(dd, numbins, "DD");
  }
  if(calc_DR == "yes"){
    wr.store(dr, numbins, "DR");
  }
  if(calc_RR == "yes"){
    wr.store(rr, numbins, "RR");
  }
  if(calc_xi == "yes"){
    vector<float>xi;
    for(int i = 0; i < numbins; i++){
      xi.push_back(0.);
    }
    get_xi(x_data.size(), x_rand.size(), dd, dr, rr, xi);
    wr.store(xi, numbins, "xi");
  }
  wr.write2file(out_file);

  prog.end();

  return 0;
}
