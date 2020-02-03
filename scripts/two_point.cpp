#include <iostream>
#include <algorithm>
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

  string mode, calc_DD, calc_DR, calc_RR, calc_xi, uselog_str, useweight_str, minimum_str, maximum_str, numbins_str, mu_bins_str;
  bool mode_found = false, calc_DD_found = false, calc_DR_found = false, calc_RR_found = false, calc_xi_found = false;
  bool useweight_found = false, uselog_found = false, minimum_found = false, maximum_found = false, numbins_found = false, mu_bins_found = false;

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

  whichparam = "useweight";
  extract_param(paramfile, whichparam, &useweight_str, &useweight_found);
  for_each(useweight_str.begin(), useweight_str.end(), to_lowercase);
  summarise_param(whichparam, useweight_str, useweight_found);

  whichparam = "minimum";
  extract_param(paramfile, whichparam, &minimum_str, &minimum_found);
  summarise_param(whichparam, minimum_str, minimum_found);

  whichparam = "maximum";
  extract_param(paramfile, whichparam, &maximum_str, &maximum_found);
  summarise_param(whichparam, maximum_str, maximum_found);

  whichparam = "numbins";
  extract_param(paramfile, whichparam, &numbins_str, &numbins_found);
  summarise_param(whichparam, numbins_str, numbins_found);

  whichparam = "mu_bins";
  extract_param(paramfile, whichparam, &mu_bins_str, &mu_bins_found);
  summarise_param(whichparam, mu_bins_str, mu_bins_found);

  int numbins, mu_bins;
  double minimum, maximum;
  bool uselog, useweight;

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

  if(useweight_str == "yes"){
    useweight = true;
  }
  else if(useweight_str == "no"){
    useweight = false;
  }
  else{
    cout << "Error: useweight is ill defined -- found in paramfile = " << useweight_str << endl;
    cout << "Should be either yes/no" << endl;
    cout << "Will assume 'no'" << endl;
    useweight = false;
  }

  minimum = stof(minimum_str);
  maximum = stof(maximum_str);
  numbins = stoi(numbins_str);
  mu_bins = stoi(mu_bins_str);

  // --------------------------------------------------------------------------------- //

  prog.start_process("Reading data", 1);

  vector<double>pos_data;

  int row_length, column_length;

  if(data_file_found == true){
    row_length = read_ascii_table(data_file, pos_data);
    column_length = pos_data.size()/row_length;
  }

  vector<double>x_data;
  vector<double>y_data;
  vector<double>z_data;
  vector<double>phi_data;
  vector<double>theta_data;
  vector<double>w_data;

  if(data_file_found == true){
    if(mode == "2d" || mode == "3d" || mode == "poly"){
      extract_from_table(pos_data, 0, row_length, x_data);
      extract_from_table(pos_data, 1, row_length, y_data);
      if(mode == "2d" && useweight == true){
        extract_from_table(pos_data, 2, row_length, w_data);
      }
      if(mode == "3d" || mode == "poly"){
        extract_from_table(pos_data, 2, row_length, z_data);
        if(mode == "3d" && useweight == true){
          extract_from_table(pos_data, 3, row_length, w_data);
        }
        if(mode == "poly" && useweight == true){
          extract_from_table(pos_data, 3, row_length, w_data);
        }
      }
    }
    else if(mode == "tomo"){
      extract_from_table(pos_data, 0, row_length, phi_data);
      extract_from_table(pos_data, 1, row_length, theta_data);
      if(mode == "tomo" && useweight == true){
        extract_from_table(pos_data, 2, row_length, w_data);
      }
    }
  }

  if((data_file_found == true) && (calc_DD == "yes")){
    cout << "Loaded data file." << endl;
  }

  // --------------------------------------------------------------------------------- //

  vector<double>pos_rand;

  if(rand_file_found == true){
    row_length = read_ascii_table(rand_file, pos_rand);
    column_length = pos_rand.size()/row_length;
  }

  vector<double>x_rand;
  vector<double>y_rand;
  vector<double>z_rand;
  vector<double>phi_rand;
  vector<double>theta_rand;
  vector<double>w_rand;

  if(data_file_found == true){
    if(mode == "2d" || mode == "3d" || mode == "poly"){
      extract_from_table(pos_rand, 0, row_length, x_rand);
      extract_from_table(pos_rand, 1, row_length, y_rand);
      if(mode == "2d" && useweight == true){
        extract_from_table(pos_rand, 2, row_length, w_rand);
      }
      if(mode == "3d" || mode == "poly"){
        extract_from_table(pos_rand, 2, row_length, z_rand);
        if(mode == "3d" && useweight == true){
          extract_from_table(pos_rand, 3, row_length, w_rand);
        }
        if(mode == "poly" && useweight == true){
          extract_from_table(pos_rand, 3, row_length, w_rand);
        }
      }
    }
    else if(mode == "tomo"){
      extract_from_table(pos_rand, 0, row_length, phi_rand);
      extract_from_table(pos_rand, 1, row_length, theta_rand);
      if(mode == "tomo" && useweight == true){
        extract_from_table(pos_rand, 2, row_length, w_rand);
      }
    }
  }

  if((rand_file_found == true) && ((calc_DR == "yes") || (calc_RR == "yes"))){
    cout << "Loaded randoms file." << endl;
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating DD", 1);

  vector<double>dd;

  if(mode == "poly"){
    for(int i = 0; i < numbins*mu_bins; i++){
      dd.push_back(0.);
    }
  }
  else{
    for(int i = 0; i < numbins; i++){
      dd.push_back(0.);
    }
  }

  if(calc_DD == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_dd_2d(dd, minimum, maximum, numbins, uselog, x_data, y_data);
      }
      else{
        get_dd_2d_w(dd, minimum, maximum, numbins, uselog, x_data, y_data, w_data);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_dd_3d(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data);
      }
      else{
        get_dd_3d_w(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data, w_data);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_dd_tomo(dd, minimum, maximum, numbins, uselog, phi_data, theta_data);
      }
      else{
        get_dd_tomo_w(dd, minimum, maximum, numbins, uselog, phi_data, theta_data, w_data);
      }
    }
    else if(mode == "poly"){
      if(useweight == false){
        get_dd_poly(dd, minimum, maximum, numbins, mu_bins, uselog, x_data, y_data, z_data);
      }
      else{
        get_dd_poly_w(dd, minimum, maximum, numbins, mu_bins, uselog, x_data, y_data, z_data, w_data);
      }
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating DR", 1);

  vector<double>dr;

  if(mode == "poly"){
    for(int i = 0; i < numbins*mu_bins; i++){
      dr.push_back(0.);
    }
  }
  else{
    for(int i = 0; i < numbins; i++){
      dr.push_back(0.);
    }
  }

  if(calc_DR == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_dr_2d(dr, minimum, maximum, numbins, uselog, x_data, y_data, x_rand, y_rand);
      }
      else{
        get_dr_2d_w(dr, minimum, maximum, numbins, uselog, x_data, y_data, x_rand, y_rand, w_data, w_rand);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_dr_3d(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand);
      }
      else{
        get_dr_3d_w(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand, w_data, w_rand);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_dr_tomo(dr, minimum, maximum, numbins, uselog, phi_data, theta_data, phi_rand, theta_rand);
      }
      else{
        get_dr_tomo_w(dr, minimum, maximum, numbins, uselog, phi_data, theta_data, phi_rand, theta_rand, w_data, w_rand);
      }
    }
    else if(mode == "poly"){
      if(useweight == false){
        get_dr_poly(dr, minimum, maximum, numbins, mu_bins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand);
      }
      else{
        get_dr_poly_w(dr, minimum, maximum, numbins, mu_bins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand, w_data, w_rand);
      }
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Calculating RR", 1);

  vector<double>rr;

  if(mode == "poly"){
    for(int i = 0; i < numbins*mu_bins; i++){
      rr.push_back(0.);
    }
  }
  else{
    for(int i = 0; i < numbins; i++){
      rr.push_back(0.);
    }
  }

  if(calc_RR == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_dd_2d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand);
      }
      else{
        get_dd_2d_w(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, w_rand);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_dd_3d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand);
      }
      else{
        get_dd_3d_w(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand, w_rand);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_dd_tomo(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand);
      }
      else{
        get_dd_tomo_w(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand, w_rand);
      }
    }
    else if(mode == "poly"){
      if(useweight == false){
        get_dd_poly(rr, minimum, maximum, numbins, mu_bins, uselog, x_rand, y_rand, z_rand);
      }
      else{
        get_dd_poly_w(rr, minimum, maximum, numbins, mu_bins, uselog, x_rand, y_rand, z_rand, w_rand);
      }
    }
  }

  // --------------------------------------------------------------------------------- //

  prog.start_process("Save data", 1);

  double _r[numbins], dmu;
  long int r_len;

  if(mode == "poly"){
    r_len = numbins*mu_bins;
  }
  else{
    r_len = numbins;
  }

  double r[r_len], mu[r_len];

  if(mode == "poly"){
    get_r(_r, minimum, maximum, numbins, uselog);
    dmu = 1./((float)mu_bins);
    for(int i = 0; i < mu_bins; i++){
      for(int j = 0; j < numbins; j++){
        r[i*numbins + j] = _r[j];
        mu[i*numbins + j] = 0.5*dmu + dmu*i;
      }
    }
  }
  else{
    get_r(r, minimum, maximum, numbins, uselog);
  }

  Writer wr;
  wr.add2header("TwoFast Two-Point Output");
  if(mode == "poly"){
    wr.store(r, numbins*mu_bins, "r");
    wr.store(mu, numbins*mu_bins, "mu");
  }
  else{
    wr.store(r, numbins, "r");
  }
  if(calc_DD == "yes"){
    if(mode == "poly"){
      wr.store(dd, numbins*mu_bins, "DD");
    }
    else{
      wr.store(dd, numbins, "DD");
    }
  }
  if(calc_DR == "yes"){
    if(mode == "poly"){
      wr.store(dr, numbins*mu_bins, "DR");
    }
    else{
      wr.store(dr, numbins, "DR");
    }
  }
  if(calc_RR == "yes"){
    if(mode == "poly"){
      wr.store(rr, numbins*mu_bins, "RR");
    }
    else{
      wr.store(rr, numbins, "RR");
    }
  }
  /* Should get rid of this or make smarter */
  if(calc_xi == "yes"){
    vector<double>xi;
    if(mode == "poly"){
      for(int i = 0; i < numbins*mu_bins; i++){
        xi.push_back(0.);
      }
      get_xi(x_data.size(), x_rand.size(), dd, dr, rr, xi);
      wr.store(xi, numbins*mu_bins, "xi");
    }
    else{
      for(int i = 0; i < numbins; i++){
        xi.push_back(0.);
      }
      get_xi(x_data.size(), x_rand.size(), dd, dr, rr, xi);
      wr.store(xi, numbins, "xi");
    }
  }
  wr.write2file(out_file);

  prog.end();

  return 0;
}
