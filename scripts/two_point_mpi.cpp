#include "mpi.h"
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

  // MPI START ----------------------------------------------------------------------- //

  MPI_Status status;

  MPI_Init(&argc, &argv); // starts MPI

  int myid, processors;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &processors); // get number of processes

  progress prog;

  if(myid == 0){
    prog.start("TwoFast: 2-Point");
  }

  // READING PARAMFILE --------------------------------------------------------------- //

  if(myid == 0){
    prog.start_process("Reading Parameter File", 1);
  }

  string prefix, paramfile, whichparam, data_file, rand_file, out_file;
  bool data_file_found = false, rand_file_found = false, out_file_found = false;

  string mode, calc_DD, calc_DR, calc_RR, calc_xi, uselog_str, useweight_str, minimum_str, maximum_str, numbins_str;
  bool mode_found = false, calc_DD_found = false, calc_DR_found = false, calc_RR_found = false, calc_xi_found = false;
  bool uselog_found = false, useweight_found = false, minimum_found = false, maximum_found = false, numbins_found = false;

  prefix = "processor = " + to_string(myid+1) + "/" + to_string(processors) + ":  ";
  paramfile = argv[1];

  if(myid == 0){
    cout << "Paramfile = " << paramfile << endl;
  }

  whichparam = "data_file";
  extract_param(paramfile, whichparam, &data_file, &data_file_found);
  if(myid == 0){
    summarise_param(whichparam, data_file, data_file_found);
  }

  whichparam = "rand_file";
  extract_param(paramfile, whichparam, &rand_file, &rand_file_found);
  if(myid == 0){
    summarise_param(whichparam, rand_file, rand_file_found);
  }

  whichparam = "out_file";
  extract_param(paramfile, whichparam, &out_file, &out_file_found);
  if(myid == 0){
    summarise_param(whichparam, out_file, out_file_found);
  }

  whichparam = "mode";
  extract_param(paramfile, whichparam, &mode, &mode_found);
  for_each(mode.begin(), mode.end(), to_lowercase);
  if(myid == 0){
  summarise_param(whichparam, mode, mode_found);
  }

  whichparam = "calc_DD";
  extract_param(paramfile, whichparam, &calc_DD, &calc_DD_found);
  for_each(calc_DD.begin(), calc_DD.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, calc_DD, calc_DD_found);
  }

  whichparam = "calc_DR";
  extract_param(paramfile, whichparam, &calc_DR, &calc_DR_found);
  for_each(calc_DR.begin(), calc_DR.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, calc_DR, calc_DR_found);
  }

  whichparam = "calc_RR";
  extract_param(paramfile, whichparam, &calc_RR, &calc_RR_found);
  for_each(calc_RR.begin(), calc_RR.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, calc_RR, calc_RR_found);
  }

  whichparam = "calc_xi";
  extract_param(paramfile, whichparam, &calc_xi, &calc_xi_found);
  for_each(calc_xi.begin(), calc_xi.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, calc_xi, calc_xi_found);
  }

  whichparam = "uselog";
  extract_param(paramfile, whichparam, &uselog_str, &uselog_found);
  for_each(uselog_str.begin(), uselog_str.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, uselog_str, uselog_found);
  }

  whichparam = "useweight";
  extract_param(paramfile, whichparam, &useweight_str, &useweight_found);
  for_each(useweight_str.begin(), useweight_str.end(), to_lowercase);
  if(myid == 0){
    summarise_param(whichparam, useweight_str, useweight_found);
  }

  whichparam = "minimum";
  extract_param(paramfile, whichparam, &minimum_str, &minimum_found);
  if(myid == 0){
    summarise_param(whichparam, minimum_str, minimum_found);
  }

  whichparam = "maximum";
  extract_param(paramfile, whichparam, &maximum_str, &maximum_found);
  if(myid == 0){
    summarise_param(whichparam, maximum_str, maximum_found);
  }

  whichparam = "numbins";
  extract_param(paramfile, whichparam, &numbins_str, &numbins_found);
  if(myid == 0){
    summarise_param(whichparam, numbins_str, numbins_found);
  }

  int numbins;
  double minimum, maximum;
  bool uselog, useweight;

  if(uselog_str == "yes"){
    uselog = true;
  }
  else if(uselog_str == "no"){
    uselog = false;
  }
  else{
    if(myid == 0){
      cout << "Error: uselog is ill defined -- found in paramfile = " << uselog_str << endl;
      cout << "Should be either yes/no" << endl;
      cout << "Will assume 'no'" << endl;
      uselog = false;
    }
  }

  if(useweight_str == "yes"){
    useweight = true;
  }
  else if(useweight_str == "no"){
    useweight = false;
  }
  else{
    if(myid == 0){
      cout << "Error: useweight is ill defined -- found in paramfile = " << useweight_str << endl;
      cout << "Should be either yes/no" << endl;
      cout << "Will assume 'no'" << endl;
      useweight = false;
    }
  }

  minimum = stof(minimum_str);
  maximum = stof(maximum_str);
  numbins = stof(numbins_str);

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Reading data", 1);
  }

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
    if(mode == "2d" || mode == "3d"){
      extract_from_table(pos_data, 0, row_length, x_data);
      extract_from_table(pos_data, 1, row_length, y_data);
      if(mode == "2d" && useweight == true){
        extract_from_table(pos_data, 2, row_length, w_data);
      }
      if(mode == "3d"){
        extract_from_table(pos_data, 2, row_length, z_data);
        if(mode == "3d" && useweight == true){
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
    if(myid == 0){
      cout << "Loaded data file." << endl;
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

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
    if(mode == "2d" || mode == "3d"){
      extract_from_table(pos_rand, 0, row_length, x_rand);
      extract_from_table(pos_rand, 1, row_length, y_rand);
      if(mode == "2d" && useweight == true){
        extract_from_table(pos_rand, 2, row_length, w_rand);
      }
      if(mode == "3d"){
        extract_from_table(pos_rand, 2, row_length, z_rand);
        if(mode == "3d" && useweight == true){
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
    if(myid == 0){
      cout << "Loaded randoms file." << endl;
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating DD", 1);
  }

  vector<double>dd;

  for(int i = 0; i < numbins; i++){
    dd.push_back(0.);
  }

  long int total_size_dd, partition_begin_dd, partition_end_dd;

  if(mode == "2d" || mode == "3d"){
    if(x_data.size() % 2 == 0){
      total_size_dd = (x_data.size()/2)*(x_data.size()-1);
    }
    else{
      total_size_dd = x_data.size()*((x_data.size()-1)/2);
    }
  }
  else if(mode == "tomo"){
    if(phi_data.size() % 2 == 0){
      total_size_dd = (phi_data.size()/2)*(phi_data.size()-1);
    }
    else{
      total_size_dd = phi_data.size()*((phi_data.size()-1)/2);
    }
  }

  partition(total_size_dd, processors, &myid, &partition_begin_dd, &partition_end_dd);

  long int part_begin_dd_i, part_begin_dd_j, part_end_dd_i, part_end_dd_j;

  get_triangle_ij(&partition_begin_dd, &part_begin_dd_i, &part_begin_dd_j);
  get_triangle_ij(&partition_end_dd, &part_end_dd_i, &part_end_dd_j);

  part_end_dd_j -= 1;

  if(part_end_dd_j < 0){
    part_end_dd_i -= 1;
    part_end_dd_j = part_end_dd_i - 1;
  }

  if(calc_DD == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_mpi_dd_2d(dd, minimum, maximum, numbins, uselog, x_data, y_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
      else{
        get_mpi_dd_2d_w(dd, minimum, maximum, numbins, uselog, x_data, y_data, w_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_mpi_dd_3d(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
      else{
        get_mpi_dd_3d_w(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data, w_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_mpi_dd_tomo(dd, minimum, maximum, numbins, uselog, phi_data, theta_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
      else{
        get_mpi_dd_tomo_w(dd, minimum, maximum, numbins, uselog, phi_data, theta_data, w_data,
          &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
          &part_end_dd_i, &part_end_dd_j, &prefix);
      }
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating DR", 1);
  }

  vector<double>dr;

  for(int i = 0; i < numbins; i++){
    dr.push_back(0.);
  }

  long int partition_begin_dr, partition_end_dr;

  if(mode == "2d" || mode == "3d"){
    partition(x_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);
  }
  else if(mode == "tomo"){
    partition(phi_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);
  }

  if(calc_DR == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_mpi_dr_2d(dr, minimum, maximum, numbins, uselog, x_data, y_data,
          x_rand, y_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
      else{
        get_mpi_dr_2d_w(dr, minimum, maximum, numbins, uselog, x_data, y_data,
          x_rand, y_rand, w_data, w_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_mpi_dr_3d(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
      else{
        get_mpi_dr_3d_w(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
          x_rand, y_rand, z_rand, w_data, w_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_mpi_dr_tomo(dr, minimum, maximum, numbins, uselog, phi_data, theta_data,
          phi_rand, theta_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
      else{
        get_mpi_dr_tomo_w(dr, minimum, maximum, numbins, uselog, phi_data, theta_data,
          phi_rand, theta_rand, w_data, w_rand, &partition_begin_dr, &partition_end_dr, &prefix);
      }
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating RR", 1);
  }

  vector<double>rr;

  for(int i = 0; i < numbins; i++){
    rr.push_back(0.);
  }

  long int total_size_rr, partition_begin_rr, partition_end_rr;

  if(mode == "2d" || mode == "3d"){
    if(x_rand.size() % 2 == 0){
      total_size_rr = (x_rand.size()/2)*(x_rand.size()-1);
    }
    else{
      total_size_rr = x_rand.size()*((x_rand.size()-1)/2);
    }
  }
  else if(mode == "tomo"){
    if(phi_data.size() % 2 == 0){
      total_size_rr = (phi_rand.size()/2)*(phi_rand.size()-1);
    }
    else{
      total_size_rr = phi_rand.size()*((phi_rand.size()-1)/2);
    }
  }

  partition(total_size_rr, processors, &myid, &partition_begin_rr, &partition_end_rr);

  long int part_begin_rr_i, part_begin_rr_j, part_end_rr_i, part_end_rr_j;

  get_triangle_ij(&partition_begin_rr, &part_begin_rr_i, &part_begin_rr_j);
  get_triangle_ij(&partition_end_rr, &part_end_rr_i, &part_end_rr_j);

  part_end_rr_j -= 1;

  if(part_end_rr_j < 0){
    part_end_rr_i -= 1;
    part_end_rr_j = part_end_rr_i - 1;
  }

  if(calc_RR == "yes"){
    if(mode == "2d"){
      if(useweight == false){
        get_mpi_dd_2d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
      else{
        get_mpi_dd_2d_w(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, w_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
    }
    else if(mode == "3d"){
      if(useweight == false){
        get_mpi_dd_3d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
      else{
        get_mpi_dd_3d_w(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand, w_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
    }
    else if(mode == "tomo"){
      if(useweight == false){
        get_mpi_dd_tomo(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
      else{
        get_mpi_dd_tomo_w(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand, w_rand,
          &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
          &part_end_rr_i, &part_end_rr_j, &prefix);
      }
    }
  }
  
  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Save data", 1);
  }

  double r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);

  int dest = 0, tag1 = 0, tag2 = 1, tag3 = 2;
  if(myid == 0){
    vector<double>dd_total;
    vector<double>dr_total;
    vector<double>rr_total;
    for(long int i = 0; i < numbins; i++){
      dd_total.push_back(dd[i]);
      dr_total.push_back(dr[i]);
      rr_total.push_back(rr[i]);
    }
    for(int source = 1; source < processors; source++){
      MPI_Recv(&dd[0], dd.size(), MPI_DOUBLE, source, tag1, MPI_COMM_WORLD, &status);
      for(long int i = 0; i < dd.size(); i++){
        dd_total[i] += dd[i];
      }
      MPI_Recv(&dr[0], dr.size(), MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
      for(long int i = 0; i < dr.size(); i++){
        dr_total[i] += dr[i];
      }
      MPI_Recv(&rr[0], rr.size(), MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
      for(long int i = 0; i < rr.size(); i++){
        rr_total[i] += rr[i];
      }
    }
    Writer wr;
    wr.add2header("Levy flight Xi estimation");
    wr.store(r, numbins, "r");
    if(calc_DD == "yes"){
      wr.store(dd_total, numbins, "DD");
    }
    if(calc_DR == "yes"){
      wr.store(dr_total, numbins, "DR");
    }
    if(calc_RR == "yes"){
      wr.store(rr_total, numbins, "RR");
    }
    if(calc_xi == "yes"){
      vector<double>xi;
      for(long int i = 0; i < numbins; i++){
        xi.push_back(0.);
      }
      get_xi(x_data.size(), x_rand.size(), dd_total, dr_total, rr_total, xi);
      wr.store(xi, numbins, "xi");
    }
    wr.write2file(out_file);
  }
  else{
    MPI_Send(&dd[0], dd.size(), MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD);
    MPI_Send(&dr[0], dr.size(), MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
    MPI_Send(&rr[0], rr.size(), MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD);
  }

  if(myid == 0){
    prog.end();
  }

  // MPI END ------------------------------------------------------------------------- //

  MPI_Finalize(); // Finalize the MPI environment

  return 0;
}
