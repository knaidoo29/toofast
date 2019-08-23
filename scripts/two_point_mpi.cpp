#include "mpi.h"
#include <iostream>
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

  string mode, calc_DD, calc_DR, calc_RR, calc_xi, uselog_str, minimum_str, maximum_str, numbins_str;
  bool mode_found = false, calc_DD_found = false, calc_DR_found = false, calc_RR_found = false, calc_xi_found = false;
  bool uselog_found = false, minimum_found = false, maximum_found = false, numbins_found = false;

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
  float minimum, maximum;
  bool uselog;

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

  minimum = stof(minimum_str);
  maximum = stof(maximum_str);
  numbins = stof(numbins_str);

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Reading data", 1);
  }

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
    if(myid == 0){
      cout << "Loaded data file." << endl;
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

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
    if(myid == 0){
      cout << "Loaded randoms file." << endl;
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating DD", 1);
  }

  vector<float>dd;

  for(int i = 0; i < numbins; i++){
    dd.push_back(0.);
  }

  long int total_size, partition_begin_dd, partition_end_dd;
  
  if(mode == "2d" || mode == "3d"){
    total_size = x_data.size()*(x_data.size()-1)/2;
  }
  else if(mode == "tomo"){
    total_size = phi_data.size()*(phi_data.size()-1)/2;
  }

  partition(total_size, processors, &myid, &partition_begin_dd, &partition_end_dd);

  long int part_begin_dd_i, part_begin_dd_j, part_end_dd_i, part_end_dd_j;

  get_triangle_ij(&partition_begin_dd, &part_begin_dd_i, &part_begin_dd_j);
  get_triangle_ij(&partition_end_dd, &part_end_dd_i, &part_end_dd_j);

  if(part_end_dd_j == 0){
    part_end_dd_i = part_end_dd_i - 1;
    part_end_dd_j = part_end_dd_i - 1;
  }
  else{
    part_end_dd_j = part_end_dd_j - 1;
  }

  partition_end_dd = partition_end_dd - 1;

  if(calc_DD == "yes"){
    if(mode == "2d"){
      get_mpi_dd_2d(dd, minimum, maximum, numbins, uselog, x_data, y_data,
        &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
        &part_end_dd_i, &part_end_dd_j, &prefix);
    }
    else if(mode == "3d"){
      get_mpi_dd_3d(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
        &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
        &part_end_dd_i, &part_end_dd_j, &prefix);
    }
    else if(mode == "tomo"){
      get_mpi_dd_tomo(dd, minimum, maximum, numbins, uselog, phi_data, theta_data,
        &partition_begin_dd, &partition_end_dd, &part_begin_dd_i, &part_begin_dd_j,
        &part_end_dd_i, &part_end_dd_j, &prefix);
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating DR", 1);
  }

  vector<float>dr;

  for(int i = 0; i < numbins; i++){
    dr.push_back(0.);
  }

  long int partition_begin_dr, partition_end_dr;

  if(mode == "2d" || mode == "3d"){
    total_size = x_data.size();
    partition(x_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);
  }
  else if(mode == "tomo"){
    total_size = phi_data.size();
    partition(phi_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);
  }

  if(calc_DR == "yes"){
    if(mode == "2d"){
      get_mpi_dr_2d(dr, minimum, maximum, numbins, uselog, x_data, y_data,
        x_rand, y_rand, &partition_begin_dr, &partition_end_dr, &prefix);
    }
    else if(mode == "3d"){
      get_mpi_dr_3d(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
        x_rand, y_rand, z_rand, &partition_begin_dr, &partition_end_dr, &prefix);
    }
    else if(mode == "tomo"){
      get_mpi_dr_tomo(dr, minimum, maximum, numbins, uselog, phi_data, theta_data,
        phi_rand, theta_rand, &partition_begin_dr, &partition_end_dr, &prefix);
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Calculating RR", 1);
  }

  vector<float>rr;

  for(int i = 0; i < numbins; i++){
    rr.push_back(0.);
  }

  long int partition_begin_rr, partition_end_rr;

  if(mode == "2d" || mode == "3d"){
    total_size = x_rand.size()*(x_rand.size()-1)/2;
    partition(total_size, processors, &myid, &partition_begin_rr, &partition_end_rr);
  }
  else if(mode == "tomo"){
    total_size = phi_rand.size()*(phi_rand.size()-1)/2;
    partition(total_size, processors, &myid, &partition_begin_rr, &partition_end_rr);
  }

  long int part_begin_rr_i, part_begin_rr_j, part_end_rr_i, part_end_rr_j;

  get_triangle_ij(&partition_begin_rr, &part_begin_rr_i, &part_begin_rr_j);
  get_triangle_ij(&partition_end_rr, &part_end_rr_i, &part_end_rr_j);

  if(part_end_rr_j == 0){
    part_end_rr_i -= 1;
    part_end_rr_j = part_end_rr_i - 1;
  }
  else{
    part_end_rr_j -= 1;
  }

  if(calc_RR == "yes"){
    if(mode == "2d"){
      get_mpi_dd_2d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand,
        &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
        &part_end_rr_i, &part_end_rr_j, &prefix);
    }
    else if(mode == "3d"){
      get_mpi_dd_3d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand,
        &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
        &part_end_rr_i, &part_end_rr_j, &prefix);
    }
    else if(mode == "tomo"){
      get_mpi_dd_tomo(rr, minimum, maximum, numbins, uselog, phi_rand, theta_rand,
        &partition_begin_rr, &partition_end_rr, &part_begin_rr_i, &part_begin_rr_j,
        &part_end_rr_i, &part_end_rr_j, &prefix);
    }
  }

  // --------------------------------------------------------------------------------- //

  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    prog.start_process("Save data", 1);
  }

  float r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);

  int dest = 0, tag = 0;
  if(myid == 0){
    vector<float>dd_total;
    vector<float>dr_total;
    vector<float>rr_total;
    for(int i = 0; i < numbins; i++){
      dd_total.push_back(dd[i]);
      dr_total.push_back(dr[i]);
      rr_total.push_back(rr[i]);
    }
    for(int source = 1; source < processors; source++){
      MPI_Recv(&dd[0], dd.size(), MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
      for(int i = 0; i < dd.size(); i++){
        dd_total[i] += dd[i];
      }
      MPI_Recv(&dr[0], dr.size(), MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
      for(int i = 0; i < dr.size(); i++){
        dr_total[i] += dr[i];
      }
      MPI_Recv(&rr[0], rr.size(), MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
      for(int i = 0; i < rr.size(); i++){
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
      vector<float>xi;
      for(int i = 0; i < numbins; i++){
        xi.push_back(0.);
      }
      get_xi(x_data.size(), x_rand.size(), dd_total, dr_total, rr_total, xi);
      wr.store(xi, numbins, "xi");
    }
    wr.write2file(out_file);
  }
  else{
    MPI_Send(&dd[0], dd.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    MPI_Send(&dr[0], dr.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    MPI_Send(&rr[0], rr.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  }

  if(myid == 0){
    prog.end();
  }

  // MPI END ------------------------------------------------------------------------- //

  MPI_Finalize(); // Finalize the MPI environment

  return 0;
}
