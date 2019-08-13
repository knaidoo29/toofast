#include "mpi.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <cctype>
#include "toofast.h"

using namespace std;

// FUNCTION DECLARATION -------------------------------------------------------------- //

void extract_param(string paramfile, string whichparam, string *param, bool *found);

void summarise_param(string whichparam, string param, bool found);

void partition(long int total, int processors, int *myid, long int *begin, long int *end);

void get_triangle_ij(long int *count, long int *i, long int *j);

void get_mpi_dd_3d(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dr_3d(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_data, vector<float> &y_data, vector<float> &z_data,
  vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

void to_lowercase(char& c);

void get_xi(long int num_data, long int num_rand, vector<float> &dd,
  vector<float> &dr, vector<float> &rr, vector<float> &xi);

// MAIN ============================================================================== //

int main(int argc, char** argv){

  // MPI START ----------------------------------------------------------------------- //

  MPI_Status status;

  MPI_Init(&argc, &argv); // starts MPI

  int myid, processors;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &processors); // get number of processes

  progress prog;

  if(myid == 0){
    prog.start("Testing Two Point Correlation Function");
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

  if(data_file_found == true){
    extract_from_table(pos_data, 0, row_length, x_data);
    extract_from_table(pos_data, 1, row_length, y_data);
    extract_from_table(pos_data, 2, row_length, z_data);
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

  if(data_file_found == true){
    extract_from_table(pos_rand, 0, row_length, x_rand);
    extract_from_table(pos_rand, 1, row_length, y_rand);
    extract_from_table(pos_rand, 2, row_length, z_rand);
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

  total_size = x_data.size()*(x_data.size()-1)/2;

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
    if(mode == "3d"){
      get_mpi_dd_3d(dd, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
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

  total_size = x_data.size();

  partition(x_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);

  if(calc_DR == "yes"){
    if(mode == "3d"){
      get_mpi_dr_3d(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
        x_rand, y_rand, z_rand, &partition_begin_dr, &partition_end_dr, &prefix);
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

  total_size = x_rand.size()*(x_rand.size()-1)/2;

  partition(total_size, processors, &myid, &partition_begin_rr, &partition_end_rr);

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
    if(mode == "3d"){
      get_mpi_dd_3d(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand,
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

// END MAIN ========================================================================= //

// FUNCTIONS ------------------------------------------------------------------------- //

void extract_param(string paramfile, string whichparam, string *param, bool *found){
  /* Extracts a parameter from a parameter file.

  Parameters
  ----------
  string paramfile:
    Name of the parameter file.
  string whichparam:
    Name of the parameter we want to extract from the file.
  string param
    Sets the param value according to what was found in the parameter file.
  bool found
  */
  string _line;
  int check = 0;
  ifstream _file (paramfile);
  if(_file.is_open()){
    while(getline(_file, _line)){
      if(_line.find(whichparam) == 0){
        _line.erase(_line.find(whichparam), whichparam.length());
        remove(_line.begin(), _line.end(), ' ');
        if(_line.find('%') < _line.size()){
          _line.erase(_line.find('%'), _line.size()-_line.find('%'));
        }
        *param = _line;
        if(_line.size() == 0){
          *found = false;
        }
        else{
          *found = true;
        }
        check = 1;
      }
    }
  }
  if(check == 0){
    cout << "Parameter " << whichparam << " was not found in the parameter file " << paramfile << endl;
  }
}

void summarise_param(string whichparam, string param, bool found){
  if(found == true){
    cout << whichparam << " = " << param << "  (Found = Yes)" << endl;
  }
  else{
    cout << whichparam << " = " << param << "  (Found = No)" << endl;
  }
}

void partition(long int total, int processors, int *myid, long int *begin, long int *end){

  long int partition_size = total/processors, partition_excess = total%processors;

  *begin = partition_size*(*myid);

  if(*myid == processors-1){
    *end = partition_size*(*myid+1)+partition_excess;
  }
  else{
    *end = partition_size*(*myid+1);
  }
}

void get_triangle_ij(long int *count, long int *i, long int *j){
  float row;
  row = (1. + sqrt(1. + 8.*((float) *count)))/2.;
  *i = (int) (row - fmod(row, 1.));
  *j = (int) (*i * fmod((row-1.), 1.));
}

void get_mpi_dd_3d(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z, long int *partition_begin,
  long int *partition_end, long int *partition_begin_i, long int *partition_begin_j,
  long int *partition_end_i, long int *partition_end_j, string *prefix){

    long int k = 0, ind;
    float dx, log10min, log10max, log10dx, dist;

    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    long int j_start, j_end;

    double part_total = (double) (*partition_end - *partition_begin);

    long int _part_total, _temp_begin, _temp_end;

    _temp_begin = (*partition_begin_i+1)*(*partition_begin_i)/2;
    _temp_end = (*partition_end_i+1)*(*partition_end_i)/2;

    _temp_begin -= *partition_begin_i - *partition_begin_j + 1;
    _temp_end -= *partition_end_i - *partition_end_j;

    _part_total = _temp_end - _temp_begin;

    for(long int i = *partition_begin_i; i <= *partition_end_i; i++){
      if((i == *partition_begin_i) && (i == *partition_end_i)){
        j_start = *partition_begin_j;
        j_end = *partition_end_j;
      }
      else if(i == *partition_begin_i){
        j_start = *partition_begin_j;
        j_end = i-1;
      }
      else if(i == *partition_end_i){
        j_start = 0;
        j_end = *partition_end_j;
      }
      else{
        j_start = 0;
        j_end = i-1;
      }
      for(long int j = j_start; j <= j_end; j++){
        dist = get_distance_3d(x[i], y[i], z[i], x[j], y[j], z[j]);
        if (dist >= minimum and dist <= maximum){
          if(uselog == true){
            ind = floor((log10(dist) - log10min)/log10dx);
          }
          else{
            ind = floor((dist - minimum)/dx);
          }
          dd[ind] += 1.;
        }
        k += 1;
        if(k % (_part_total/10) == 0){
          cout << *prefix << "Auto Pairs = " << k/(_part_total/10) << '/' << 10 << endl;
        }
      }
    }
    //cout << *prefix << k << endl;
    for(int i = 0; i < numbins; i++){
      dd[i] *= 2.;
    }
}

void get_mpi_dr_3d(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_data, vector<float> &y_data, vector<float> &z_data,
  vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand,
  long int *partition_begin, long int *partition_end, string *prefix){

    long int ind;
    float dx, log10min, log10max, log10dx, dist;

    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    for(long int i = *partition_begin; i < *partition_end; i++){
      for(long int j = 0; j < x_rand.size(); j++){
        dist = get_distance_3d(x_data[i], y_data[i], z_data[i], x_rand[j], y_rand[j], z_rand[j]);
        if (dist >= minimum and dist <= maximum){
          if(uselog == true){
            ind = floor((log10(dist) - log10min)/log10dx);
          }
          else{
            ind = floor((dist - minimum)/dx);
          }
          dr[ind] += 1.;
        }
      }
      if((i - *partition_begin) % ((*partition_end - *partition_begin)/10) == 0){
        cout << *prefix << "Cross Pairs = " << (i - *partition_begin)/((*partition_end - *partition_begin)/10) + 1 << '/' << 10 << endl;
      }
    }
}

void to_lowercase(char& c){
  c = tolower(static_cast<unsigned char>(c));
}

void get_xi(long int num_data, long int num_rand, vector<float> &dd,
  vector<float> &dr, vector<float> &rr, vector<float> &xi){
  float nd, nr;
  nd = (float) num_data;
  nr = (float) num_rand;
  for(int i = 0; i < dd.size(); i++){
    xi[i] = 1. + pow(nr/nd, 2.)*(dd[i]/rr[i]) - 2.*(nr/nd)*(dr[i]/rr[i]);
  }
}
