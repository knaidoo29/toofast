#include "mpi.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "toofast.h"

using namespace std;

void partition(long int total, int processors, int *myid, long int *begin, long int *end);

void get_mpi_dd(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  long int *partition_begin, long int *partition_end, vector<float> &x, vector<float> &y, vector<float> &z);

void get_mpi_dr(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  long int *partition_begin, long int *partition_end, vector<float> &x_data, vector<float> &y_data,
  vector<float> &z_data, vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand);

int main(int argc, char** argv){

  // MPI START -------------------------------------------------------------- //

  MPI_Status status;

  MPI_Init(&argc, &argv); /* starts MPI */

  int myid, processors;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* get current process id */
  MPI_Comm_size(MPI_COMM_WORLD, &processors); /* get number of processes */

  // ------------------------------------------------------------------------ //

  string filename_data, filename_rand, filename_output;
  int numbins;
  float minimum, maximum;
  bool uselog;

  filename_data = "data.txt";
  filename_rand = "randoms.txt";
  filename_output = "xi_" + to_string(myid) + ".txt";

  minimum = 0.01;
  maximum = 75.;
  numbins = 40;

  uselog = true;

  // ------------------------------------------------------------------------ //

  cout << "Processor " << myid << " : Reading data " << endl;

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

  cout << "Processor " << myid << " : Reading randoms " << endl;

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

  cout << "Processor " << myid << " : Getting DD " << endl;

  vector<float>dd;

  for(int i = 0; i < numbins; i++){
    dd.push_back(0.);
  }

  long int total_size, partition_begin_dd, partition_end_dd;

  total_size = x_data.size()*(x_data.size()-1)/2;

  partition(total_size, processors, &myid, &partition_begin_dd, &partition_end_dd);

  get_mpi_dd(dd, minimum, maximum, numbins, uselog, &partition_begin_dd, &partition_end_dd, x_data, y_data, z_data);

  // ------------------------------------------------------------------------ //

  cout << "Processor " << myid << " : Getting DR " << endl;

  vector<float>dr;

  for(int i = 0; i < numbins; i++){
    dr.push_back(0.);
  }

  long int partition_begin_dr, partition_end_dr;

  total_size = x_data.size();

  partition(x_data.size(), processors, &myid, &partition_begin_dr, &partition_end_dr);

  get_mpi_dr(dr, minimum, maximum, numbins, uselog, &partition_begin_dr, &partition_end_dr ,
    x_data, y_data, z_data, x_rand, y_rand, z_rand);

  // ------------------------------------------------------------------------ //

  cout << "Processor " << myid << " : Getting RR " << endl;

  vector<float>rr;

  for(int i = 0; i < numbins; i++){
    rr.push_back(0.);
  }

  long int partition_begin_rr, partition_end_rr;

  total_size = x_rand.size()*(x_rand.size()-1)/2;

  partition(total_size, processors, &myid, &partition_begin_rr, &partition_end_rr);

  cout << "rr start = " << partition_begin_rr << " end = " << partition_end_rr << endl;

  get_mpi_dd(rr, minimum, maximum, numbins, uselog, &partition_begin_rr, &partition_end_rr, x_rand, y_rand, z_rand);

  // ------------------------------------------------------------------------ //

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

  //prog.start_process("Save data", 2);

  cout << "Processor " << myid << " : Saving data " << endl;

  Writer wr;
  wr.add2header("Levy flight Xi estimation");
  wr.store(r, numbins, "r");
  wr.store(dd, numbins, "DD");
  wr.store(rr, numbins, "RR");
  wr.store(dr, numbins, "DR");
  //wr.store(xi, numbins, "xi");
  wr.write2file(filename_output);
  wr.clean();

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
    filename_output = "xi_total.txt";
    //Writer wr;
    wr.add2header("Levy flight Xi estimation");
    wr.store(r, numbins, "r");
    wr.store(dd_total, numbins, "DD");
    wr.store(rr_total, numbins, "RR");
    wr.store(dr_total, numbins, "DR");
    //wr.store(xi, numbins, "xi");
    wr.write2file(filename_output);
  }
  else{
    MPI_Send(&dd[0], dd.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    MPI_Send(&dr[0], dr.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    MPI_Send(&rr[0], rr.size(), MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  }

  //prog.end();
  cout << "Processor " << myid << " : Finished " << endl;

  // MPI END ---------------------------------------------------------------- //

  MPI_Finalize(); /* Finalize the MPI environment */

  return 0;
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

void get_mpi_dd(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  long int *partition_begin, long int *partition_end, vector<float> &x, vector<float> &y, vector<float> &z){

    long int k = 0, counter = 0, ind;
    float dx, log10min, log10max, log10dx, dist;

    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    for(long int i = 0; i < x.size(); i++){
      for(long int j = 0; j < i; j++){
        if(k >= *partition_begin){
          if(k < *partition_end){
            counter += 1;
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
          }
        }
        k += 1;
      }
    }
    for(int i = 0; i < numbins; i++){
      dd[i] *= 2.;
    }

}

void get_mpi_dr(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  long int *partition_begin, long int *partition_end, vector<float> &x_data, vector<float> &y_data,
  vector<float> &z_data, vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand){

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
    }
}
