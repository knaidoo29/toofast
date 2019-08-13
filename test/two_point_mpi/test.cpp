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

void get_triangle_ij(long int *count, long int *i, long int *j){
  float row;
  row = (1. + sqrt(1. + 8.*((float) *count)))/2.;
  *i = (int) (row - fmod(row, 1.));
  *j = (int) (*i * fmod((row-1.), 1.));
}

int main(int argc, char** argv){

  // MPI START -------------------------------------------------------------- //

  MPI_Status status;

  MPI_Init(&argc, &argv); // starts MPI

  int myid, processors;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); // get current process id
  MPI_Comm_size(MPI_COMM_WORLD, &processors); // get number of processes

  // MAIN CODE -------------------------------------------------------------- //

  long int num_data = 10;
  long int total = num_data*(num_data - 1)/2;

  cout << total << endl;

  long int total_excess = total%processors;
  long int total_div = (total - total_excess)/processors;

  long int part_begin, part_end;

  part_begin = myid*total_div;
  part_end = (myid+1)*total_div - 1;
  if(myid + 1 == processors){
    part_end += total_excess;
  }

  cout << part_begin << '\t' << part_end << endl;

  long int part_begin_i, part_begin_j, part_end_i, part_end_j;

  get_triangle_ij(&part_begin, &part_begin_i, &part_begin_j);

  cout << myid << '\t' << part_begin_i << '\t' << part_begin_j << endl;

  get_triangle_ij(&part_end, &part_end_i, &part_end_j);

  cout << myid << '\t' << part_end_i << '\t' << part_end_j << endl;

  // MPI END ---------------------------------------------------------------- //

  MPI_Finalize(); // Finalize the MPI environment

  return 0;
}
