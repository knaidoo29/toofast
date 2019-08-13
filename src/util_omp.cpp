#include "../include/twofast.h"
#include <math.h>

using namespace std;

void partition(long int total, int processors, int *myid, long int *begin, long int *end){
  /* Partitions data according to the processor name.

  Parameters
  ----------
  long int total:
    Size of the full data to be partitioned.
  int processors:
    The number of processors being used for the MPI process.
  int myid:
    The current processors ID.
  long int begin:
    The integer for the end of the divided for-loop iterator.
  long int end:
    The integer for the end of the divided for-loop iterator.
  */
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
  /* Gets the indexes of a lower triangle matrix of size n*n for a given count value.

  Parameters
  ----------
  long int count:
    Count value for which we need the indexes of.
  long int i:
    The index for the row.
  long int j:
    The index for the column.
  */
  float row;
  row = (1. + sqrt(1. + 8.*((float) *count)))/2.;
  *i = (int) (row - fmod(row, 1.));
  *j = (int) (*i * fmod((row-1.), 1.));
}
