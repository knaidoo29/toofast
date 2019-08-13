//#ifndef UTIL_OMP_H
//#define UTIL_OMP_H

void partition(long int total, int processors, int *myid, long int *begin, long int *end);

void get_triangle_ij(long int *count, long int *i, long int *j);

//#endif