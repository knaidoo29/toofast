//#ifndef TWO_POINT_H
//#define TWO_POINT_H

#include <vector>

using namespace std;

float get_distance_3d(float x1, float y1, float z1, float x2, float y2, float z2);

void get_r(float r[], float minimum, float maximum, int numbins, bool uselog);

void get_xi(long int num_data, long int num_rand, vector<float> &dd, vector<float> &dr, 
  vector<float> &rr, vector<float> &xi);

void get_dd_3d(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z);

void get_dr_3d(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_data, vector<float> &y_data, vector<float> &z_data,
  vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand);

void get_mpi_dd_3d(vector<float> &dd, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dr_3d(vector<float> &dr, float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_data, vector<float> &y_data, vector<float> &z_data,
  vector<float> &x_rand, vector<float> &y_rand, vector<float> &z_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

//#endif