//#ifndef TWO_POINT_H
//#define TWO_POINT_H

#include <vector>
#include <string>

using namespace std;

double get_distance_2d(double x1, double y1, double x2, double y2);

double get_distance_3d(double x1, double y1, double z1, double x2, double y2, double z2);

double get_distance_tomo(double phi1, double theta1, double phi2, double theta2);

void get_r(double r[], double minimum, double maximum, int numbins, bool uselog);

void get_xi(long int num_data, long int num_rand, vector<double> &dd, vector<double> &dr,
  vector<double> &rr, vector<double> &xi);

void get_dd_2d(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y);

void get_dd_2d_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &w);

void get_dr_2d(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &x_rand, vector<double> &y_rand);

void get_dr_2d_w(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &x_rand, vector<double> &y_rand,
  vector<double> &w_data, vector<double> &w_rand);

void get_dd_3d(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z);

void get_dd_3d_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &w);

void get_dr_3d(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand);

void get_dr_3d_w(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand,
  vector<double> &w_data, vector<double> &w_rand);

void get_dd_tomo(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi, vector<double> &theta);

void get_dd_tomo_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi, vector<double> &theta, vector<double> &w);

void get_dr_tomo(vector<double> &dr,double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi_data, vector<double> &theta_data, vector<double> &phi_rand, vector<double> &theta_rand);

void get_dr_tomo_w(vector<double> &dr,double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi_data, vector<double> &theta_data, vector<double> &phi_rand, vector<double> &theta_rand,
  vector<double> &w_data, vector<double> &w_rand);

void get_dd_poly(vector<double> &dd, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z);

void get_dd_poly_w(vector<double> &dd, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &w);

void get_dr_poly(vector<double> &dr, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand);

void get_dr_poly_w(vector<double> &dr, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand,
  vector<double> &w_data, vector<double> &w_rand);

void get_mpi_dd_2d(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dd_2d_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &w, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dr_2d(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &x_rand, vector<double> &y_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dr_2d_w(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &x_rand, vector<double> &y_rand,
  vector<double> &w_data, vector<double> &w_rand, long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dd_3d(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dd_3d_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &w, long int *partition_begin,
  long int *partition_end, long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dr_3d(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dr_3d_w(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand,
  vector<double> &w_data, vector<double> &w_rand, long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dd_tomo(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi, vector<double> &theta, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dd_tomo_w(vector<double> &dd, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi, vector<double> &theta, vector<double> &w, long int *partition_begin, long int *partition_end,
  long int *partition_begin_i, long int *partition_begin_j, long int *partition_end_i,
  long int *partition_end_j, string *prefix);

void get_mpi_dr_tomo(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi_data, vector<double> &theta_data, vector<double> &phi_rand, vector<double> &theta_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dr_tomo_w(vector<double> &dr, double minimum, double maximum, int numbins, bool uselog,
  vector<double> &phi_data, vector<double> &theta_data, vector<double> &phi_rand, vector<double> &theta_rand,
  vector<double> &w_data, vector<double> &w_rand, long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dd_poly(vector<double> &dd, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, long int *partition_begin,
  long int *partition_end, long int *partition_begin_i, long int *partition_begin_j,
  long int *partition_end_i, long int *partition_end_j, string *prefix);

void get_mpi_dd_poly_w(vector<double> &dd, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &w, long int *partition_begin,
  long int *partition_end, long int *partition_begin_i, long int *partition_begin_j,
  long int *partition_end_i, long int *partition_end_j, string *prefix);

void get_mpi_dr_poly(vector<double> &dr, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

void get_mpi_dr_poly_w(vector<double> &dr, double minimum, double maximum, int numbins, int mu_bins, bool uselog,
  vector<double> &x_data, vector<double> &y_data, vector<double> &z_data,
  vector<double> &x_rand, vector<double> &y_rand, vector<double> &z_rand, vector<double> &w_data, vector<double> &w_rand,
  long int *partition_begin, long int *partition_end, string *prefix);

//#endif
