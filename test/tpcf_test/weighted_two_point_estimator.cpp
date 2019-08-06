#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "toofast.h"

using namespace std;

void get_grid(vector<float> &x_grid, vector<float> &y_grid, vector<float> &z_grid,
  float lmin, float lmax, int ngrid){

    float dgrid;
    dgrid = (lmax - lmin)/((float)(ngrid));

    for(int i = 0; i < ngrid; i++){
      for(int j = 0; j < ngrid; j++){
        for(int k = 0; k < ngrid; k++){
          x_grid.push_back(lmin + (dgrid/2.) + ((float)(i))*dgrid);
          y_grid.push_back(lmin + (dgrid/2.) + ((float)(j))*dgrid);
          z_grid.push_back(lmin + (dgrid/2.) + ((float)(k))*dgrid);
        }
      }
    }
}

void get_weight(vector<float> &x, vector<float> &y, vector<float> &z,
  vector<float> &x_grid, vector<float> &y_grid, vector<float> &z_grid,
  vector<float> &weights, float lmin, float lmax, int ngrid){

    float dgrid;
    long int i_ind, j_ind, k_ind;

    dgrid = (lmax - lmin)/((float)(ngrid));

    for(long int i = 0; i < ngrid*ngrid*ngrid; i++){
      weights.push_back(0.);
    }

    for(long int i = 0; i < x.size(); i++){
      i_ind = floor((x[i]-lmin)/dgrid);
      j_ind = floor((y[i]-lmin)/dgrid);
      k_ind = floor((z[i]-lmin)/dgrid);
      weights[k_ind + ngrid*(j_ind + ngrid*i_ind)] += 1.;
    }
}


void get_weighted_dd(float dd[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_grid, vector<float> &y_grid, vector<float> &z_grid, vector<float> &weights){

    long int size, ind;
    long int k = 0, kmax;
    float dx, log10min, log10max, log10dx, dist;

    size = x_grid.size();
    kmax = (size*(size - 1))/2;
    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    for(int i = 0; i < numbins; i++){
      dd[i] = 0.;
    }

    for(long int i = 1; i < size; i++){
      for(long int j = i; j < size; j++){
        dist = get_distance_3d(x_grid[i], y_grid[i], z_grid[i], x_grid[j], y_grid[j], z_grid[j]);
        if (dist >= minimum and dist <= maximum){
          if(uselog == true){
            ind = floor((log10(dist) - log10min)/log10dx);
          }
          else{
            ind = floor((dist - minimum)/dx);
          }
          dd[ind] += weights[i]*weights[j];
        }
        progress_bar(k, kmax);
        k += 1;
      }
    }

    for(int i = 0; i < numbins; i++){
      dd[i] *= 2.;
    }
}


void get_weighted_dr(float dr[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x_grid, vector<float> &y_grid, vector<float> &z_grid,
  vector<float> &weights1, vector<float> &weights2){

    long int size1, size2, ind;
    long int k = 0, kmax;
    float dx, log10min, log10max, log10dx, dist;

    size1 = x_grid.size();
    size2 = size1;
    kmax = size1*size2;
    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    for(int i = 0; i < numbins; i++){
      dr[i] = 0.;
    }

    for(long int i = 0; i < size1; i++){
      for(long int j = 0; j < size2; j++){
        dist = get_distance_3d(x_grid[i], y_grid[i], z_grid[i], x_grid[j], y_grid[j], z_grid[j]);
        if (dist >= minimum and dist <= maximum){
          if(uselog == true){
            ind = floor((log10(dist) - log10min)/log10dx);
          }
          else{
            ind = floor((dist - minimum)/dx);
          }
          dr[ind] += weights1[i]*weights2[j];
        }
        progress_bar(k, kmax);
        k += 1;
      }
    }
}

int main(){

  // ------------------------------------------------------------------------ //

  string filename_data, filename_rand, filename_output;
  int numbins, ngrid;
  float minimum, maximum, lmin, lmax;
  bool uselog;

  // Paramfile -- Proxy

  filename_data = "data.txt";
  filename_rand = "randoms.txt";
  filename_output = "xi_weighted.txt";

  minimum = 75./40.;
  maximum = 75.;
  numbins = 40;

  uselog = true;

  ngrid = 40;
  lmin = 0.;
  lmax = 75.;

  // ------------------------------------------------------------------------ //

  string title = "Testing Two Point Correlation Function";
  progress prog;
  prog.start(title);

  // ------------------------------------------------------------------------ //

  prog.start_process("Read Data", 2);

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

  prog.start_process("Read Randoms", 2);

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

  prog.start_process("Create Grid", 2);

  vector<float>x_grid;
  vector<float>y_grid;
  vector<float>z_grid;

  get_grid(x_grid, y_grid, z_grid, lmin, lmax, ngrid);

  // ------------------------------------------------------------------------ //

  prog.start_process("Get Weights", 2);

  vector<float>weights_d;
  vector<float>weights_r;

  get_weight(x_data, y_data, z_data, x_grid, y_grid, z_grid,
    weights_d, lmin, lmax, ngrid);

  get_weight(x_rand, y_rand, z_rand, x_grid, y_grid, z_grid,
    weights_r, lmin, lmax, ngrid);

  // ------------------------------------------------------------------------ //

  prog.start_process("Calculate Xi(r)", 2);

  prog.start_process("Get DD", 3);

  float dd[numbins];

  get_weighted_dd(dd, minimum, maximum, numbins, uselog, x_grid, y_grid, z_grid, weights_d);

  // ------------------------------------------------------------------------ //

  prog.start_process("Get RR", 3);

  float rr[numbins];

  get_weighted_dd(rr, minimum, maximum, numbins, uselog, x_grid, y_grid, z_grid, weights_r);

  // ------------------------------------------------------------------------ //

  prog.start_process("Get DR", 3);

  float dr[numbins];

  get_weighted_dr(dr, minimum, maximum, numbins, uselog, x_grid, y_grid, z_grid, weights_d, weights_d);

  // ------------------------------------------------------------------------ //

  prog.start_process("Print xi(r)", 3);

  float r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);

  float xi[numbins];

  get_xi(dd, dr, rr, numbins, x_data.size(), x_rand.size(), xi);

  cout << '\t' << '\t' << "r" << '\t' << "Xi(r)" << endl;

  for(int i = 0; i < numbins; i++){
    cout << '\t' << '\t' <<  r[i] << '\t' << xi[i] << endl;
  }

  // ------------------------------------------------------------------------ //

  prog.start_process("Save data", 2);

  Writer wr;
  wr.add2header("Levy flight weigthed Xi estimation");
  wr.store(r, numbins, "r");
  wr.store(dd, numbins, "DD");
  wr.store(rr, numbins, "RR");
  wr.store(dr, numbins, "DR");
  wr.store(xi, numbins, "xi");
  wr.write2file(filename_output);

  prog.end();

  return 0;
}
