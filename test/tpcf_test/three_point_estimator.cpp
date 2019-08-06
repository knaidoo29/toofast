#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "toofast.h"

using namespace std;

int get_triangles(vector<float> &r1, vector<float> &r2, vector<float> &r3,
  float minimum, float maximum, int numbins, bool uselog){

  float dx, log10min, log10max, log10dx;

  dx = (maximum - minimum)/((float)(numbins));
  log10min = log10(minimum);
  log10max = log10(maximum);
  log10dx = (log10max - log10min)/((float)(numbins));

  float r[numbins];

  for(int i = 0; i < numbins; i++){
    if(uselog == true){
      r[i] = pow(10., log10min + log10dx/2. + float(i)*log10dx);
    }
    else{
      r[i] = minimum + dx/2. + float(i)*dx;
    }
  }

  float r1_val, r2_val, r3_val;
  int n_triangles = 0;

  for(int i = 0; i < numbins; i++){
    for(int j = 0; j < numbins; j++){
      for(int k = 0; k < numbins; k++){
        r1_val = r[i];
        r2_val = r[j];
        r3_val = r[k];
        if(r1_val <= r2_val and r2_val <= r3_val and r1_val + r2_val >= r3_val){
          r1.push_back(r1_val);
          r2.push_back(r2_val);
          r3.push_back(r3_val);
          n_triangles += 1;
        }
      }
    }
  }
  return n_triangles;
}

int get_triangle_index(vector<float> &r1, vector<float> &r2, vector<float> &r3,
  float minimum, float maximum, int numbins, bool uselog, float r1_val,
  float r2_val, float r3_val){

    float dx, log10min, log10max, log10dx;

    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    int ind;
    float diff_r1, diff_r2, diff_r3;

    for(int i = 0; i < r1.size(); i++){
      if(uselog == false){
        diff_r1 = abs((r1_val - r1[i])/dx);
        diff_r2 = abs((r2_val - r2[i])/dx);
        diff_r3 = abs((r3_val - r3[i])/dx);
        if(diff_r1 <= 0.5 and diff_r2 <= 0.5 and diff_r3 <= 0.5){
          ind = i;
        }
      }
      else{
        diff_r1 = abs((log10(r1_val) - log10(r1[i]))/log10dx);
        diff_r2 = abs((log10(r2_val) - log10(r2[i]))/log10dx);
        diff_r3 = abs((log10(r3_val) - log10(r3[i]))/log10dx);
        if(diff_r1 <= 0.5 and diff_r2 <= 0.5 and diff_r3 <= 0.5){
          ind = i;
        }
      }
    }

    return ind;
}

long int get_ddd_count(int total){
  long int count = 0;
  for(long int i = 0; i < total; i++){
    count += i*(i-1)/2;
  }
  return count;
}

long int get_ddr_count(int total){
  long int count;
  count = total*total*(total-1)/2;
  return count;
}

void get_ddd(vector<float> &r1, vector<float> &r2, vector<float> &r3,
  float minimum, float maximum, int numbins, bool uselog, vector<float> &x,
  vector<float> &y, vector<float> &z, int n_triangles, float ddd[]){

  for(int i = 0; i < n_triangles; i++){
    ddd[i] = 0.;
  }

  int ind, size;
  long int count = 0, count_total;
  float r123[3], r12, r23, r31;

  size = x.size()/50;
  count_total = get_ddd_count(size);

  for(long int i = 0; i < size; i++){
    for(long int j = i + 1; j < size; j++){
      for(long int k = j + 1; k < size; k++){
        r12 = get_distance_3d(x[i], y[i], z[i], x[j], y[j], z[j]);
        r23 = get_distance_3d(x[j], y[j], z[j], x[k], y[k], z[k]);
        r31 = get_distance_3d(x[k], y[k], z[k], x[i], y[i], z[i]);
        r123[0] = r12;
        r123[1] = r23;
        r123[2] = r31;
        sort(r123,r123+3);
        r12 = r123[0];
        r23 = r123[1];
        r31 = r123[2];
        ind = get_triangle_index(r1, r2, r3, minimum, maximum, numbins, uselog, r12, r23, r31);
        ddd[ind] += 1.;
        progress_bar(count, count_total, true);
        if(count >= 32340){
          cout << ind << '\t' << count << '\t' << r12 << '\t' << r23 << '\t' << r31 << endl;  
        }
        count += 1;
      }
    }
  }

  for(int i = 0; i < n_triangles; i++){
    ddd[i] *= 3.;
  }

}

int main(){

  // ------------------------------------------------------------------------ //

  string filename_data, filename_rand, filename_output;
  int numbins;
  float minimum, maximum;
  bool uselog;

  // Paramfile -- Proxy

  filename_data = "data.txt";
  filename_rand = "randoms.txt";
  filename_output = "xi.txt";

  minimum = 0.01;
  maximum = 75.;
  numbins = 10;

  uselog = false;

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

  vector<float>r1;
  vector<float>r2;
  vector<float>r3;

  int n_triangles;

  n_triangles = get_triangles(r1, r2, r3, minimum, maximum, numbins, uselog);

  cout << '\n' << '\t' << "Number of triangles = " << n_triangles << endl;

  bool print_triangles = false;

  if(print_triangles == true){
    cout << '\n';
    cout << '\t' << "r1" << '\t' << "r2" << '\t' << "r3" << endl;
    for(int i = 0; i < n_triangles; i++){
      cout << '\t' << r1[i] << '\t' << r2[i] << '\t' << r3[i] << endl;
    }
  }

  prog.start_process("Calculate Xi(r)", 2);

  prog.start_process("Get DDD", 3);

  float ddd[n_triangles];

  get_ddd(r1, r2, r3, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
    n_triangles, ddd);

  // ------------------------------------------------------------------------ //
  /*
  prog.start_process("Get RR", 3);

  float rr[numbins];

  get_dd(rr, minimum, maximum, numbins, uselog, x_rand, y_rand, z_rand);
  */
  // ------------------------------------------------------------------------ //
  /*
  prog.start_process("Get DR", 3);

  float dr[numbins];

  get_dr(dr, minimum, maximum, numbins, uselog, x_data, y_data, z_data,
    x_rand, y_rand, z_rand);
  */
  // ------------------------------------------------------------------------ //
  /*
  prog.start_process("Print xi(r)", 3);

  float r[numbins];

  get_r(r, minimum, maximum, numbins, uselog);

  float xi[numbins];

  get_xi(dd, dr, rr, numbins, x_data.size(), x_rand.size(), xi);

  cout << '\t' << '\t' << "r" << '\t' << "Xi(r)" << endl;

  for(int i = 0; i < numbins; i++){
    cout << '\t' << '\t' <<  r[i] << '\t' << xi[i] << endl;
  }
  */
  // ------------------------------------------------------------------------ //
  /*
  prog.start_process("Save data", 2);

  Writer wr;
  wr.add2header("Levy flight Xi estimation");
  wr.store(r, numbins, "r");
  wr.store(dd, numbins, "DD");
  wr.store(rr, numbins, "RR");
  wr.store(dr, numbins, "DR");
  wr.store(xi, numbins, "xi");
  wr.write2file(filename_output);
  */
  prog.end();

  return 0;
}
