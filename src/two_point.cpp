#include "../include/toofast.h"

float get_distance_3d(float x1, float y1, float z1, float x2, float y2, float z2){

  float dist;

  dist = sqrt(pow(x1 - x2, 2.) + pow(y1 - y2, 2.) + pow(z1 - z2, 2.));
  return dist;
}

void get_r(float r[], float minimum, float maximum, int numbins, bool uselog){

  float dx, log10min, log10max, log10dx;

  dx = (maximum - minimum)/((float)(numbins));
  log10min = log10(minimum);
  log10max = log10(maximum);
  log10dx = (log10max - log10min)/((float)(numbins));

  for(int i = 0; i < numbins; i++){
    if(uselog == true){
      r[i] = pow(10., log10min + log10dx/2. + float(i)*log10dx);
    }
    else{
      r[i] = minimum + dx/2. + float(i)*dx;
    }
  }
}

void get_dd(float dd[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z){

    long int size, ind;
    long int k = 0, kmax;
    float dx, log10min, log10max, log10dx, dist;

    size = x.size();
    kmax = (size*(size - 1))/2;
    dx = (maximum - minimum)/((float)(numbins));
    log10min = log10(minimum);
    log10max = log10(maximum);
    log10dx = (log10max - log10min)/((float)(numbins));

    for(int i = 0; i < numbins; i++){
      dd[i] = 0.;
    }

    for(long int i = 0; i < size; i++){
      for(long int j = i+1; j < size; j++){
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
        progress_bar(k, kmax);
        k += 1;
      }
    }

    for(int i = 0; i < numbins; i++){
      dd[i] *= 2.;
    }
}

void get_dr(float dr[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x1, vector<float> &y1, vector<float> &z1,
  vector<float> &x2, vector<float> &y2, vector<float> &z2){

    long int size1, size2, ind;
    long int k = 0, kmax;
    float dx, log10min, log10max, log10dx, dist;

    size1 = x1.size();
    size2 = x2.size();
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
        dist = get_distance_3d(x1[i], y1[i], z1[i], x2[j], y2[j], z2[j]);
        if (dist >= minimum and dist <= maximum){
          if(uselog == true){
            ind = floor((log10(dist) - log10min)/log10dx);
          }
          else{
            ind = floor((dist - minimum)/dx);
          }
          dr[ind] += 1.;
        }
        progress_bar(k, kmax);
        k += 1;
      }
    }
}

void get_xi(float dd[], float dr[], float rr[], int numbins, int nd, int nr, float xi[]){

  for(int i = 0; i < numbins; i++){
    xi[i] = 1. + pow(((float)(nr))/((float)(nd)), 2.)*dd[i]/rr[i] - 2.*(((float)(nr))/((float)(nd)))*dr[i]/rr[i];
  }
}
