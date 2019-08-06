#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

using namespace std;

float get_distance_3d(float x1, float y1, float z1, float x2, float y2, float z2);

void get_r(float r[], float minimum, float maximum, int numbins, bool uselog);

void get_dd(float dd[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x, vector<float> &y, vector<float> &z);

void get_dr(float dr[], float minimum, float maximum, int numbins, bool uselog,
  vector<float> &x1, vector<float> &y1, vector<float> &z1,
  vector<float> &x2, vector<float> &y2, vector<float> &z2);

void get_xi(float dd[], float dr[], float rr[], int numbins, int nd, int nr, float xi[]);
