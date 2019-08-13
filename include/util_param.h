//#ifndef UTIL_PARAM_H
//#define UTIL_PARAM_H

#include <string>

using namespace std;

void extract_param(string paramfile, string whichparam, string *param, bool *found);

void summarise_param(string whichparam, string param, bool found);

void to_lowercase(char& c);

//#endif