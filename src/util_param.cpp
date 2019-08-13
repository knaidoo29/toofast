#include "../include/twofast.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

void extract_param(string paramfile, string whichparam, string *param, bool *found){
  /* Extracts a parameter from a parameter file.

  Parameters
  ----------
  string paramfile:
    Name of the parameter file.
  string whichparam:
    Name of the parameter we want to extract from the file.
  string param:
    Sets the param value according to what was found in the parameter file.
  bool found:
    Whether the parameter was found or not.
  */
  string _line;
  int check = 0;
  ifstream _file (paramfile);
  if(_file.is_open()){
    while(getline(_file, _line)){
      if(_line.find(whichparam) == 0){
        _line.erase(_line.find(whichparam), whichparam.length());
        remove(_line.begin(), _line.end(), ' ');
        if(_line.find('%') < _line.size()){
          _line.erase(_line.find('%'), _line.size()-_line.find('%'));
        }
        *param = _line;
        if(_line.size() == 0){
          *found = false;
        }
        else{
          *found = true;
        }
        check = 1;
      }
    }
  }
  if(check == 0){
    cout << "Parameter " << whichparam << " was not found in the parameter file " << paramfile << endl;
  }
}

void summarise_param(string whichparam, string param, bool found){
  /* Summarise whether a parameter was found or not.

  Parameters
  ----------
  string whichparam:
    Which parameter we're looking at.
  string param:
    Value for the parameter.
  bool found:
    Whether the parameter was found or not.
  */  
  if(found == true){
    cout << whichparam << " = " << param << "  (Found = Yes)" << endl;
  }
  else{
    cout << whichparam << " = " << param << "  (Found = No)" << endl;
  }
}

void to_lowercase(char& c){
  /* Converts string into lowercase.

  Parameters
  ----------
  char c:
    A string with a mix of upper and lower case symbols.
  */
  c = tolower(static_cast<unsigned char>(c));
}
