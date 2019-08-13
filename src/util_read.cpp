#include "../include/twofast.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

void print_file(string filename){
  /* Prints the contents of a file line-by-line.

  Parameters
  ----------
  string filename:
    The name of the file.
  */
  string line;
  ifstream _file (filename);
  cout << "Reading ascii file: '" << filename << "':" << endl;
  if(_file.is_open()){
    while(getline(_file, line)){
      cout << line << endl;
    }
    _file.close();
  }
  else{
    cout << "Unable to open file " << filename << endl;
  }
}

void Reader::readfile(string filename){
  /* Reads a file.

  Parameters
  ----------
  string filename:
    The name of the file.
  */
  Reader::readheader(filename);
  for(int i = 0; i < types.size(); i++){
    if(i == 0){
      store_as = types[i];
    }
    else{
      if(store_as != types[i]){
        if(store_as == "int" && types[i] == "float"){
          store_as = types[i];
        }
        if(store_as == "int" && types[i] == "double"){
          store_as = types[i];
        }
        if(store_as == "float" && types[i] == "double"){
          store_as = types[i];
        }
      }
    }
  }
  cout << store_as << endl;
  print_vector(header);
  if(store_as == "int"){
    read_ascii_table(filename, data_int, false);
  }
  if(store_as == "float"){
    read_ascii_table(filename, data_float, false);
  }
  if(store_as == "double"){
    read_ascii_table(filename, data_double, false);
  }
}

void Reader::readheader(string filename){
  /* Reads the header of a file, where the header is assumed to be any line beginning with '#'.

  Parameters
  ----------
  string filename:
    The name of the file.
  */
  string line;
  int check_var, check_type;
  string var_str, type_str;
  var_str = "# VARIABLES: ";
  type_str = "# TYPES: ";
  ifstream _file (filename);
  cout << "Reading ascii file: '" << filename << "':" << endl;
  if(_file.is_open()){
    while(getline(_file, line)){
      if(line[0] == '#'){
        cout << line << endl;
        check_var = line.find(var_str);
        if(check_var != -1){
          line.erase(0, var_str.size());
          replace(line.begin(), line.end(), '\t', ' ');
          stringstream str2float(line);
          string temp_str;
          while(str2float >> temp_str){
            temp_str.erase(temp_str.size()-1, 1);
            header.push_back(temp_str);
          }
        }
        check_var = -1;
        check_type = line.find(type_str);
        if(check_type != -1){
          line.erase(0, type_str.size());
          replace(line.begin(), line.end(), '\t', ' ');
          stringstream str2float(line);
          string temp_str;
          while(str2float >> temp_str){
            temp_str.erase(temp_str.size()-1, 1);
            types.push_back(temp_str);
          }
        }
        check_type = -1;
      }
    }
  }
}
