#include "../include/twofast.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

void Writer::add2header(string exp){
  /* Adds to the header of the file.

  Parameters
  ----------
  string exp:
    An explanation to be added to the header of a file.
  */
  explanation.push_back(exp);
}

void Writer::write2file(string filename){
  /* Writes to a file called "filename".

  Parameters
  ----------
  string filename:
    Name of file to write to.
  */
  ofstream _file;
  _file.open(filename);
  int size = structure[0];
  int num = structure.size();
  Writer::restructure();
  string line;
  line = "# ASCII DATA FILE: COLUMNS = " + to_string(size) + ", ROWS= " + to_string(num) + "\n";
  _file << line;
  if(explanation.size() != 0){
    for(int i = 0; i < explanation.size(); i++){
      line = "# " + explanation[i] + "\n";;
      _file << line;
    }
  }
  line = "# VARIABLES: ";
  for(int i = 0; i < num; i++){
    line += header[i]+",\t";
  }
  _file << line + "\n";
  line = "# TYPES: ";
  for(int i = 0; i < num; i++){
    line += types[i]+",\t";
  }
  _file << line + "\n";
  for(int i = 0; i < size; i++){
    for(int j = 0; j < num; j++){
      if(j == 0){
        if(store_as == "int"){
          line = to_string(data_int[i + j*size]) + "\t";
        }
        if(store_as == "float"){
          line = to_string(data_float[i + j*size]) + "\t";
        }
        if(store_as == "double"){
          line = to_string(data_double[i + j*size]) + "\t";
        }
      }
      else{
        if(store_as == "int"){
          line += to_string(data_int[i + j*size]) + "\t";
        }
        if(store_as == "float"){
          line += to_string(data_float[i + j*size]) + "\t";
        }
        if(store_as == "double"){
          line += to_string(data_double[i + j*size]) + "\t";
        }
      }
    }
    _file << line + "\n";
  }
  _file.close();
}

void Writer::restructure(){
  /* Forgotten what this does... restructures something but don't remember what...
  this is why you document stuff as you write it. Maybe to do with storing arrays of
  different lengths. Apparently anticipated a problem that I didn't even remember it could handle.
  */
  int count;
  for(int i = 0; i < structure.size(); i++){
    if(i == 0){
      count = structure[0];
    }
    else{
      count += structure[i];
      structure[i] = count;
    }
  }
}

void Writer::redistribute(string conversion){
  /* Converts data types.

  Parameters
  ----------
  string conversion:
    Conversion of data type to be written as a string. This always increases precision,
    i.e. int to float/double or float to double.
  */
  if(conversion == "int2float"){
    for(int i = 0; i < data_int.size(); i++){
      data_float.push_back((float) data_int[i]);
    }
    data_int.erase(data_int.begin(), data_int.begin()+data_int.size());
  }
  if(conversion == "int2double"){
    for(int i = 0; i < data_int.size(); i++){
      data_double.push_back((double) data_int[i]);
    }
    data_int.erase(data_int.begin(), data_int.begin()+data_int.size());
  }
  if(conversion == "float2double"){
    for(int i = 0; i < data_float.size(); i++){
      data_double.push_back((double) data_float[i]);
    }
    data_float.erase(data_float.begin(), data_float.begin()+data_float.size());
  }
}

void Writer::clean(){
  /* Writer class internal cleaner function. */
  if(data_int.size() != 0){
    data_int.erase(data_int.begin(), data_int.begin()+data_int.size());
  }
  if(data_float.size() != 0){
    data_float.erase(data_float.begin(), data_float.begin()+data_float.size());
  }
  if(data_double.size() != 0){
    data_double.erase(data_double.begin(), data_double.begin()+data_double.size());
  }
  if(structure.size() != 0){
    structure.erase(structure.begin(), structure.begin()+structure.size());
  }
  if(header.size() != 0){
    header.erase(header.begin(), header.begin()+header.size());
  }
  if(explanation.size() != 0){
    explanation.erase(explanation.begin(), explanation.begin()+explanation.size());
  }
  if(types.size() != 0){
    types.erase(types.begin(), types.begin()+types.size());
  }
  store_as = "undefined";
}
