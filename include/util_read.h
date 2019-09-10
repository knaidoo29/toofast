//#ifndef UTIL_READ_H
//#define UTIL_READ_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

void print_file(string filename);

template<typename TYPE>
void read_ascii_table(string filename, int rows, int columns, TYPE array[], bool printout=false){
  /* Reads ascii table files, storing the data in an array and finding out how many rows and columns
  the file has.

  Parameters
  ----------
  string filename:
    Ascii filename.
  int rows:
    Number of rows in the file (to be assigned).
  int columns:
    Number of columns in the file (to be assigned).
  type array:
    Array to store the values of the ascii file.
  bool printout:
    Whether to print statements when opening ascii files.
  */
  string line;
  int i = 0;

  ifstream _file (filename);
  if(printout == true){
    cout << "\n\tReading ascii file: '" << filename << "':" << endl;
  }
  if(_file.is_open()){
    while(getline(_file, line)){
      replace(line.begin(), line.end(), '\t', ' ');
      stringstream str2float(line);
      TYPE temp_float;
      while(str2float >> temp_float){
        array[i] = temp_float;
        i++;
      }
    }
    _file.close();
  }
  else{
    if(printout == true){
      cout << "Unable to open file " << filename << endl;
    }
  }
}

template<typename TYPE>
int read_ascii_table(string filename, vector<TYPE> &array, bool printout=false){
  /* Reads ascii table files, storing the data in a vector and finding out how many rows and columns
  the file has.

  Parameters
  ----------
  string filename:
    Ascii filename.
  type array:
    Vector array to store the values of the ascii file.
  bool printout:
    Whether to print statements when opening ascii files.
  */
  string line;
  int count_col, count_row;
  count_col = 0;
  count_row = 0;

  ifstream _file (filename);
  if(printout == true){
    cout << "\n\tReading ascii file: '" << filename << "':" << endl;
  }
  if(_file.is_open()){
    while(getline(_file, line)){
      if(line[0] != '#'){
        replace(line.begin(), line.end(), '\t', ' ');
        stringstream str2float(line);
        TYPE temp_float;
        while(str2float >> temp_float){
          array.push_back(temp_float);
          if(count_row == 0){
            count_col++;
          }
        }
        count_row++;
      }
    }
    if(printout == true){
      cout << "\t- Number of columns = " << count_col << endl;
      cout << "\t- Number of rows = " << count_row << endl;
    }
    _file.close();
  }
  else{
    if(printout == true){
      cout << "Unable to open file " << filename << endl;
    }
  }
  return count_col;
}


template<typename TYPE>
void extract_from_table(vector<TYPE> &table, int col, int col_length, vector<TYPE> &array){
  /* Extracts a column from a vector.

  Paramater
  ---------
  type table:
    The table of stored data.
  int col:
    Column we want data from.
  int col_length:
    Length of the column we want to extract from the vector table.
  type array:
    Vector array containing the information in the selected column of the table
  */
  for(int i = 0; i < table.size()/col_length; i++){
    array.push_back(table[i*col_length + col]);
  }
}

class Reader{
  /* A class for reading ascii table data. */
  public:
    void readfile(string filename);
    template<typename TYPE>
    void extract(string description, TYPE array[]);
    template<typename TYPE>
    void extract(string description, vector<TYPE> &array);
    vector<string> header;
    vector<string> types;
  private:
    int rows, cols;
    string store_as;
    vector<int> data_int;
    vector<float> data_float;
    vector<double> data_double;
    void readheader(string filename);
};

template<typename TYPE>
void Reader::extract(string description, vector<TYPE> &array){
  /* A function for extracting an array from a given header information.

  Parameters
  ----------
  string description:
    The header description of the column of interest.
  type array:
    The array where the extracted ascii data is stored.
  */
  int j=0;
  bool found;
  while(header[j] != description && j<header.size()){
    j++;
  }
  if(header[j] == description){
    found = true;
  }
  else{
    found = false;
    cout << "ERROR: Variable does not exist" << endl;
    cout << "Variables that are stored are :";
  }
  if(store_as == "int"){
    for(int i=0; i < data_int.size()/header.size(); i++){
      array.push_back((TYPE) data_int[header.size()*i + j]);
    }
  }
  if(store_as == "float"){
    for(int i=0; i < data_float.size()/header.size(); i++){
      array.push_back((TYPE) data_float[header.size()*i + j]);
    }
  }
  if(store_as == "double"){
    for(int i=0; i < data_double.size()/header.size(); i++){
      array.push_back((TYPE) data_double[header.size()*i + j]);
    }
  }
}

//#endif
