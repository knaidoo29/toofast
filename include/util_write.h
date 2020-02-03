//#ifndef UTIL_WRITE_H
//#define UTIL_WRITE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

class Writer{
  /* A class for writing to a file. */
  public:
    Writer(){
      store_as = "undefined";
    }
    template<typename TYPE>
    void store(TYPE array[], int size, string description, string store_type="float");
    template<typename TYPE>
    void store(vector<TYPE> &array, int size, string description, string store_type="float");
    void add2header(string exp);
    void write2file(string filename);
    void clean();
  private:
    vector<int> data_int;
    vector<float> data_float;
    vector<double> data_double;
    vector<int> structure;
    vector<string> header;
    vector<string> explanation;
    vector<string> types;
    string store_as;
    void restructure();
    void redistribute(string conversion);
};

template<typename TYPE>
void Writer::store(TYPE array[], int size, string description, string store_type){
  /* Stores an array to a vector to be written to a file later.

  Parameters
  ----------
  type array:
    Array to be stored by the writter class.
  int size:
    Size of the array.
  string description:
    A description of the array being stored.
  string store_type:
    How to store the array.
  */
  if(store_as == "undefined"){
    store_as = store_type;
  }
  else{
    if(store_as == "int" && store_type != "int"){
      if(store_type == "float"){
        Writer::redistribute("int2float");
        store_as = store_type;
      }
      if(store_type == "double"){
        Writer::redistribute("int2double");
        store_as = store_type;
      }
    }
    if(store_as == "float" && store_type != "float"){
      if(store_type == "double"){
        Writer::redistribute("float2double");
        store_as = store_type;
      }
    }
  }
  if(store_type != "int" && store_type != "float" && store_type != "double"){
    cout << "ERROR : Type '" << store_type << "' is unsupported." << endl;
  }
  else{
    types.push_back(store_type);
  }
  structure.push_back(size);
  header.push_back(description);
  for(int i = 0; i < size; i++){
    if(store_as == "int"){
      data_int.push_back((int) array[i]);
    }
    if(store_as == "float"){
      data_float.push_back((float) array[i]);
    }
    if(store_as == "double"){
      data_double.push_back((double) array[i]);
    }
  }
}


template<typename TYPE>
void Writer::store(vector<TYPE> &array, int size, string description, string store_type){
  /* Stores a vector to a vector to be written to a file later.

  Parameters
  ----------
  type array:
    A vector to be stored by the writter class.
  int size:
    Size of the array.
  string description:
    A description of the array being stored.
  string store_type:
    How to store the array.
  */
  if(store_as == "undefined"){
    store_as = store_type;
  }
  else{
    if(store_as == "int" && store_type != "int"){
      if(store_type == "float"){
        Writer::redistribute("int2float");
        store_as = store_type;
      }
      if(store_type == "double"){
        Writer::redistribute("int2double");
        store_as = store_type;
      }
    }
    if(store_as == "float" && store_type != "float"){
      if(store_type == "double"){
        Writer::redistribute("float2double");
        store_as = store_type;
      }
    }
  }
  if(store_type != "int" && store_type != "float" && store_type != "double"){
    cout << "ERROR : Type '" << store_type << "' is unsupported." << endl;
  }
  else{
    types.push_back(store_type);
  }
  structure.push_back(size);
  header.push_back(description);
  for(int i = 0; i < size; i++){
    if(store_as == "int"){
      data_int.push_back((int) array[i]);
    }
    if(store_as == "float"){
      data_float.push_back((float) array[i]);
    }
    if(store_as == "double"){
      data_double.push_back((double) array[i]);
    }
  }
}


template<class InputIterator, class InputIterator2>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              int precision1=10, int precision2=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2; ++begin1, ++begin2)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              int precision1=10, int precision2=10, int precision3=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3; ++begin1, ++begin2, ++begin3)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4; ++begin1, ++begin2, ++begin3, ++begin4)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\n';
}

template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4, class InputIterator5>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              InputIterator5 begin5, InputIterator5 end5,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10, int precision5=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4 and begin5 != end5; ++begin1, ++begin2, ++begin3, ++begin4, ++begin5)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\t'
      << std::setprecision(precision5) << *begin5 << '\n';
}


template<class InputIterator, class InputIterator2, class InputIterator3, class InputIterator4, class InputIterator5, class InputIterator6>
void writedat(string filename,
              InputIterator begin1, InputIterator end1,
              InputIterator2 begin2, InputIterator2 end2,
              InputIterator3 begin3, InputIterator3 end3,
              InputIterator4 begin4, InputIterator4 end4,
              InputIterator5 begin5, InputIterator5 end5,
              InputIterator6 begin6, InputIterator6 end6,
              int precision1=10, int precision2=10, int precision3=10, int precision4=10, int precision5=10, int precision6=10){
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  f.open(filename);
  for ( ; begin1 != end1 and begin2 != end2 and begin3 != end3 and begin4 != end4 and begin5 != end5 and begin6 != end6; ++begin1, ++begin2, ++begin3, ++begin4, ++begin5, ++begin6)
    f << std::setprecision(precision1) << *begin1 << '\t'
      << std::setprecision(precision2) << *begin2 << '\t'
      << std::setprecision(precision3) << *begin3 << '\t'
      << std::setprecision(precision4) << *begin4 << '\t'
      << std::setprecision(precision5) << *begin5 << '\t'
      << std::setprecision(precision6) << *begin6 << '\n';
}

//#endif
