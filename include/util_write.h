//#ifndef UTIL_WRITE_H
//#define UTIL_WRITE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

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

//#endif
