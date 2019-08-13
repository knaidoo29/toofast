//#ifndef UTIL_PRINT_H
//#define UTIL_PRINT_H

#include <iostream>
#include <vector>

using namespace std;

template<typename TYPE>
void print_array(TYPE array[], int size, bool list=true){
    // Prints the array.
    cout << "\n";
    for (int i = 0; i < size; i++)
    {
        if(list == true) {
            // Prints as a list.
            if(i < size-1) {
                cout << array[i] << ", ";
            }
            else{
                cout << array[i] << endl;
            }
        }
        else{
            // Prints each element on each line.
            cout << array[i] << endl;
        }
    }
}

template<class TYPE>
void print_vector(vector<TYPE> &array, bool list=true){
    // Prints the array.
    cout << "\n";
    for (int i = 0; i < array.size(); i++)
    {
        if(list == true) {
            // Prints as a list.
            if(i < array.size()-1) {
                cout << array[i] << ", ";
            }
            else{
                cout << array[i] << endl;
            }
        }
        else{
            // Prints each element on each line.
            cout << array[i] << endl;
        }
    }
}

//#endif