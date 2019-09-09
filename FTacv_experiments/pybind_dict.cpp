#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/pybind11.h>
#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/stl.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
namespace py = pybind11;
using namespace std;
template <typename V>
V get(py::dict m, const std::string &key, const V &defval) {
    cout<<defval<<"\n";
    return m[key.c_str()].cast<V>();
}
int dict_reader(py::dict params){
  float test_val=get(params,std::string("test_val"),35.0);
  cout<<test_val<<"\n";
  for (auto item : params){
    cout<<item.first<<"\n";
    cout<<item.second<<"\n";
    }
  return 0;
}
int main () {
  return 0;
}
PYBIND11_MODULE(pybind_dicts, m) {
  	m.def("dict_reader", &dict_reader, "read a dict value");
}
