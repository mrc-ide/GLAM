
#include <cpp11.hpp>
#include "misc.h"

using namespace cpp11;

[[cpp11::register]]
int cpp_test() {
  
  
  std::vector<int> x = seq_int(1, 5);
  print_vector(x);
  
  return 5;
}
