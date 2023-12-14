#include <cpp11.hpp>
using namespace cpp11;

[[cpp11::register]]
int cpp_test() {
  return 5;
}
