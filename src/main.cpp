
#include <cpp11.hpp>
#include "misc.h"

using namespace cpp11;

[[cpp11::register]]
int mcmc_cpp(cpp11::sexp rng_ptr) {
  print("foo");
  
  return 5;
}
