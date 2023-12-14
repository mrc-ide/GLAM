
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "misc.h"

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
list mcmc_cpp(cpp11::sexp rng_ptr) {
  
  print("Running C++ code");
  
  // dust initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  
  // make a random draw
  double rand1 = dust::random::random_real<double>(state);
  print(rand1);
  
  // Return outputs in a list
  return writable::list({
    "foo"_nm = -9,
    "rng_ptr"_nm = rng_ptr
  });
}
