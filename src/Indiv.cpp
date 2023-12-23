
#include <chrono>
#include <thread>


#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "Indiv.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// initialise
void Indiv::init(std::vector<std::vector<bool>> data_bool,
                 std::vector<double> obs_times,
                 cpp11::sexp rng_ptr,
                 int n_infections,
                 std::vector<double> infection_times) {
  
  // copy over values
  data = data_bool;
  n_haplos = data.size();
  this->obs_times = obs_times;
  n_obs = data[0].size();
  this->n_infections = n_infections;
  this->infection_times = infection_times;
  
  // sanity checks on inputs
  if (infection_times.size() != n_infections) {
    stop("Error in Indiv: infection_times.size() does not match n_infections");
  }
  if (obs_times.size() != n_obs) {
    stop("Error in Indiv: obs_times.size() does not match n_obs");
  }
  
  // main alleles array. Initialise with all alleles introduced (must be at least one introduced)
  infection_alleles = std::vector<std::vector<bool>>(n_infections, std::vector<bool>(n_haplos, true));
  
  // initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  rng_state = rng->state(0);
  
}

//------------------------------------------------
// update number of infections
void Indiv::update_n_infections() {
  
  // dummy update step. Randomly drop or add infection
  if (dust::random::random_real<double>(rng_state) < 0.5) {
    print("foo");
  } else {
    print("bar");
  }
  
}

//------------------------------------------------
// update timings of all infections
void Indiv::update_infection_times() {
  
}

//------------------------------------------------
// return n_infections
int Indiv::get_n_infections() {
  return n_infections;
}

//------------------------------------------------
// return vector of infection times
std::vector<double> Indiv::get_infection_times() {
  return infection_times;
}
