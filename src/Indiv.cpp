
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
void Indiv::init(System &sys,
                 std::vector<std::vector<bool>> data_bool,
                 std::vector<double> obs_times,
                 double start_time,
                 double end_time,
                 int max_infections,
                 int n_infections,
                 std::vector<double> infection_times) {
  
  // initialise RNG
  rng_state = sys.rng_state;
  
  // copy over values
  data = data_bool;
  n_haplos = data.size();
  this->obs_times = obs_times;
  n_obs = data[0].size();
  this->n_infections = n_infections;
  this->infection_times = infection_times;
  this->start_time = start_time;
  this->end_time = end_time;
  this->max_infections = max_infections;
  
  // sanity checks on inputs
  if (infection_times.size() != n_infections) {
    stop("Error in Indiv: infection_times.size() does not match n_infections");
  }
  if (obs_times.size() != n_obs) {
    stop("Error in Indiv: obs_times.size() does not match n_obs");
  }
  
  // main alleles array. Initialise with all alleles introduced (must be at least one introduced)
  infection_alleles = std::vector<std::vector<bool>>(n_infections, std::vector<bool>(n_haplos, true));
  
}

//------------------------------------------------
// update number of infections
void Indiv::update_n_infections() {
  
  // dummy update step. Randomly drop or add infection
  if (runif1(rng_state) < 0.5) {  // add infection
    if (n_infections == max_infections) {
      return;
    }
    
    double new_time = runif1(rng_state, start_time, end_time);
    std::vector<bool> allele_vec(n_haplos, true);
    
    infection_alleles.push_back(allele_vec);
    infection_times.push_back(new_time);
    n_infections++;
    
  } else {  // drop infection
    if (n_infections == 0) {
      return;
    }
    
    int tmp1 = sample2(rng_state, 0, n_infections - 1);
    
    infection_alleles.erase(infection_alleles.begin() + tmp1);
    infection_times.erase(infection_times.begin() + tmp1);
    n_infections--;
    
  }
  
}

//------------------------------------------------
// update timings of all infections
void Indiv::update_infection_times() {
  
  // dummy update
  for (int i = 0; i < n_infections; ++i) {
    infection_times[i] = runif1(rng_state, start_time, end_time);
  }
  
}

//------------------------------------------------
// basic log-likelihood, no marginalisation
double Indiv::loglike_basic(double lambda, double theta, double decay_rate, double sens) {
  return 0.0;
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
