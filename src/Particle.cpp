
#include <chrono>
#include <thread>


#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "Particle.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// constructor
void Particle::init(double lambda,
                    double theta,
                    double decay_rate,
                    double sens,
                    std::vector<int> n_infections,
                    std::vector<double> proposal_sd,
                    double beta,
                    double start_time,
                    double end_time,
                    cpp11::sexp rng_ptr) {
  
  // initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  rng_state = rng->state(0);
  
  // copy over known values
  this->lambda = lambda;
  this->theta = theta;
  this->decay_rate = decay_rate;
  this->sens = sens;
  this->n_infections = n_infections;
  this->start_time = start_time;
  this->end_time = end_time;
  
  proposal_sd_vec = proposal_sd;
  n_proposal_sd = proposal_sd_vec.size();
  
  // draw starting infection times
  n_samp = n_infections.size();
  infection_times = std::vector<std::vector<double>>(n_samp);
  for (int i = 0; i < n_samp; ++i) {
    infection_times[i] = std::vector<double>(n_infections[i]);
    for (int j = 0; j < n_infections[i]; ++j) {
      infection_times[i][j] = start_time + (end_time - start_time) * dust::random::random_real<double>(rng_state);
    }
  }
  
}

//------------------------------------------------
// update
void Particle::update() {
  
  lambda = dust::random::random_real<double>(rng_state);
  theta = dust::random::random_real<double>(rng_state);
  
  for (int i = 0; i < n_samp; ++i) {
    for (int j = 0; j < n_infections[i]; ++j) {
      infection_times[i][j] = start_time + (end_time - start_time) * dust::random::random_real<double>(rng_state);
    }
  }
  
}
