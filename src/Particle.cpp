
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
  
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  
  // copy over known values
  this->lambda = lambda;
  this->theta = theta;
  this->decay_rate = decay_rate;
  this->sens = sens;
  this->n_infections = n_infections;
  
  proposal_sd_vec = proposal_sd;
  n_proposal_sd = proposal_sd_vec.size();
  
  // draw starting infection times
  n_samp = n_infections.size();
  infection_times = std::vector<std::vector<double>>(n_samp);
  for (int i = 0; i < n_samp; ++i) {
    infection_times[i] = std::vector<double>(n_infections[i]);
    for (int j = 0; j < n_infections[i]; ++j) {
      infection_times[i][j] = start_time + (end_time - start_time) * dust::random::random_real<double>(state);
    }
  }
  
  // dust initialise RNG
  this->rng_ptr = rng_ptr;
}

//------------------------------------------------
// update
void Particle::update() {
  
  //std::this_thread::sleep_for(std::chrono::milliseconds(1));
  
  //auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  //dust::random::xoshiro256plus& state = rng->state(0);
  //auto& state = rng->state(0);
  
  //lambda = dust::random::random_real<double>(state);
  //theta = dust::random::random_real<double>(state);
  
  double start_time = 0.0;
  double end_time = 10.0;
  
  for (int i = 0; i < n_samp; ++i) {
    for (int j = 0; j < n_infections[i]; ++j) {
      //infection_times[i][j] = start_time + (end_time - start_time) * dust::random::random_real<double>(state);
      //infection_times[i][j] = 10.0 * dust::random::random_real<double>(state);
    }
  }
  
}
