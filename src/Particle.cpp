
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
                    std::vector<double> proposal_sd,
                    double beta,
                    cpp11::sexp rng_ptr) {
  
  this->lambda = lambda;
  this->theta = theta;
  
  proposal_sd_vec = proposal_sd;
  n_proposal_sd = proposal_sd_vec.size();
  
  // dust initialise RNG
  this->rng_ptr = rng_ptr;
  
  print(lambda, theta);
}

//------------------------------------------------
// update
void Particle::update() {
  
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  
  lambda = dust::random::random_real<double>(state);
  
}
