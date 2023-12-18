
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"

//------------------------------------------------
// class defining particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // parameter values
  double lambda;
  double theta;
  double decay_rate;
  double sens;
  std::vector<int> n_infections;
  int n_samp;
  std::vector<std::vector<double>> infection_times;
  double start_time;
  double end_time;
  
  // proposal sd
  std::vector<double> proposal_sd_vec;
  int n_proposal_sd;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(double lambda,
            double theta,
            double decay_rate,
            double sens,
            std::vector<int> n_infections,
            std::vector<double> proposal_sd,
            double beta,
            double start_time,
            double end_time,
            cpp11::sexp rng_ptr);
  void update();
  
};
