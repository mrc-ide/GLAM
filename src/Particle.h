
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"
#include "probability.h"
#include "Indiv.h"
#include "System.h"

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
  int max_infections;
  
  // proposal sd
  std::vector<double> proposal_sd_vec;
  int n_proposal_sd;
  
  // Indiv objects
  std::vector<Indiv> indiv_vec;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(System &sys,
            double lambda,
            double theta,
            double decay_rate,
            double sens,
            std::vector<int> n_infections,
            std::vector<std::vector<double>> infection_times,
            std::vector<double> proposal_sd,
            double beta,
            double start_time,
            double end_time,
            int max_infections);
  
  void update();
  void update_lambda();
  void update_theta();
  void update_decay_rate();
  void update_sens();
  
};
