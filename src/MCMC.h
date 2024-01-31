
#pragma once

#include <cpp11.hpp>
#include "Particle.h"
#include "misc.h"
#include "System.h"

//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  // PUBLIC OBJECTS
  
  // parameter values
  
  // pointer to system object
  System * sys;
  
  int n_rungs;
  double start_time;
  double end_time;
  int max_infections;
  
  std::vector<Particle> particle_vec;
  
  // proposal sd
  std::vector<std::vector<double>> proposal_sd_mat;
  int n_proposal_sd;
  
  // counters
  int iteration_counter;
  std::vector<std::vector<int>> acceptance_out;
  std::vector<int> swap_acceptance_out;
  
  // objects for storing results
  std::vector<double> lambda_store;
  std::vector<double> theta_store;
  std::vector<double> decay_rate_store;
  std::vector<double> sens_store;
  std::vector<std::vector<int>> n_infections_store;
  std::vector<std::vector<std::vector<double>>> infection_times_store;
  cpp11::writable::list param_list_out;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(System &sys,
            cpp11::list data_list,
            cpp11::list obs_time_list,
            cpp11::list param_list,
            cpp11::list proposal_sd,
            const int iteration_counter_init,
            const cpp11::doubles beta,
            const double start_time,
            const double end_time,
            int max_infections);
  
  void run_mcmc(bool burnin, int interations);
  
};
