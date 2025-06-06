
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"

//------------------------------------------------
// class for passing fixed values around easily
class System {
  
public:
  // PUBLIC OBJECTS
  
  cpp11::list data_list;
  cpp11::list obs_time_list;
  cpp11::doubles obs_time_start;
  cpp11::doubles obs_time_end;
  cpp11::list param_list;
  cpp11::list param_update_list;
  cpp11::list proposal_sd;
  cpp11::function lambda_prior;
  cpp11::function theta_prior;
  cpp11::function decay_rate_prior;
  cpp11::function sens_prior;
  int iteration_counter_init;
  cpp11::doubles beta;
  int max_infections;
  dust::random::xoshiro256plus& rng_state;
  bool interactive;
  
  int n_samp;
  int n_rungs;
  bool lambda_fixed;
  bool theta_fixed;
  bool decay_rate_fixed;
  bool sens_fixed;
  bool n_infections_fixed;
  bool infection_times_fixed;
  bool w_list_fixed;
  
  std::vector<std::vector<std::vector<bool>>> data_bool;
  std::vector<std::vector<double>> obs_time_vec;
  std::vector<double> haplo_freqs;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  //System(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  System(dust::random::xoshiro256plus& rng_state) :
    lambda_prior(R_NilValue),
    theta_prior(R_NilValue),
    decay_rate_prior(R_NilValue),
    sens_prior(R_NilValue),
    rng_state(rng_state) {}
  
  void init(cpp11::list data_list,
            cpp11::list obs_time_list,
            cpp11::doubles obs_time_start,
            cpp11::doubles obs_time_end,
            const cpp11::doubles haplo_freqs,
            cpp11::list param_list,
            cpp11::list param_update_list,
            cpp11::list proposal_sd,
            cpp11::function lambda_prior,
            cpp11::function theta_prior,
            cpp11::function decay_rate_prior,
            cpp11::function sens_prior,
            int iteration_counter_init,
            cpp11::doubles beta,
            int max_infections,
            cpp11::sexp rng_ptr,
            bool interactive);
  
};
