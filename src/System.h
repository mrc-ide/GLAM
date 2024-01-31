
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
  cpp11::list param_list;
  cpp11::list proposal_sd;
  int iteration_counter_init;
  cpp11::doubles beta;
  double start_time;
  double end_time;
  int max_infections;
  dust::random::xoshiro256plus& rng_state;
  
  int n_rungs;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  void init(cpp11::list data_list,
            cpp11::list obs_time_list,
            cpp11::list param_list,
            cpp11::list proposal_sd,
            int iteration_counter_init,
            cpp11::doubles beta,
            double start_time,
            double end_time,
            int max_infections,
            cpp11::sexp rng_ptr);
  
};