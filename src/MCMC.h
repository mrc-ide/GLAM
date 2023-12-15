
#pragma once

#include <cpp11.hpp>
#include "misc.h"

//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  // PUBLIC OBJECTS
  
  // parameter values
  
  int n_rungs;
  
  // proposal sd
  std::vector<std::vector<double>> proposal_sd_mat;
  int n_proposal_sd;
  
  // counters
  int iteration_counter;
  std::vector<std::vector<int>> acceptance_out;
  std::vector<int> swap_acceptance_out;
  
  // RNG
  cpp11::sexp rng_ptr;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(cpp11::list param_list,
       cpp11::list proposal_sd,
       const int iteration_counter_init,
       const cpp11::doubles beta,
       cpp11::sexp rng_ptr);
  
  // member functions
  void run_mcmc(bool burnin, int interations);
  
};
