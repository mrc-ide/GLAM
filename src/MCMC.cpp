
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "MCMC.h"
//#include "Particle.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// constructor
MCMC::MCMC(cpp11::list param_list,
           cpp11::list proposal_sd,
           const int iteration_counter_init,
           const cpp11::doubles beta,
           cpp11::sexp rng_ptr) {
  
  n_rungs = beta.size();
  
  // extract proposal sd
  proposal_sd_mat = list_to_mat_double(proposal_sd);
  n_proposal_sd = proposal_sd_mat[0].size();
  
  // initialise counters
  iteration_counter = iteration_counter_init + 1;
  acceptance_out = std::vector<std::vector<int>>(n_rungs, std::vector<int>(n_proposal_sd));
  swap_acceptance_out = std::vector<int>(n_rungs - 1);
  
  // dust initialise RNG
  this->rng_ptr = rng_ptr;
  
}

//------------------------------------------------
// run MCMC loop
void MCMC::run_mcmc(bool burnin, int iterations) {
  
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  
  int start_i = 0;
  if (burnin && iteration_counter == 1) {
    start_i = 1;
    iteration_counter++;
  }
  for (int i = start_i; i < iterations; ++i) {
    
    // mock updates
    for (int r = 0; r < n_rungs; ++r) {
      for (int j = 0; j < n_proposal_sd; ++j) {
        double rand1 = dust::random::random_real<double>(state);
        if (rand1 < 0.5) {
          acceptance_out[r][j]++;
        }
      }
    }
    
    // mock Metropolis coupling
    for (int r = 0; r < (n_rungs - 1); ++r) {
      swap_acceptance_out[r]++;
    }
    
    iteration_counter++;
    
  }
  
}
