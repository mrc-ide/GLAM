
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "misc.h"

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
list mcmc_cpp(const int iterations,
              const bool burnin,
              const int iteration_counter_init,
              list proposal_sd,
              const doubles beta,
              cpp11::sexp rng_ptr) {
  
  print("Running C++ code");
  
  // start timer
  std::chrono::high_resolution_clock::time_point t0 =  std::chrono::high_resolution_clock::now();
  
  // extract proposal sd
  std::vector<std::vector<double>> proposal_sd_array = list_to_mat_double(proposal_sd);
  int n_proposal_sd = proposal_sd_array[0].size();
  
  const int n_rung = beta.size();
  
  // dust initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  auto& state = rng->state(0);
  
  // counters
  int iteration_counter = iteration_counter_init + 1;
  writable::integers_matrix<> acceptance_out(n_rung, n_proposal_sd);
  cpp11_init(acceptance_out, 0);
  writable::integers swap_acceptance_out(n_rung - 1);
  cpp11_init(swap_acceptance_out, 0);
  
  // dummy MCMC
  int start_i = 0;
  if (burnin && iteration_counter_init == 0) {
    start_i = 1;
    iteration_counter++;
  }
  for (int i = start_i; i < iterations; ++i) {
    
    // mock updates
    for (int r = 0; r < n_rung; ++r) {
      for (int j = 0; j < n_proposal_sd; ++j) {
        double rand1 = dust::random::random_real<double>(state);
        if (rand1 < 0.5) {
          acceptance_out(r,j)++;
        }
      }
    }
    
    // mock Metropolis coupling
    for (int r = 0; r < (n_rung - 1); ++r) {
      swap_acceptance_out[r]++;
    }
    
    iteration_counter++;
    
  }  // close main MCMC loop
  
  // calculate elapsed time
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
  double dur = time_span.count();
  
  // return outputs in a list
  return writable::list({
    "acceptance_out"_nm = acceptance_out,
    "swap_acceptance_out"_nm = swap_acceptance_out,
    "dur"_nm = dur,
    "rng_ptr"_nm = rng_ptr
  });
}
