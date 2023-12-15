
#include <cpp11.hpp>
#include "MCMC.h"
#include "misc.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
[[cpp11::register]]
list mcmc_cpp(const int iterations,
              const bool burnin,
              list param_list,
              list proposal_sd,
              const int iteration_counter_init,
              const doubles beta,
              cpp11::sexp rng_ptr) {
  
  print("Running C++ code");
  
  // start timer
  std::chrono::high_resolution_clock::time_point t0 =  std::chrono::high_resolution_clock::now();
  
  // initialise MCMC
  MCMC mcmc(param_list,
            proposal_sd,
            iteration_counter_init,
            beta,
            rng_ptr);
  
  // run main loop
  mcmc.run_mcmc(true, 100);
  
  // get output objects into cpp11 format
  writable::list acceptance_out;
  for (int i = 0; i < mcmc.acceptance_out.size(); ++i) {
    acceptance_out.push_back({""_nm = mcmc.acceptance_out[i]});
  }
  
  // calculate elapsed time
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
  double dur = time_span.count();
  
  // return outputs in a list
  return writable::list({
    "acceptance_out"_nm = acceptance_out,
    "swap_acceptance_out"_nm = mcmc.swap_acceptance_out,
    "dur"_nm = dur,
    "rng_ptr"_nm = rng_ptr
  });
}
