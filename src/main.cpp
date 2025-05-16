
#include <cpp11.hpp>
#include "MCMC.h"
#include "Particle.h"
#include "misc.h"
#include "System.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// main MCMC function
[[cpp11::register]]
list mcmc_cpp(cpp11::list data_list,
              cpp11::list obs_time_list,
              const doubles obs_time_start,
              const doubles obs_time_end,
              const doubles haplo_freqs,
              const int iterations,
              const bool burnin,
              list param_list,
              list param_update_list,
              list proposal_sd,
              cpp11::function lambda_prior,
              cpp11::function theta_prior,
              cpp11::function decay_rate_prior,
              cpp11::function sens_prior,
              const int iteration_counter_init,
              const doubles beta,
              int max_infections,
              cpp11::sexp rng_ptr,
              const bool interactive) {
  
  // start timer
  std::chrono::high_resolution_clock::time_point t0 =  std::chrono::high_resolution_clock::now();
  
  // initialise RNG state
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  dust::random::xoshiro256plus& rng_state = rng->state(0);
  
  // create system object to hold all input values. This makes it easier when
  // passing objects around between classes
  System sys(rng_state);
  sys.init(data_list,
           obs_time_list,
           obs_time_start,
           obs_time_end,
           haplo_freqs,
           param_list,
           param_update_list,
           proposal_sd,
           lambda_prior,
           theta_prior,
           decay_rate_prior,
           sens_prior,
           iteration_counter_init,
           beta,
           max_infections,
           rng_ptr,
           interactive);
  
  // initialise MCMC
  MCMC mcmc(sys.rng_state);
  
  mcmc.init(sys,
            param_list,
            proposal_sd,
            beta);
  
  
  
  // run main loop
  mcmc.run_mcmc(burnin, iterations);
  
  // get output objects into cpp11 format
  writable::list n_infections = mat_int_to_list(mcmc.n_infections_store);
  writable::list infection_times = array_double_to_list(mcmc.infection_times_store);
  
  // calculate elapsed time
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
  double dur = time_span.count();
  
  // return outputs in a list
  return writable::list({
    "lambda"_nm = mcmc.lambda_store,
    "theta"_nm = mcmc.theta_store,
    "decay_rate"_nm = mcmc.decay_rate_store,
    "sens"_nm = mcmc.sens_store,
    "n_infections"_nm = n_infections,
    "infection_times"_nm = infection_times,
    "param_list_out"_nm = mcmc.param_list_out,
    "prop_sd_list_out"_nm = mcmc.prop_sd_list_out,
    "swap_acceptance_out"_nm = mcmc.swap_acceptance_out,
    "dur"_nm = dur,
    "rng_ptr"_nm = rng_ptr
  });
}

//------------------------------------------------
// debug algorithm 1
[[cpp11::register]]
void debug_algo1_cpp(cpp11::list data_list,
                     cpp11::list obs_time_list,
                     const doubles obs_time_start,
                     const doubles obs_time_end,
                     const doubles haplo_freqs,
                     list param_list,
                     list param_update_list,
                     list proposal_sd,
                     cpp11::function lambda_prior,
                     cpp11::function theta_prior,
                     cpp11::function decay_rate_prior,
                     cpp11::function sens_prior,
                     const doubles beta,
                     int max_infections,
                     cpp11::sexp rng_ptr) {
  
  // initialise RNG state
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  dust::random::xoshiro256plus& rng_state = rng->state(0);
  
  // create system object to hold all input values. This makes it easier when
  // passing objects around between classes
  System sys(rng_state);
  sys.init(data_list,
           obs_time_list,
           obs_time_start,
           obs_time_end,
           haplo_freqs,
           param_list,
           param_update_list,
           proposal_sd,
           lambda_prior,
           theta_prior,
           decay_rate_prior,
           sens_prior,
           0,
           beta,
           max_infections,
           rng_ptr,
           true);
  
  // initialise MCMC, which will run the debugging step
  MCMC mcmc(sys.rng_state);
  mcmc.init(sys,
            param_list,
            proposal_sd,
            beta);
}
