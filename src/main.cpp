
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
              const doubles haplo_freqs,
              const int iterations,
              const bool burnin,
              list param_list,
              list param_update_list,
              list proposal_sd,
              const int iteration_counter_init,
              const doubles beta,
              double start_time,
              double end_time,
              int max_infections,
              cpp11::sexp rng_ptr) {
  
  // start timer
  std::chrono::high_resolution_clock::time_point t0 =  std::chrono::high_resolution_clock::now();
  
  // initialise RNG state
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  dust::random::xoshiro256plus& rng_state = rng->state(0);
  
  
  // TODO - delete
  // minimal failing example
  
  std::vector<double> z(3, 8.5);
  cpp11::writable::list ret1;
  cpp11::writable::doubles v0(z.begin(), z.end());
  ret1.push_back(v0);
  
  stop("early stop 2");
  
  //list tmp1 = param_list[0];
  //list tmp4 = tmp1["infection_times"];
  
  Particle p = Particle(rng_state);
  //p.infection_times = list_to_mat_double(tmp4);
  
  std::vector<std::vector<double>> x(1);
  x[0] = std::vector<double>(3, 9.5);
  p.infection_times = x;
  
  print_matrix(p.infection_times);
  
  cpp11::writable::list ret0;
  //for (int i = 0; i < x.size(); ++i) {
    cpp11::doubles v = cpp11::as_sexp(x[0]);
    //ret0.push_back(v);
  //}
  
  //cpp11::writable::list tmp = mat_double_to_list(p.infection_times);
  //cpp11::writable::list tmp = mat_double_to_list(x);
  
  //std::cout << "TYPEOF(tmp): " << TYPEOF(tmp) << std::endl;
  
  stop("early stop");
  
  
  //stop("foo1");
  
  // create system object to hold all input values. This makes it easier when
  // passing objects around between classes
  System sys(rng_state);
  sys.init(data_list,
           obs_time_list,
           haplo_freqs,
           param_list,
           param_update_list,
           proposal_sd,
           iteration_counter_init,
           beta,
           start_time,
           end_time,
           max_infections,
           rng_ptr);
  
  //stop("foo1");
  
  // initialise MCMC
  MCMC mcmc(sys.rng_state);
  
  //stop("foo1");
  
  mcmc.init(sys,
            param_list,
            proposal_sd,
            beta);
  
  //stop("foo1");
  
  // run main loop
  mcmc.run_mcmc(true, iterations);
  
  stop("foo1");
  
  // get output objects into cpp11 format
  writable::list acceptance_out = mat_int_to_list(mcmc.acceptance_out);
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
    "acceptance_out"_nm = acceptance_out,
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
                     const doubles haplo_freqs,
                     list param_list,
                     list param_update_list,
                     list proposal_sd,
                     const doubles beta,
                     double start_time,
                     double end_time,
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
           haplo_freqs,
           param_list,
           param_update_list,
           proposal_sd,
           0,
           beta,
           start_time,
           end_time,
           max_infections,
           rng_ptr);
  
  // initialise MCMC, which will run the debugging step
  MCMC mcmc(sys.rng_state);
  mcmc.init(sys,
            param_list,
            proposal_sd,
            beta);
}
