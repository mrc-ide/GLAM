
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include <RProgress.h>
#include "MCMC.h"
//#include "Particle.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// initialise
void MCMC::init(System &sys,
                cpp11::list param_list,
                cpp11::list proposal_sd,
                const cpp11::doubles beta) {
  
  // pointer to system object
  this->sys = &sys;
  
  // initialise RNG
  rng_state = sys.rng_state;
  
  n_rungs = beta.size();
  
  // extract proposal sd
  proposal_sd_mat = list_to_mat_double(proposal_sd);
  n_proposal_sd = proposal_sd_mat[0].size();
  
  // initialise particles
  particle_vec = std::vector<Particle>(n_rungs, Particle(rng_state));
  for (int r = 0; r < n_rungs; ++r) {
    
    // extract parameters
    cpp11::list tmp1 = param_list[r];
    doubles tmp2 = tmp1["lambda"];
    double lambda = tmp2[0];
    
    tmp2 = tmp1["theta"];
    double theta = tmp2[0];
    
    tmp2 = tmp1["decay_rate"];
    double decay_rate = tmp2[0];
    
    tmp2 = tmp1["sens"];
    double sens = tmp2[0];
    
    integers tmp3 = tmp1["n_infections"];
    std::vector<int> n_infections = integers_to_vec(tmp3);
    
    list tmp4 = tmp1["infection_times"];
    std::vector<std::vector<double>> infection_times = list_to_mat_double(tmp4);
    
    list tmp5 = tmp1["w_list"];
    int n_samp = tmp5.size();
    std::vector<std::vector<std::vector<bool>>> w_array(n_samp);
    for (int i = 0; i < n_samp; ++i) {
      list tmp5b = tmp5[i];
      w_array[i] = list_to_mat_bool(tmp5b);
    }
    
    // pass parameters into particle
    particle_vec[r].init(sys,
                         lambda,
                         theta,
                         decay_rate,
                         sens,
                         n_infections,
                         infection_times,
                         w_array,
                         proposal_sd_mat[r],
                         beta[r]);
  }
  
  // initialise counters
  iteration_counter = sys.iteration_counter_init + 1;
  acceptance_out = std::vector<std::vector<int>>(n_rungs, std::vector<int>(n_proposal_sd));
  swap_acceptance_out = std::vector<int>(n_rungs - 1);
  
}

//------------------------------------------------
// run MCMC loop
void MCMC::run_mcmc(bool burnin, int iterations) {
  
  // objects for storing results
  lambda_store = std::vector<double>(iterations);
  theta_store = std::vector<double>(iterations);
  decay_rate_store = std::vector<double>(iterations);
  sens_store = std::vector<double>(iterations);
  n_infections_store = std::vector<std::vector<int>>(iterations);
  infection_times_store = std::vector<std::vector<std::vector<double>>>(iterations);
  
  // initialise progress bar
  RProgress::RProgress progress("Progress [:bar] Time remaining: :eta");
  progress.set_total(iterations);
  
  // store values on first iteration
  int start_i = 0;
  if (burnin && iteration_counter == 1) {
    
    lambda_store[0] = particle_vec[0].lambda;
    theta_store[0] = particle_vec[0].theta;
    decay_rate_store[0] = particle_vec[0].decay_rate;
    sens_store[0] = particle_vec[0].sens;
    n_infections_store[0] = particle_vec[0].n_infections;
    infection_times_store[0] = particle_vec[0].infection_times;
    
    start_i = 1;
    iteration_counter++;
  }
  
  // run loop
  for (int i = start_i; i < iterations; ++i) {
    progress.tick();
    
    // update particles
    for (int r = 0; r < n_rungs; ++r) {
      particle_vec[r].update();
    }
    
    // TODO - Metropolis coupling
    for (int r = 0; r < (n_rungs - 1); ++r) {
      swap_acceptance_out[r]++;
    }
    
    // store results
    lambda_store[i] = particle_vec[0].lambda;
    theta_store[i] = particle_vec[0].theta;
    decay_rate_store[i] = particle_vec[0].decay_rate;
    sens_store[i] = particle_vec[0].sens;
    n_infections_store[i] = particle_vec[0].n_infections;
    infection_times_store[i] = particle_vec[0].infection_times;
    
    iteration_counter++;
  }
  
  // store final state
  for (int r = 0; r < n_rungs; ++r) {
    writable::list tmp({
      "lambda"_nm = particle_vec[r].lambda,
      "theta"_nm = particle_vec[r].theta,
      "decay_rate"_nm = particle_vec[r].decay_rate,
      "sens"_nm = particle_vec[r].sens,
      "n_infections"_nm = particle_vec[r].n_infections,
      "infection_times"_nm = mat_double_to_list(particle_vec[r].infection_times),
      "w_list"_nm = particle_vec[r].get_w_list()
    });
    param_list_out.push_back(tmp);
  }
  
}
