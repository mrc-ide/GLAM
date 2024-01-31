
#include <chrono>
#include <thread>


#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "Particle.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// constructor
void Particle::init(System &sys,
                    double lambda,
                    double theta,
                    double decay_rate,
                    double sens,
                    std::vector<int> n_infections,
                    std::vector<std::vector<double>> infection_times,
                    std::vector<double> proposal_sd,
                    double beta,
                    double start_time,
                    double end_time,
                    int max_infections) {
  
  // initialise RNG
  rng_state = sys.rng_state;
  
  // copy over known values
  this->lambda = lambda;
  this->theta = theta;
  this->decay_rate = decay_rate;
  this->sens = sens;
  this->n_infections = n_infections;
  this->infection_times = infection_times;
  this->start_time = start_time;
  this->end_time = end_time;
  this->max_infections = max_infections;
  n_samp = n_infections.size();
  
  // proposal_sd_vec goes through {lambda, theta, decay_rate, sens, infection_time}
  proposal_sd_vec = proposal_sd;
  n_proposal_sd = proposal_sd_vec.size();
  
  // initialise Indiv objects
  indiv_vec = std::vector<Indiv>(n_samp, Indiv(rng_state));
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].init(sys,
                      sys.data_bool[i],
                      sys.obs_time_vec[i],
                      start_time,
                      end_time,
                      max_infections,
                      n_infections[i],
                      infection_times[i]);
  }
  
}

//------------------------------------------------
// update
void Particle::update() {
  
  // split-merge update steps on all individuals
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].update_n_infections();
  }
  
  // update all infection times
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].update_infection_times();
  }
  
  // store all infection times
  for (int i = 0; i < n_samp; ++i) {
    n_infections[i] = indiv_vec[i].get_n_infections();
    infection_times[i] = indiv_vec[i].get_infection_times();
  }
  
  // update global parameters
  update_lambda();
  update_theta();
  update_decay_rate();
  update_sens();
  
}

//------------------------------------------------
// update
void Particle::update_lambda() {
  
  // dummy update
  double lambda_prop = rnorm1_pos(rng_state, lambda, proposal_sd_vec[0]);
  
  double loglike = 0.0;
  double loglike_prop = 0.0;
  for (int i = 0; i < n_samp; ++i) {
    loglike += indiv_vec[i].loglike_basic(lambda, theta, decay_rate, sens);
    loglike_prop += indiv_vec[i].loglike_basic(lambda_prop, theta, decay_rate, sens);
  }
  
  double MH = loglike_prop - loglike;
  
  if (log(runif1(rng_state)) < MH) {
    lambda = lambda_prop;
  }
  
}

//------------------------------------------------
// update
void Particle::update_theta() {
  theta = dust::random::random_real<double>(rng_state);
}

//------------------------------------------------
// update
void Particle::update_decay_rate() {
  decay_rate = dust::random::random_real<double>(rng_state);
}
//------------------------------------------------
// update
void Particle::update_sens() {
  sens = dust::random::random_real<double>(rng_state);
}
