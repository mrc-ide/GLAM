
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"
#include "probability.h"
#include "Indiv.h"
#include "System.h"

//------------------------------------------------
// class defining particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * sys;
  
  // parameter values
  double lambda;
  double theta;
  double decay_rate;
  double sens;
  std::vector<int> n_infections;
  int n_samp;
  std::vector<std::vector<double>> infection_times;
  int max_infections;
  
  // current loglikelihood from the forward algorithm only (conditional on other
  // objects)
  double loglike_forward;
  
  // proposal sd
  double lambda_prop_sd;
  double theta_prop_sd;
  double decay_rate_prop_sd;
  double sens_prop_sd;
  double t_inf_prop_sd;
  
  // Indiv objects
  std::vector<Indiv> indiv_vec;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(System &sys,
            double lambda,
            double theta,
            double decay_rate,
            double sens,
            std::vector<int> n_infections,
            std::vector<std::vector<double>> infection_times,
            std::vector<std::vector<std::vector<bool>>> w_array,
            std::vector<double> proposal_sd,
            double beta);
  
  void update(bool burnin, int iter);
  double get_loglike_forward(double lambda_, double theta_, double decay_rate_, double sens_);
  double get_loglike_w(double theta_);
  void update_lambda(bool burnin, int iter);
  void update_theta(bool burnin, int iter);
  void update_decay_rate(bool burnin, int iter);
  void update_sens(bool burnin, int iter);
  cpp11::list get_w_list();
  double get_logprior_lambda(double x);
  double get_logprior_theta(double x);
  double get_logprior_decay_rate(double x);
  double get_logprior_sens(double x);
  
};
