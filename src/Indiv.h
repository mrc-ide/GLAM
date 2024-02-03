
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"
#include "probability.h"
#include "System.h"

//#define DEBUG_ALGO1 // uncomment to activate debugging for algorithm1

//------------------------------------------------
// class defining single individual for inference
class Indiv {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * sys;
  
  std::vector<std::vector<bool>> data;
  int n_haplos;
  std::vector<double> obs_times;
  std::vector<double> haplo_freqs;
  int n_obs;
  int n_infections;
  std::vector<double> infection_times;
  double start_time;
  double end_time;
  int max_infections;
  
  std::vector<std::vector<bool>> infection_alleles;
  std::vector<double> S_vec;
  std::vector<double> F_vec;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Indiv(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(System &sys,
            std::vector<std::vector<bool>> data_bool,
            std::vector<double> obs_times,
            int n_infections,
            std::vector<double> infection_times);
  
  void update_n_infections();
  void update_infection_times(double lambda, double theta, double decay_rate, double sens);
  double loglike_marginal_k(int k, double lambda, double theta, double decay_rate, double sens,
                            std::vector<double> &inf_times);
  void update_w_mat(double lambda, double theta, double decay_rate, double sens);
  void update_w_mat_k(int k, double lambda, double theta, double decay_rate, double sens);
  double loglike_basic(double lambda, double theta, double decay_rate, double sens);
  double algorithm1(int haplo_i, double lambda, double theta, double decay_rate, double sens,
                    std::vector<double> &inf_times, int override_k, bool override_value);
  int get_n_infections();
  std::vector<double> get_infection_times();
  
};
