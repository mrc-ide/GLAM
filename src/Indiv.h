
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"

//------------------------------------------------
// class defining single individual for inference
class Indiv {
  
public:
  // PUBLIC OBJECTS
  
  std::vector<std::vector<bool>> data;
  int n_haplos;
  std::vector<double> obs_times;
  int n_obs;
  int n_infections;
  std::vector<double> infection_times;
  
  std::vector<std::vector<bool>> infection_alleles;
  
  // RNG
  dust::random::xoshiro256plus& rng_state;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Indiv(dust::random::xoshiro256plus& rng_state) : rng_state(rng_state) {};
  
  // member functions
  void init(std::vector<std::vector<bool>> data_bool,
            std::vector<double> obs_times,
            cpp11::sexp rng_ptr,
            int n_infections,
            std::vector<double> infection_times);
  void update_n_infections();
  void update_infection_times();
  int get_n_infections();
  std::vector<double> get_infection_times();
  
};
