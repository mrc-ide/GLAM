
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "misc.h"

//------------------------------------------------
// class defining particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // parameter values
  double lambda;
  double theta;
  
  // proposal sd
  std::vector<double> proposal_sd_vec;
  int n_proposal_sd;
  
  // RNG
  cpp11::sexp rng_ptr;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // member functions
  void init(double lambda,
            double theta,
            std::vector<double> proposal_sd,
            double beta,
            cpp11::sexp rng_ptr);
  void update();
  
};
