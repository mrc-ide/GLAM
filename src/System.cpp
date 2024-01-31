
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "System.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// initialise
void System::init(cpp11::list data_list,
                 cpp11::list obs_time_list,
                 cpp11::list param_list,
                 cpp11::list proposal_sd,
                 int iteration_counter_init,
                 cpp11::doubles beta,
                 double start_time,
                 double end_time,
                 int max_infections,
                 cpp11::sexp rng_ptr) {
  
  // copy over values
  this->data_list = data_list;
  this->obs_time_list = obs_time_list;
  this->param_list = param_list;
  this->proposal_sd = proposal_sd;
  this->iteration_counter_init = iteration_counter_init;
  this->beta = beta;
  this->start_time = start_time;
  this->end_time = end_time;
  this->max_infections = max_infections;
  
  n_samp = data_list.size();
  n_rungs = beta.size();
  
  // initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  rng_state = rng->state(0);
  
  // convert data to vector of boolean matrices
  data_bool = std::vector<std::vector<std::vector<bool>>>(n_samp);
  for (int i = 0; i < n_samp; ++i) {
    
    // convert data to boolean matrix for this individual
    list tmp = data_list[i];
    data_bool[i] = std::vector<std::vector<bool>>(tmp.size());
    for (int j = 0; j < tmp.size(); ++j) {
      doubles tmp_j = tmp[j];
      data_bool[i][j] = std::vector<bool>(tmp_j.size());
      for (int k = 0; k < tmp_j.size(); ++k) {
        data_bool[i][j][k] = tmp_j[k];
      }
    }
  }
  
}
