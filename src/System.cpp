
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include "System.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// initialise
void System::init(cpp11::list data_list,
                 cpp11::list obs_time_list,
                 const cpp11::doubles haplo_freqs,
                 cpp11::list param_list,
                 cpp11::list param_update_list,
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
  this->param_update_list = param_update_list;
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
  
  // get observation times in vector of vector doubles
  obs_time_vec = std::vector<std::vector<double>>(n_samp);
  for (int i = 0; i < n_samp; ++i) {
    
    doubles tmp2 = obs_time_list[i];
    obs_time_vec[i] = std::vector<double>(tmp2.size());
    for (int j = 0; j < tmp2.size(); ++j) {
      obs_time_vec[i][j] = tmp2[j];
    }
  }
  
  // get haplo_freqs into std vector
  this->haplo_freqs = std::vector<double>(haplo_freqs.size());
  for (int i = 0; i < haplo_freqs.size(); ++i) {
    this->haplo_freqs[i] = haplo_freqs[i];
  }
  
  // get which params to update
  logicals tmp = param_update_list["lambda_fixed"];
  lambda_fixed = tmp[0];
  tmp = param_update_list["theta_fixed"];
  theta_fixed = tmp[0];
  tmp = param_update_list["decay_rate_fixed"];
  decay_rate_fixed = tmp[0];
  tmp = param_update_list["sens_fixed"];
  sens_fixed = tmp[0];
  tmp = param_update_list["n_infections_fixed"];
  n_infections_fixed = tmp[0];
  tmp = param_update_list["infection_times_fixed"];
  infection_times_fixed = tmp[0];
  tmp = param_update_list["w_list_fixed"];
  w_list_fixed = tmp[0];
  
}
