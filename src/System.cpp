
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
  
  n_rungs = beta.size();
  
  // proposal sd
  proposal_sd_mat = list_to_mat_double(proposal_sd);
  n_proposal_sd = proposal_sd_mat[0].size();
  
  
  // initialise RNG
  auto rng = dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  rng_state = rng->state(0);
  
}

//------------------------------------------------
// get one row of the proposal_sd_mat matrix
std::vector<double> System::get_proposal_sd_vec(int rung_index) {
  return proposal_sd_mat[rung_index];
}
/*
//------------------------------------------------
// get observation data for this individual
std::vector<std::vector<bool>> System::get_data_bool(int ind_index) {
  
  // convert data to boolean matrix for this individual
  list tmp = data_list[ind_index];
  std::vector<std::vector<bool>> data_bool(tmp.size());
  for (int i = 0; i < tmp.size(); ++i) {
    doubles tmp_i = tmp[i];
    data_bool[i] = std::vector<bool>(tmp_i.size());
    for (int j = 0; j < tmp_i.size(); ++j) {
      data_bool[i][j] = tmp_i[j];
    }
  }
  
  return data_bool;
}
*/
