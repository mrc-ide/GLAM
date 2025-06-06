
#include <chrono>
#include <thread>
#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/normal.hpp>
#include "Indiv.h"

using namespace cpp11;
namespace writable = cpp11::writable;

//------------------------------------------------
// initialise
void Indiv::init(System &sys,
                 std::vector<std::vector<bool>> data_bool,
                 std::vector<double> obs_times,
                 double start_time,
                 double end_time,
                 int n_infections,
                 std::vector<double> infection_times,
                 std::vector<std::vector<bool>> infection_alleles) {
  
  // pointer to system object
  this->sys = &sys;
  
  // initialise RNG
  rng_state = sys.rng_state;
  
  // copy over values
  data = data_bool;
  n_haplos = data.size();
  this->obs_times = obs_times;
  this->start_time = start_time;
  this->end_time = end_time;
  haplo_freqs = sys.haplo_freqs;
  n_obs = data[0].size();
  this->n_infections = n_infections;
  this->infection_times = infection_times;
  this->infection_alleles = infection_alleles;
  this->max_infections = sys.max_infections;
  
  // sanity checks on inputs
  if (infection_times.size() != n_infections) {
    stop("Error in Indiv: infection_times.size() does not match n_infections");
  }
  if (obs_times.size() != n_obs) {
    stop("Error in Indiv: obs_times.size() does not match n_obs");
  }
  
  // vectors for storing success and failure from algorithm1
  S_vec = std::vector<double>(n_haplos);
  F_vec = std::vector<double>(n_haplos);
}

//------------------------------------------------
// update number of infections
void Indiv::update_n_infections(double lambda, double theta, double decay_rate, double sens) {
  if (sys->n_infections_fixed) {
    return;
  }
  
  // dummy update step. Randomly drop or add infection
  if (runif1(rng_state) < 0.5) {  // add infection
    if (n_infections == max_infections) {
      return;
    }
    
    // get log-likelihood of current state
    double loglike_now = get_loglike_forward(lambda, theta, decay_rate, sens);
    
    // generate a new value over the entire interval
    double new_time = runif1(rng_state, start_time, end_time);
    
    // find out where it belongs in the infection_times vector (iterator form)
    // so we can insert while maintaining sorting in increasing value
    auto it_time = std::lower_bound(infection_times.begin(), infection_times.end(), new_time);
    
    // convert to int
    int idx = static_cast<int>(std::distance(infection_times.begin(), it_time));
    
    // insert into infection_times
    infection_times.insert(it_time, new_time);
    
    // insert a row of dummy values into infection_alleles
    std::vector<bool> new_alleles(n_haplos, true);
    infection_alleles.insert(infection_alleles.begin() + idx, std::move(new_alleles));
    
    n_infections++;
    
    // get log-likelihood of proposed state, marginalized over idx
    double loglike_prop = get_loglike_marginal_k(idx, lambda, theta, decay_rate, sens, infection_times);
    
    // get log ratio of proposed vs current prior
    double prior_ratio = log(lambda);
    
    // get log ratio of current vs. proposed
    double propose_ratio = log(end_time - start_time) - log(n_infections);
    
    // calculate Metropolis-Hastings ratio
    double MH = (loglike_prop - loglike_now) + prior_ratio + propose_ratio;
    bool accept_move = (log(runif1(rng_state)) < MH);
    
    // accept or reject move
    if (accept_move) {
      
      // Gibbs sample introduced haplos given new timings
      update_w_mat_k(idx, lambda, theta, decay_rate, sens);
    } else {
      
      // revert back to previous state
      infection_times.erase(infection_times.begin() + idx);
      infection_alleles.erase(infection_alleles.begin() + idx);
      n_infections--;
    }
    
  } else {  // drop infection
    if (n_infections == 0) {
      return;
    }
    
    // draw which value to drop with equal probability
    int idx = sample2(rng_state, 0, n_infections - 1);
    
    // get log-likelihood of current state, marginalized over idx
    double loglike_now = get_loglike_marginal_k(idx, lambda, theta, decay_rate, sens, infection_times);
    
    // store the alleles at this row in case we need to reinstate them
    double store_time = infection_times[idx];
    std::vector<bool> store_alleles = infection_alleles[idx];
    
    // drop values
    infection_times.erase(infection_times.begin() + idx);
    infection_alleles.erase(infection_alleles.begin() + idx);
    n_infections--;
    
    // get log-likelihood of proposed state
    double loglike_prop = get_loglike_forward(lambda, theta, decay_rate, sens);
    
    // get log ratio of proposed vs current prior
    double prior_ratio = -log(lambda);
    
    // get log ratio of current vs. proposed
    double propose_ratio = log(n_infections + 1) - log(end_time - start_time);
    
    // calculate Metropolis-Hastings ratio
    double MH = (loglike_prop - loglike_now) + prior_ratio + propose_ratio;
    bool accept_move = (log(runif1(rng_state)) < MH);
    
    // accept or reject move
    if (!accept_move) {
      
      // revert back to previous state
      infection_times.insert(infection_times.begin() + idx, store_time);
      infection_alleles.insert(infection_alleles.begin() + idx, store_alleles);
      n_infections++;
    }
    
  }
  
}

//------------------------------------------------
// update timings of all infections
void Indiv::update_infection_times(double lambda, double theta, double decay_rate, double sens,
                                   bool burnin, double &t_inf_prop_sd, int iter) {
  
  if (sys->infection_times_fixed) {
    return;
  }
  
  std::vector<double> infection_times_prop = infection_times;
  for (int k = 0; k < n_infections; ++k) {
    
    // propose new value within interval
    double bound_lower = start_time;
    double bound_upper = end_time;
    if (k > 0) {
      bound_lower = infection_times[k - 1];
    }
    if (k < (n_infections - 1)) {
      bound_upper = infection_times[k + 1];
    }
    infection_times_prop[k] = rnorm1_interval(rng_state, infection_times[k], t_inf_prop_sd, bound_lower, bound_upper);
    
    // calculate loglikelihood of current and proposed state
    double loglike_now = get_loglike_marginal_k(k, lambda, theta, decay_rate, sens, infection_times);
    double loglike_prop = get_loglike_marginal_k(k, lambda, theta, decay_rate, sens, infection_times_prop);
    
    // calculate Metropolis-Hastings ratio
    double MH = loglike_prop - loglike_now;
    bool accept_move = (log(runif1(rng_state)) < MH);
    
    // accept or reject move
    if (accept_move) {
      infection_times[k] = infection_times_prop[k];
      
      // Gibbs sample introduced haplos given new timings
      if (!sys->w_list_fixed) {
        update_w_mat_k(k, lambda, theta, decay_rate, sens);
      }
      
      // Robbins Monroe step
      if (burnin && (t_inf_prop_sd < (end_time - start_time))) {
        t_inf_prop_sd = exp(log(t_inf_prop_sd) + (1 - 0.23) / iter);
      }
    } else {
      infection_times_prop[k] = infection_times[k];
      
      // Robbins Monroe step
      if (burnin) {
        t_inf_prop_sd = exp(log(t_inf_prop_sd) - 0.23 / iter);
      }
    }
    
  }
  
}

//------------------------------------------------
// log-likelihood marginalised over the kth infection
double Indiv::get_loglike_marginal_k(int k, double lambda, double theta, double decay_rate, double sens,
                                     std::vector<double> &inf_times) {
  
  // populate vectors of (logged) success and failure probabilities
  for (int j = 0; j < n_haplos; ++j) {
    S_vec[j] = algorithm1(j, lambda, theta, decay_rate, sens, inf_times, k, true);
    F_vec[j] = algorithm1(j, lambda, theta, decay_rate, sens, inf_times, k, false);
  }
  
  double l1 = 0;
  double l2 = -theta;
  for (int j = 0; j < n_haplos; ++j) {
    double q = 1.0 - exp(-theta * haplo_freqs[j]);
    l1 += log(q) + S_vec[j] + log(1.0 + (1.0/q - 1.0)*exp(F_vec[j] - S_vec[j]));
    l2 += F_vec[j];
  }
  double ret = l1 + log(1.0 - exp(l2 - l1));
  ret -= log(1.0 - exp(-theta));
  
  return ret;
}

//------------------------------------------------
// Gibbs sampler on W matrix
void Indiv::update_w_mat(double lambda, double theta, double decay_rate, double sens) {
  if (sys->w_list_fixed) {
    return;
  }
  
  for (int k = 0; k < n_infections; ++k) {
    update_w_mat_k(k, lambda, theta, decay_rate, sens);
  }
}

//------------------------------------------------
// Gibbs sampler on kth infection of W matrix
void Indiv::update_w_mat_k(int k, double lambda, double theta, double decay_rate, double sens) {
  
  // populate vectors of (logged) success and failure probabilities
  for (int j = 0; j < n_haplos; ++j) {
    S_vec[j] = algorithm1(j, lambda, theta, decay_rate, sens, infection_times, k, true);
    F_vec[j] = algorithm1(j, lambda, theta, decay_rate, sens, infection_times, k, false);
  }
  
  // update each haplo in turn
  bool no_success = true;
  for (int j = 0; j < n_haplos; ++j) {
    
    // calculate basic probability
    double q = 1.0 - exp(-theta * haplo_freqs[j]);
    double prob_success = q / (q + (1 - q)*exp(F_vec[j] - S_vec[j]));
    
    // modify probability if no success seen so far
    if (no_success) {
      double log_a = 0.0;
      double log_b = 0.0;
      for (int l = j; l < n_haplos; ++l) {
        double ql = 1.0 - exp(-theta * haplo_freqs[l]);
        log_a += log(ql) + S_vec[l] + log(1.0 + (1.0/q - 1)*exp(F_vec[l] - S_vec[j]));
        log_b += log(1 - ql) + F_vec[l];
      }
      prob_success /= (1.0 - exp(log_b - log_a));
    }
    
    // draw new value of whether haplo introduced
    infection_alleles[k][j] = rbernoulli1(rng_state, prob_success);
    if (infection_alleles[k][j]) {
      no_success = false;
    }
  }
}

//------------------------------------------------
// log-probability of w matrix (infection_alleles) given vector of (log)
// probabilities of being positive or negative for each haplotype, plus a log
// normalisation constant that is applied to each infection to account for
// zero-truncation
double Indiv::get_loglike_w(std::vector<double> &log_prob_pos, std::vector<double> &log_prob_neg, double log_norm) {
  
  double ret = 0.0;
  for (int i = 0; i < n_infections; ++i) {
    ret += log_norm;
    for (int j = 0; j < n_haplos; ++j) {
      if (infection_alleles[i][j]) {
        ret += log_prob_pos[j];
      } else {
        ret += log_prob_neg[j];
      }
    }
  }
  return ret;
}

//------------------------------------------------
// basic log-likelihood from forward algorithm, no marginalisation
double Indiv::get_loglike_forward(double lambda, double theta, double decay_rate, double sens) {
  
  double ret = 0.0;
  for (int j = 0; j < n_haplos; ++j) {
    ret += algorithm1(j, lambda, theta, decay_rate, sens, infection_times, -1, true);
  }
  return ret;
}

//------------------------------------------------
// algorithm1 (see math notes)
double Indiv::algorithm1(int haplo_i, double lambda, double theta, double decay_rate, double sens,
                         std::vector<double> &inf_times, int override_k, bool override_value) {
  
  // extract key parameters and calculate starting values of A and B based on
  // equilibrium arguments
  int n_inf = inf_times.size();
  double p = sys->haplo_freqs[haplo_i];
  double q = 1.0 - exp(-theta*p);
  double prob_equilib = lambda*q / (lambda*q + decay_rate);
  double prob_given_pos = data[haplo_i][0]*sens + (1 - data[haplo_i][0])*(1 - sens);
  double prob_given_neg = 1 - data[haplo_i][0];
  double A = prob_equilib * prob_given_pos;
  double B = (1 - prob_equilib) * prob_given_neg;
  double log_running = 0;
  
#ifdef DEBUG_ALGO1
  std::stringstream ss;
  ss << "init:        " << A << " " << B;
  cpp11::message(ss.str());
#endif
  
  int i = 1; // indexes sampling times, starting from second
  int c = 0; // indexes infections
  double tau = start_time;
  bool next_is_infection = false;
  while (i < obs_times.size()) {
    
    // work out if next event is an infection or an observation
    if ((n_inf == 0) || (c == n_inf)) {
      next_is_infection = false;
    } else {
      if (obs_times[i] < inf_times[c]) {
        next_is_infection = false;
      } else {
        next_is_infection = true;
      }
    }
    
    // implement changes
    if (next_is_infection) {
      
      bool w = infection_alleles[c][haplo_i];
      if (c == override_k) {
        w = override_value;
      }
      double T11 = exp(-decay_rate * (inf_times[c] - tau));
      double A_new = w*(A + B) + (1 - w)*A*T11;
      double B_new = (1 - w)*(A*(1 - T11) + B);
      log_running += log(A_new + B_new);
      A = A_new / (A_new + B_new);
      B = B_new / (A_new + B_new);
      tau = inf_times[c];
      c++;
      
    } else {
      
      bool data_i = data[haplo_i][i];
      double T11 = exp(-decay_rate * (obs_times[i] - tau));
      double prob_given_pos = data_i*sens + (1 - data_i)*(1 - sens);
      double prob_given_neg = 1 - data_i;
      double A_new = A*T11*prob_given_pos;
      double B_new = (A*(1 - T11) + B)*prob_given_neg;
      log_running += log(A_new + B_new);
      A = A_new / (A_new + B_new);
      B = B_new / (A_new + B_new);
      tau = obs_times[i];
      i++;
      
    }
    
#ifdef DEBUG_ALGO1
    std::stringstream ss;
    if (next_is_infection) {
      ss << "infection:   ";
    } else {
      ss << "observation: ";
    }
    ss << A << " " << B;
    cpp11::message(ss.str());
#endif
  }
  
  return log_running;
}

//------------------------------------------------
// return n_infections
int Indiv::get_n_infections() {
  return n_infections;
}

//------------------------------------------------
// return vector of infection times
std::vector<double> Indiv::get_infection_times() {
  return infection_times;
}
