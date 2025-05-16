
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
                    std::vector<std::vector<std::vector<bool>>> w_array,
                    std::vector<double> proposal_sd,
                    double beta) {
  
  // pointer to system object
  this->sys = &sys;
  
  // initialise RNG
  rng_state = sys.rng_state;
  
  // copy over known values
  this->lambda = lambda;
  this->theta = theta;
  this->decay_rate = decay_rate;
  this->sens = sens;
  this->n_infections = n_infections;
  this->infection_times = infection_times;
  this->max_infections = sys.max_infections;
  n_samp = n_infections.size();
  
  loglike_forward = 0.0;
  
  // proposal standard deviations
  lambda_prop_sd = proposal_sd[0];
  theta_prop_sd = proposal_sd[1];
  decay_rate_prop_sd = proposal_sd[2];
  sens_prop_sd = proposal_sd[3];
  t_inf_prop_sd = proposal_sd[4];
  
  // initialise Indiv objects
  indiv_vec = std::vector<Indiv>(n_samp, Indiv(rng_state));
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].init(sys,
                      sys.data_bool[i],
                      sys.obs_time_vec[i],
                      sys.obs_time_start[i],
                      sys.obs_time_end[i],
                      n_infections[i],
                      infection_times[i],
                      w_array[i]);
  }
  
#ifdef DEBUG_ALGO1
  
  // set combination to debug
  int indiv_i = 0;
  int haplo_i = 0;
  
  // print out parameter values
  print("PARAMETERS:");
  print("lambda:", lambda);
  print("theta:", theta);
  print("decay_rate:", decay_rate);
  print("sens:", sens);
  print("indiv:", indiv_i);
  print("haplo:", haplo_i);
  print("infection_times:");
  print_vector(infection_times[indiv_i]);
  print("infection_alleles:");
  print_matrix(indiv_vec[indiv_i].infection_alleles);
  
  // run algorithm and print running output
  print("\nALGORITHM:");
  double tmp = indiv_vec[indiv_i].algorithm1(haplo_i, lambda, theta, decay_rate, sens,
                                             infection_times[indiv_i], -1, false);
  
  // print final value
  print("result:", tmp);
  stop("end debugging");
#endif
  
}

//------------------------------------------------
// update
void Particle::update(bool burnin, int iter) {
  
  //std::this_thread::sleep_for(std::chrono::milliseconds(1));
  
  // split-merge update steps on all individuals
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].update_n_infections(lambda, theta, decay_rate, sens);
  }
  
  // update all infection times
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].update_infection_times(lambda, theta, decay_rate, sens, burnin, t_inf_prop_sd, iter);
  }
  
  // Gibbs sample W matrix
  for (int i = 0; i < n_samp; ++i) {
    indiv_vec[i].update_w_mat(lambda, theta, decay_rate, sens);
  }
  
  // store all infection times
  for (int i = 0; i < n_samp; ++i) {
    n_infections[i] = indiv_vec[i].get_n_infections();
    infection_times[i] = indiv_vec[i].get_infection_times();
  }
  
  // calculate basic log-likelihood in current state. This is from the forward
  // algorithm only, and is still conditional on many parameters/objects
  loglike_forward = get_loglike_forward(lambda, theta, decay_rate, sens);
  
  // update global parameters
  update_lambda(burnin, iter);
  update_theta(burnin, iter);
  update_decay_rate(burnin, iter);
  update_sens(burnin, iter);
  
}

//------------------------------------------------
// get basic loglikelihood from forward algorithm over all Indiv
double Particle::get_loglike_forward(double lambda_, double theta_, double decay_rate_, double sens_) {
  
  double ret = 0.0;
  for (int i = 0; i < n_samp; ++i) {
    ret += indiv_vec[i].get_loglike_forward(lambda_, theta_, decay_rate_, sens_);
  }
  return ret;
}

//------------------------------------------------
// get log-probability of w matrix over all Indiv
double Particle::get_loglike_w(double theta_) {
  
  int n_haplos = indiv_vec[0].n_haplos; // this is a bit sloppy
  std::vector<double> log_prob_pos(n_haplos);
  std::vector<double> log_prob_neg(n_haplos);
  for (int i = 0; i < n_haplos; ++i) {
    log_prob_pos[i] = log(1.0 - exp(-theta_*sys->haplo_freqs[i]));
    log_prob_neg[i] = -theta_*sys->haplo_freqs[i];
  }
  double log_norm = -log(1.0 - exp(-theta_));
  
  double ret = 0.0;
  for (int i = 0; i < n_samp; ++i) {
    ret += indiv_vec[i].get_loglike_w(log_prob_pos, log_prob_neg, log_norm);
  }
  return ret;
}

//------------------------------------------------
// update
void Particle::update_lambda(bool burnin, int iter) {
  if (sys->lambda_fixed) {
    return;
  }
  
  // propose move
  double lambda_prop = rnorm1_pos(rng_state, lambda, lambda_prop_sd);
  
  // calculate new log-likelihood
  double loglike_forward_prop = get_loglike_forward(lambda_prop, theta, decay_rate, sens);
  
  // account for other terms required to convert conditional probability into joint probability
  double adj_current = 0;
  double adj_prop = 0;
  for (int i = 0; i < n_samp; ++i) {
    int ni = indiv_vec[i].n_infections;
    double delta_t = sys->obs_time_end[i] - sys->obs_time_start[i];
    adj_current += ni*log(lambda) - lambda*delta_t;
    adj_prop += ni*log(lambda_prop) - lambda_prop*delta_t;
  }
  
  // get log-priors
  double logprior_current = get_logprior_lambda(lambda);
  double logprior_prop = get_logprior_lambda(lambda_prop);
  
  // Metropolis-Hastings step
  double MH = (loglike_forward_prop + logprior_prop + adj_prop) - (loglike_forward + logprior_current + adj_current);
  
  if (log(runif1(rng_state)) < MH) {
    lambda = lambda_prop;
    loglike_forward = loglike_forward_prop;
    
    // Robbins Monroe step
    if (burnin) {
      lambda_prop_sd = exp(log(lambda_prop_sd) + (1 - 0.23) / sqrt(iter));
    }
    
  } else {
    // Robbins Monroe step
    if (burnin) {
      lambda_prop_sd = exp(log(lambda_prop_sd) - 0.23 / sqrt(iter));
    }
  }
  
}

//------------------------------------------------
// update
void Particle::update_theta(bool burnin, int iter) {
  if (sys->theta_fixed) {
    return;
  }
  
  // propose move
  double theta_prop = rnorm1_pos(rng_state, theta, theta_prop_sd);
  
  // calculate new log-likelihood
  double loglike_forward_prop = get_loglike_forward(lambda, theta_prop, decay_rate, sens);
  
  // account for other terms required to convert conditional probability into joint probability
  double adj_current = get_loglike_w(theta);
  double adj_prop = get_loglike_w(theta_prop);
  
  // get log-priors
  double logprior_current = get_logprior_theta(theta);
  double logprior_prop = get_logprior_theta(theta_prop);
  
  // Metropolis-Hastings step
  double MH = (loglike_forward_prop + logprior_prop + adj_prop) - (loglike_forward + logprior_current + adj_current);
  
  if (log(runif1(rng_state)) < MH) {
    theta = theta_prop;
    loglike_forward = loglike_forward_prop;
    
    // Robbins Monroe step
    if (burnin) {
      theta_prop_sd = exp(log(theta_prop_sd) + (1 - 0.23) / sqrt(iter));
    }
  } else {
    // Robbins Monroe step
    if (burnin) {
      theta_prop_sd = exp(log(theta_prop_sd) - 0.23 / sqrt(iter));
    }
  }
}

//------------------------------------------------
// update
void Particle::update_decay_rate(bool burnin, int iter) {
  if (sys->decay_rate_fixed) {
    return;
  }
  
  // propose move
  double decay_rate_prop = rnorm1_pos(rng_state, decay_rate, decay_rate_prop_sd);
  
  // calculate new log-likelihood
  double loglike_forward_prop = get_loglike_forward(lambda, theta, decay_rate_prop, sens);
  
  // get log-priors
  double logprior_current = get_logprior_decay_rate(decay_rate);
  double logprior_prop = get_logprior_decay_rate(decay_rate_prop);
  
  // Metropolis-Hastings step
  double MH = (loglike_forward_prop + logprior_prop ) - (loglike_forward + logprior_current);
  
  if (log(runif1(rng_state)) < MH) {
    decay_rate = decay_rate_prop;
    loglike_forward = loglike_forward_prop;
    
    // Robbins Monroe step
    if (burnin) {
      decay_rate_prop_sd = exp(log(decay_rate_prop_sd) + (1 - 0.23) / sqrt(iter));
    }
  } else {
    // Robbins Monroe step
    if (burnin) {
      decay_rate_prop_sd = exp(log(decay_rate_prop_sd) - 0.23 / sqrt(iter));
    }
  }
}
//------------------------------------------------
// update
void Particle::update_sens(bool burnin, int iter) {
  if (sys->sens_fixed) {
    return;
  }
  
  // propose move
  double sens_prop = rnorm1_interval(rng_state, sens, sens_prop_sd, 0, 1);
  
  // calculate new log-likelihood
  double loglike_forward_prop = get_loglike_forward(lambda, theta, decay_rate, sens_prop);
  
  // get log-priors
  double logprior_current = get_logprior_sens(sens);
  double logprior_prop = get_logprior_sens(sens_prop);
  
  // Metropolis-Hastings step
  double MH = (loglike_forward_prop + logprior_prop ) - (loglike_forward + logprior_current);
  
  if (log(runif1(rng_state)) < MH) {
    sens = sens_prop;
    loglike_forward = loglike_forward_prop;
    
    // Robbins Monroe step
    if (burnin) {
      sens_prop_sd = exp(log(sens_prop_sd) + (1 - 0.23) / sqrt(iter));
    }
  } else {
    // Robbins Monroe step
    if (burnin) {
      sens_prop_sd = exp(log(sens_prop_sd) - 0.23 / sqrt(iter));
    }
  }
}

//------------------------------------------------
// return infection alleles over all Indiv as a nested list
list Particle::get_w_list() {
  
  cpp11::writable::list ret;
  for (int i = 0; i < n_samp; ++i) {
    cpp11::writable::list tmp;
    for (int j = 0; j < indiv_vec[i].n_infections; ++j) {
      cpp11::logicals v = cpp11::as_sexp(indiv_vec[i].infection_alleles[j]);
      tmp.push_back(v);
    }
    ret.push_back(tmp);
  }
  
  return ret;
}

//------------------------------------------------
// log-prior on lambda
double Particle::get_logprior_lambda(double x) {
  SEXP ret_raw = sys->lambda_prior(x);
  cpp11::doubles ret_cpp11(ret_raw);
  double ret = ret_cpp11[0];
  return ret;
}

//------------------------------------------------
// log-prior on theta
double Particle::get_logprior_theta(double x) {
  SEXP ret_raw = sys->theta_prior(x);
  cpp11::doubles ret_cpp11(ret_raw);
  double ret = ret_cpp11[0];
  return ret;
}

//------------------------------------------------
// log-prior on decay rate
double Particle::get_logprior_decay_rate(double x) {
  SEXP ret_raw = sys->decay_rate_prior(x);
  cpp11::doubles ret_cpp11(ret_raw);
  double ret = ret_cpp11[0];
  return ret;
}

//------------------------------------------------
// log-prior on sensitivity
double Particle::get_logprior_sens(double x) {
  SEXP ret_raw = sys->sens_prior(x);
  cpp11::doubles ret_cpp11(ret_raw);
  double ret = ret_cpp11[0];
  return ret;
}

