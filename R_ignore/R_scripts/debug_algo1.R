# debug_algo1.R
#
# Author: Bob Verity
# Date: 2024-02-02
#
# Purpose:
# Check that the C++ implementation of algorithm1 (see math notes) is correct.
# Compares outputs by *printing* from C++ code, which requires manually changing
# the Indiv.h file so the DEBUG_ALGO1 hash-define is active. Can then compare
# printed output against R calculations for the same standardised inputs.
#
# ------------------------------------------------------------------

library(tidyverse)

# ------------------------------------------------------------------
# DATA AND PARAMS

# define parameters
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 0.1
decay_rate <- 0.05
sens <- 0.95
max_infections <- 10

# simulate data
set.seed(2)
sim2 <- sim_cohort(n = 24,
                   samp_time = samp_time,
                   haplo_freqs = haplo_freqs,
                   lambda = lambda,
                   theta = theta,
                   decay_rate = decay_rate,
                   sens = sens,
                   n_inf = NULL,
                   return_full = TRUE)

# get remaining properties
n_infections <- mapply(function(x) length(x$t_inf), sim2$raw_list)
infection_times <- mapply(function(x) x$t_inf, sim2$raw_list, SIMPLIFY = FALSE)
w_list <- mapply(function(x) x$w_inf, sim2$raw_list, SIMPLIFY = FALSE)

# -------------------------------
# RUN C++ IMPLEMENTATION

# run C++ version
g <- glam_mcmc$new(df_data = sim2$df_data)

g$init(start_time = 0,
       end_time = 10,
       lambda = lambda,
       theta = theta,
       decay_rate = decay_rate,
       sens = sens,
       n_infections = n_infections,
       infection_times = infection_times,
       w_list = w_list,
       haplo_freqs = haplo_freqs)

g$debug_algo1()

# -------------------------------
# RUN R IMPLEMENTATION

algo1 <- function(haplo_i, lambda, theta, decay_rate, sens, inf_times,
                  override_k, override_value, haplo_freqs, dat_bool,
                  samp_time, w_mat) {
  
  p <- haplo_freqs[haplo_i]
  q <- 1.0 - exp(-theta*p)  # chance of haplo being introduced
  prob_equilib <- lambda*q / (lambda*q + decay_rate)
  prob_given_pos <- ifelse(dat_bool[1] == 1, sens, 1 - sens)
  prob_given_neg <- ifelse(dat_bool[1] == 1, 0, 1)
  
  # create single data.frame of all observations and infections in increasing order
  df_all <- data.frame(time = samp_time,
                             type = "obs",
                             value = dat_bool)
  if (!is.null(w_mat)) {
    df_all <- rbind(df_all,
                    data.frame(time = inf_times,
                               type = "inf",
                               value = w_mat[,haplo_i]))
  }
  df_all <- df_all |>
    arrange(time) |>
    mutate(A = NA, B = NA)
  
  # initialise
  df_all$A[1] <- prob_equilib * prob_given_pos
  df_all$B[1] <- (1 - prob_equilib) * prob_given_neg
  
  # loop through steps
  for (i in 2:nrow(df_all)) {
    
    T11 <- exp(-decay_rate*(df_all$time[i] - df_all$time[i - 1]))
    T10 <- 1.0 - T11
    
    if (df_all$type[i] == "obs") {
      df_all$A[i] <- df_all$A[i - 1] * T11 * ifelse(df_all$value[i] == 1, sens, 1 - sens)
      df_all$B[i] <- (df_all$A[i - 1] * T10 + df_all$B[i - 1]) * ifelse(df_all$value[i] == 1, 0, 1)
    } else {
      df_all$A[i] <- ifelse(df_all$value[i] == 1, df_all$A[i - 1] + df_all$B[i - 1], df_all$A[i - 1] * T11)
      df_all$B[i] <- ifelse(df_all$value[i] == 1, 0, df_all$A[i - 1] * T10 + df_all$B[i - 1])
    }
  }
  
  return(df_all)
}

# get data as list
data_list <- g$.__enclos_env__$private$list_data

# run for one combination
indiv_i <- 1
haplo_i <- 1
dat_bool <- data_list[[indiv_i]][,haplo_i][[1]]
w_mat <- w_list[[indiv_i]]

df_algo1 <- algo1(haplo_i = 1,
                  lambda = lambda,
                  theta = theta,
                  decay_rate = decay_rate,
                  sens = sens,
                  inf_times = infection_times[[indiv_i]],
                  override_k = -1,
                  override_value = FALSE,
                  haplo_freqs = haplo_freqs,
                  dat_bool = dat_bool,
                  samp_time = samp_time,
                  w_mat = w_mat)

df_algo1

