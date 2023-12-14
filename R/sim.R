
#------------------------------------------------
#' @title Simulate a single individual
#'
#' @description Simulate a single individual from the generative model. Returns
#'   outputs as a list, including observed data as well as optionally hidden
#'   true states. If hidden states are returned then \code{plot_ind()} can be
#'   used to visualise this simulation.
#'
#' @details Note that the method for simulating data is deliberately different
#'   from the way that it is represented in the inference step, although it is
#'   mathematically equivalent. This allows us to check that our different
#'   representations of the process are consistent.
#'
#' @param samp_time numeric vector of times at which observations were taken.
#'   Must be at least two unique time points increasing in value.
#' @param haplo_freqs vector of population-level haplotype frequencies.
#' @param lambda force of infection.
#' @param theta COI intensity parameter.
#' @param decay_rate rate at which haplotypes clear.
#' @param sens sensitivity of sequencing. Assumed the same for all haplotypes.
#' @param ind_name the name given to this individual in the returned data.frame.
#' @param return_full Boolean (Default \code{FALSE}). If \code{TRUE} then
#'   return true values of hidden variables alongside observed data.
#'
#' @importFrom stats rbinom rexp rpois runif
#' @export

sim_ind <- function(samp_time, haplo_freqs, lambda, theta, decay_rate, sens,
                    ind_name = "ind1", return_full = FALSE) {
  
  ret <- cpp_test()
  return(ret)
  
  # check inputs
  assert_vector_numeric(samp_time)
  assert_greq(length(samp_time), 2)
  assert_increasing(samp_time)
  assert_noduplicates(samp_time)
  assert_vector_bounded(haplo_freqs)
  assert_eq(sum(haplo_freqs), 1)
  assert_single_pos(lambda)
  assert_single_pos(theta, zero_allowed = FALSE)
  assert_single_pos(decay_rate)
  assert_single_bounded(sens)
  assert_single_string(ind_name)
  assert_single_logical(return_full)
  
  # get dimensions
  n_haplo <- length(haplo_freqs)
  n_samp <- length(samp_time)
  start_time <- samp_time[1]
  end_time <- samp_time[n_samp]
  t_period <- end_time - start_time
  
  # transform raw haplotype frequencies into the probability of transmission
  q <- 1 - exp(-theta*haplo_freqs)
  
  #------------------------------------------------
  # draw the true timings of all haplo introductions and clearance. At this
  # stage we do not worry about infections being "overwritten" by later
  # infections. In the next stage we will truncate all clearance times as needed
  # to avoid this.
  
  # draw the number of infections that occur during the observation period and
  # the timings of these infections
  n_inf <- rpois(1, lambda*t_period)
  t_inf <- sort(runif(n_inf, start_time, end_time))
  
  # consider whether this individual initialises positive for each haplotype.
  # The probability of initialising positive is assumed to be given by the
  # relative rates of acquiring haplotypes vs. losing them.
  p_equilib <- lambda*q / (lambda*q + decay_rate)
  w_init <- rbinom(n_haplo, 1, p_equilib)
  
  # draw when initial haplotypes clear
  clear_init <- rexp(n_haplo, rate = decay_rate)
  clear_init[w_init == 0] <- NA
  
  # w_inf and clear_inf are matrices storing whether an allele is introduced,
  # and it's clearance time. If there are no infections then these are set to
  # NULL
  w_inf <- clear_inf <- NULL
  if (n_inf > 0) {
    w_inf <- matrix(0, n_inf, n_haplo)
    clear_inf <- matrix(0, n_inf, n_haplo)
    
    # loop through all new infections
    for (i in 1:n_inf) {
      
      # draw the number of infected hepatocytes from zero-truncated Poisson
      hepat <- rztpois(1, theta)
      
      # draw the haplos in each infection
      w_inf[i,] <- rmultinom(1, hepat, prob = q)[,1] > 0
      
      # draw the clearance times of haplos
      clear_inf[i,] <- t_inf[i] + rexp(n_haplo, rate = decay_rate)
    }
    clear_inf[w_inf == 0] <- NA
  }
  
  #------------------------------------------------
  # truncate clearance times so they cannot go past the next infection time or
  # the end of the study period
  
  # initial clearance times
  for (i in seq_len(n_inf)) {
    w <- which(w_inf[i,] == 1)
    clear_init[w] <- pmin(clear_init[w], t_inf[i])
  }
  clear_init <- pmin(clear_init, end_time)
  
  # new infection clearance times
  t_running <- rep(end_time, n_haplo)
  for (i in n_inf:1) {
    clear_inf[i,] <- pmin(clear_inf[i,], t_running)
    t_running[w_inf[i,] == 1] <- t_inf[i]
  }
  
  #------------------------------------------------
  # work out the true and observed haplo state at all observation times
  
  # define a matrix for holding the true infection state of every haplotype
  # (columns) at every timepoint (rows)
  state_true <- matrix(0, n_samp, n_haplo)
  
  # apply initial conditions
  for (j in which(w_init == 1)) {
    state_true[,j] <- (samp_time < clear_init[j])
  }
  
  # apply subsequent infections
  for (i in seq_along(t_inf)) {
    for (j in which(w_inf[i,] == 1)) {
      state_true[,j] <- state_true[,j] + (samp_time > t_inf[i]) & (samp_time < clear_inf[i,j])
    }
  }
  
  # draw the observed state from the true state by taking into account sensitivity
  state_obs <- state_true * matrix(rbinom(n_haplo * n_samp, 1, sens), nrow = n_samp)
  
  #------------------------------------------------
  # reformat state_obs into data.frame
  
  df_data <- data.frame(ind = ind_name,
                        haplo = rep(1:n_haplo, each = n_samp),
                        time = samp_time,
                        positive = as.vector(state_obs))
  
  #------------------------------------------------
  # return as list
  
  ret <- list(df_data = df_data)
  if (return_full) {
    ret <- c(ret, list(t_inf = t_inf,
                       w_init = w_init,
                       clear_init = clear_init,
                       w_inf = w_inf,
                       clear_inf = clear_inf,
                       state_true = state_true,
                       state_obs = state_obs,
                       samp_time = samp_time))
  }
  return(ret)
}

#------------------------------------------------
#' @title Simulate a complete cohort of individuals
#'
#' @description Simulate a cohort of individuals from the generative model. Runs
#'   \code{sim_ind()} multiple times and reformats observed output into a
#'   data.frame. Also keeps hold of the true states of all individuals for
#'   reference.
#'
#' @inheritParams sim_ind
#' @param n number of individuals in cohort.
#'
#' @importFrom dplyr bind_rows
#' @export

sim_cohort <- function(n, samp_time, haplo_freqs, lambda, decay_rate, sens) {
  
  # simulate each individual
  raw_list <- list()
  for (i in 1:n) {
    raw_list[[i]] <- sim_ind(samp_time, haplo_freqs, lambda[i], decay_rate, sens)
  }
  
  # get observed data into data.frame
  df_data <- mapply(function(i) {
    data.frame(ind = i,
               haplo = seq_along(haplo_freqs),
               time = rep(samp_time, each = length(haplo_freqs)),
               positive = as.vector(raw_list[[i]]$state_obs))
  }, seq_along(raw_list), SIMPLIFY = FALSE) |>
    bind_rows()
  
  # return as list
  ret <- list(df_data = df_data,
              raw_list = raw_list)
  return(ret)
}
