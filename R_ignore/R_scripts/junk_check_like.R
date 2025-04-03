
# junk script - calculates likelihood over a continuous range of parameter values.

# get probability of W matrix. Note, this follows the generative process from a
# ZTBS, but could be made much simpler by using the final probability density,
# which has a simple form
get_wmat_prob <- function(w_mat, theta) {
  if (is.null(w_mat)) {
    return(0)
  }
  q <- 1.0 - exp(-theta*haplo_freqs)
  a <- rep(NA, n_haplos)
  for (i in seq_along(a)) {
    a[i] <- prod(1 - q[i:n_haplos])
  }
  ret <- 0
  for (inf_index in 1:nrow(w_mat)) {
    no_pos <- 1
    for (i in 1:n_haplos) {
      if (w_mat[inf_index,i] == TRUE) {
        ret <- ret + log(q[i] / (1 - no_pos*a[i]))
        no_pos <- 0
      } else {
        ret <- ret + log((1 - q[i] - no_pos*a[i]) / (1 - no_pos*a[i]))
      }
    }
  }
  return(ret)
}

# get data as list
g <- glam_mcmc$new(df_data = sim2$df_data)
data_list <- g$.__enclos_env__$private$list_data
w <- 1:24
data_list <- data_list[w]
infection_times <- mapply(function(x) x$t_inf, sim2$raw_list, SIMPLIFY = FALSE)[w]
w_list <- mapply(function(x) x$w_inf, sim2$raw_list, SIMPLIFY = FALSE)[w]
n_samp <- length(data_list)
n_haplos <- length(haplo_freqs)



get_ll <- function(theta) {
  ret_mat <- matrix(NA, n_samp, n_haplos)
  for (indiv_i in 1:n_samp) {
    w_mat <- w_list[[indiv_i]]
    
    for (haplo_i in 1:n_haplos) {
      
      dat_bool <- data_list[[indiv_i]][,haplo_i][[1]]
      df_algo1 <- algo1(haplo_i = haplo_i,
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
      
      ret_mat[indiv_i, haplo_i] <- last(df_algo1$A) + last(df_algo1$B)
    }
  }
  
  ret_vec <- rep(NA, n_samp)
  for (indiv_i in 1:n_samp) {
    w_mat <- w_list[[indiv_i]]
    ret_vec[indiv_i] <- get_wmat_prob(w_mat, theta = theta)
  }
  
  ret <- 0
  ret <- ret + sum(ret_vec)
  ret <- ret + sum(log(ret_mat))
  
  return(ret)
}

theta_vec <- seq(0.02, 0.4, 0.02)
y_vec <- rep(NA, length(theta_vec))
for (i in seq_along(theta_vec)) {
  message(sprintf("%s of %s", i, length(theta_vec)))
  y_vec[i] <- get_ll(theta_vec[i])
}
plot(theta_vec, exp(y_vec - max(y_vec)))

