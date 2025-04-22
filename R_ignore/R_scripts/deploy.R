# deploy.R
#
# Author: Bob Verity
# Date: 2023-03-02
#
# Purpose:
# Sandbox for GLAM package.
#
# ------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(dust)

samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 2
decay_rate <- 0.1
sens <- 0.9
max_infections <- 20

set.seed(2)
sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim1)

sim2 <- sim_cohort(n = 10,
                   samp_time = samp_time,
                   haplo_freqs = haplo_freqs,
                   lambda = lambda,
                   theta = theta,
                   decay_rate = decay_rate,
                   sens = sens,
                   n_inf = n_infections_true,
                   return_full = TRUE)

lambda_true <- NULL
theta_true <- NULL
decay_rate_true <- NULL
sens_true <- NULL
n_infections_true <- NULL
infection_times_true <- NULL
w_list_true <- NULL

#sim2$df_data |>
#  filter(ind == "ind2") |>
#  pivot_wider(names_from = time, values_from = positive)

#lambda_true <- lambda
#theta_true <- theta
#decay_rate_true <- decay_rate
#sens_true <- sens
n_infections_true <- mapply(function(x) length(x$t_inf), sim2$raw_list)
#infection_times_true <- mapply(function(x) x$t_inf, sim2$raw_list, SIMPLIFY = FALSE)
#w_list_true <- mapply(function(x) x$w_inf, sim2$raw_list, SIMPLIFY = FALSE)

# -----------------------

g <- glam_mcmc$new(df_data = sim2$df_data)

g$init(start_time = 0, 
       end_time = 10, 
       chains = 1, 
       rungs = 1, 
       haplo_freqs = haplo_freqs,
       lambda = lambda_true, 
       theta = theta_true, 
       decay_rate = decay_rate_true, 
       sens = sens_true,
       n_infections = n_infections_true, 
       infection_times = infection_times_true,
       max_infections = max_infections, 
       w_list = w_list_true)

g$burn(iterations = 1e2)
system.time(
  g$sample(iterations = 1e3)
  )

# -----------------------

df_global <- g$get_output_global()
df_n_infections <- g$get_output_n_infections()
df_infection_times <- g$get_output_infection_times()


# plot scalar parameters
plot(df_global$decay_rate, ylim = c(0, 2))
abline(h = decay_rate, col = 2, lwd = 3)

plot(df_global$sens, ylim = c(0, 1))
abline(h = sens, col = 2, lwd = 3)

# bivariate plots
plot(df_global$decay_rate, df_global$sens, pch = 20, col = "#00000010", xlim = c(0, 0.5), ylim = c(0, 1))
points(decay_rate, sens, pch = 4, cex = 2, col = 2, lwd = 3)

plot(df_global$lambda, df_global$theta, pch = 20, col = "#00000010", xlim = c(0, 0.5), ylim = c(0, 2))
points(lambda, theta, pch = 4, cex = 2, col = 2, lwd = 3)


# plot n_infections
df_n_infections |>
  filter(ind == 1) |>
  ggplot() + theme_bw() +
  geom_line(aes(x = iteration, y = value, col = chain)) +
  ylim(c(0, max_infections))

# plot single infection time
df_infection_times |>
  filter(individual == 1) |>
  filter(infection == 1) |>
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, col = chain)) +
  ylim(c(0, max_infections))


# density plot infection times
d <- df_infection_times |>
  group_by(individual, infection) |>
  reframe(y = density(value, from = 0, to = 10, n = 101, bw = 0.05)$y) |>
  as.data.frame()
d$x <- seq(0, 10, l = 101)

df_inf_time_true <- mapply(function(i) {
  t <- sim2$raw_list[[i]]$t_inf
  if (length(t) == 0) {
    return(NULL)
  }
  data.frame(individual = i,
             infection = seq_along(t),
             time = t)
}, seq_along(sim2$raw_list), SIMPLIFY = FALSE) |>
  bind_rows()

d |>
  ggplot() + theme_bw() +
  geom_line(aes(x = x, y = y, col = infection)) +
  geom_vline(aes(xintercept = time), linetype = "dashed", data = df_inf_time_true) +
  facet_wrap(~individual)



