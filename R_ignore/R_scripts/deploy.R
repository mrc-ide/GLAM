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

#set.seed(2)

# ------------------------------------------------------------------

# define fixed parameters
n_cohort <- 10
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.4
theta <- 1
decay_rate <- 0.2
sens <- 0.99
max_infections <- 20

# set some simulation parameters
n_infections <- rep(4, n_cohort)
#n_infections <- NULL

sim1 <- sim_cohort(n = n_cohort,
                   samp_time = samp_time,
                   haplo_freqs = haplo_freqs,
                   lambda = lambda,
                   theta = theta,
                   decay_rate = decay_rate,
                   sens = sens,
                   n_inf = n_infections,
                   return_full = TRUE)

#plot_sim(sim1)

#plot_data(sim1$df_data)

inf_times <- mapply(function(x) x$t_inf, sim1$raw_list, SIMPLIFY = FALSE)
inf_times <- NULL

# -----------------------

g <- glam_mcmc$new(df_data = sim1$df_data)

g$init(start_time = 0, 
       end_time = 10, 
       chains = 1, 
       rungs = 1, 
       haplo_freqs = haplo_freqs,
       lambda = lambda, 
       theta = NULL, 
       decay_rate = NULL, 
       sens = NULL,
       n_infections = n_infections, 
       infection_times = inf_times,
       max_infections = max_infections, 
       w_list = NULL)

g$burn(iterations = 1e3)
t0 <- Sys.time()
g$sample(iterations = 1e4)
Sys.time() - t0

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
par(mfrow = c(1,2))
plot(df_global$decay_rate, df_global$sens, pch = 20, col = "#00000010", xlim = c(0, 1), ylim = c(0, 1))
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
  t <- sim1$raw_list[[i]]$t_inf
  if (length(t) == 0) {
    return(NULL)
  }
  data.frame(individual = i,
             infection = seq_along(t),
             time = t)
}, seq_along(sim1$raw_list), SIMPLIFY = FALSE) |>
  bind_rows()

d |>
  ggplot() + theme_bw() +
  geom_line(aes(x = x, y = y, col = infection)) +
  geom_vline(aes(xintercept = time), linetype = "dashed", data = df_inf_time_true) +
  facet_wrap(~individual)



