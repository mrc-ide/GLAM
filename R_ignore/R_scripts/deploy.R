# deploy.R
#
# Author: Bob Verity
# Date: 2025-04-03
#
# Purpose:
# Sandbox for GLAM package.
#
# ------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(dust)
set.seed(2)

# -----------------------
# SIMULATE DATA

# set simulation parameters
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 2
decay_rate <- 0.1
sens <- 0.9
max_infections <- 20

# simulate a single sample and plot
sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim1)

# simulate a cohort of multiple samples
sim2 <- sim_cohort(n = 10,
                   samp_time = samp_time,
                   haplo_freqs = haplo_freqs,
                   lambda = lambda,
                   theta = theta,
                   decay_rate = decay_rate,
                   sens = sens,
                   n_inf = n_infections_true,
                   return_full = TRUE)


# -----------------------
# RUN MCMC

# create a new MCMC object and load data
g <- glam_mcmc$new(df_data = sim2$df_data)

# initialize. If parameters are set to NULL, they are estimated. If they are
# given a value, they take this fixed value. Useful for giving the inference
# method the correct true value of some parameters, to diagnose how well it
# estimates others.
g$init(start_time = 0,
       end_time = 10,
       haplo_freqs = haplo_freqs,
       lambda = lambda,
       theta = theta,
       decay_rate = decay_rate,
       sens = sens,
       n_infections = mapply(function(x) length(x$t_inf), sim2$raw_list),
       infection_times = NULL,
       max_infections = 20,
       chains = 1, rungs = 1)

# run burn-in and sampling
g$burn(iterations = 1e2) |>
  system.time()

g$sample(iterations = 1e3) |>
  system.time()

# -----------------------
# EXTRACT OUTPUTS AND PLOT RESULTS

# extract outputs into data.frames
df_global <- g$get_output_global()
df_n_infections <- g$get_output_n_infections()
df_infection_times <- g$get_output_infection_times()

# bivariate plots of scalar parameters
plot(df_global$decay_rate, df_global$sens, pch = 20, col = "#00000010", xlim = c(0, 0.5), ylim = c(0, 1))
points(decay_rate, sens, pch = 4, cex = 2, col = 2, lwd = 3)

plot(df_global$lambda, df_global$theta, pch = 20, col = "#00000010", xlim = c(0, 0.5), ylim = c(0, 2))
points(lambda, theta, pch = 4, cex = 2, col = 2, lwd = 3)


# trace plot of a single infection time
df_infection_times |>
  filter(individual == 1) |>
  filter(infection == 1) |>
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, col = chain)) +
  ylim(c(0, max_infections))


# density plot of all infection times
df_infection_times |>
  filter(phase == "sampling") |>
  group_by(individual, infection) |>
  do({
    d <- density(.$value)
    data.frame(x = d$x, y = d$y)
  }) |>
  ungroup() |>
  ggplot(aes(x = x, y = y, color = infection)) +
  geom_vline(aes(xintercept = time), col = grey(0.5), data = df_inf_time_true) +
  geom_line() +
  facet_wrap(~ individual) +
  labs(x = "Value", y = "Density", color = "Infection") +
  theme_bw() +
  theme(panel.grid.major = element_blank())


