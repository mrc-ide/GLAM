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

set.seed(2)

# ------------------------------------------------------------------

# define fixed parameters
n_cohort <- 12
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 4
decay_rate <- 0.5
sens <- 0.95
max_infections <- 10

# simulate cohort
sim1 <- sim_cohort(n = n_cohort,
                   samp_time = samp_time,
                   haplo_freqs = haplo_freqs,
                   lambda = lambda,
                   theta = theta,
                   decay_rate = decay_rate,
                   sens = sens,
                   n_inf = NULL,
                   return_full = TRUE)

# plot simulated data
plot_sim(sim1)

#plot_data(sim1$df_data)

# -----------------------

# extract true simulated number of infections
df_n_inf_true <- mapply(function(x, i) {
  data.frame(ind = i,
             n_inf = length(x$t_inf))
}, sim1$raw_list, seq_along(sim1$raw_list), SIMPLIFY = FALSE) |>
  bind_rows()

# extract true simulated infection times
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

# -----------------------

g <- glam_mcmc$new(df_data = sim1$df_data)

g$init(start_time = 0, 
       end_time = 10, 
       chains = 1, 
       rungs = 1, 
       haplo_freqs = haplo_freqs,
       lambda = NULL, 
       theta = NULL, 
       decay_rate = NULL, 
       sens = NULL,
       n_infections = NULL, 
       infection_times = NULL,
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
df_global |>
  filter(phase == "sampling") |>
  ggplot() + theme_bw() +
  geom_point(aes(x = lambda, y = theta), size = 0.2) +
  xlim(c(0, 1)) + ylim(c(0, 10)) +
  geom_point(x = lambda, y = theta, col = "red", pch = 4, size = 3, stroke = 2)

df_global |>
  filter(phase == "sampling") |>
  ggplot() + theme_bw() +
  geom_point(aes(x = decay_rate, y = sens), size = 0.2) +
  xlim(c(0, 1)) + ylim(c(0, 1)) +
  geom_point(x = decay_rate, y = sens, col = "red", pch = 4, size = 3, stroke = 2)


# plot n_infections
df_n_infections |>
  filter(phase == "sampling") |>
  ggplot() + theme_bw() +
  geom_histogram(aes(x = value), binwidth = 1, boundary = -0.5) +
  geom_vline(aes(xintercept = n_inf), linetype = "dashed", data = df_n_inf_true) +
  facet_wrap(~ind) +
  scale_x_continuous(breaks = 0:max_infections, limits = c(0, max_infections)) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Number of infections")


# density plot infection times
df_density <- df_infection_times |>
  filter(phase == "sampling") |>
  group_by(individual, infection) |>
  reframe(x = seq(0, 10, l = 101),
          y = density(value, from = 0, to = 10, n = 101, bw = 0.05)$y) |>
  as.data.frame() |>
  ungroup() |>
  left_join(df_infection_times |>
              group_by(individual, infection) |>
              summarise(w = n(),
                        .groups = "drop"),
            by = join_by(individual, infection))

df_density |>
  ggplot() + theme_bw() +
  geom_area(aes(x = x, y = w*y, fill = infection), colour = "black", size = 0.1, position = "stack") +
  geom_vline(aes(xintercept = time), linetype = "dashed", data = df_inf_time_true) +
  facet_wrap(~individual) +
  xlab("Time") + ylab("Stacked probability")



