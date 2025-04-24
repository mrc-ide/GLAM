# debug_Jess_data.R
#
# Author: Bob Verity
# Date: 2025-04-22
#
# Inputs: R_ignore/data/PBC_cid_526_for_GLAM.csv
#
# Outputs: (none)
#
# Purpose:
# I had some trouble initially running a small dataset provided by Jess Briggs. There was some issue with how the input data was being interpreted (i.e. the Type) that was causing a crash that was very difficult to debug. This bug now seems to have inexplicably gone away (!?), but keeping this deploy script in case it comes back.
#
# ------------------------------------------------------------------

library(tidyverse)

max_t <- 100

# read in and reformat Jess data
dat <- read.csv("R_ignore/data/PBC_cid_526_for_GLAM.csv")
names(dat) <- c("ind", "haplo", "time", "positive")
dat <- dat |>
  filter(time <= max_t) |>
  mutate(ind = sprintf("ind_%s", ind),
         haplo = match(haplo, unique(haplo)),
         time = as.numeric(time),
         positive = as.numeric(positive))

# simple plot of observed haplotypes
plot_data(dat)

# -----------------------

# run model with contrived parameters. Only infection times being estimated
n_haplos <- length(unique(dat$haplo))
haplo_freqs <- rep(1 / n_haplos, n_haplos)

g <- glam_mcmc$new(df_data = dat)

g$init(start_time = 0, 
       end_time = max_t, 
       chains = 1, 
       rungs = 1, 
       haplo_freqs = haplo_freqs,
       lambda = 1, 
       theta = NULL, 
       decay_rate = NULL, 
       sens = NULL,
       n_infections = 3, 
       infection_times = NULL,
       max_infections = 10,
       w_list = NULL)

g$burn(iterations = 1e3)
g$sample(1e4)

# -----------------------

# plot results
df_infection_times <- g$get_output_infection_times()

# density plot infection times
d <- df_infection_times |>
  group_by(individual, infection) |>
  reframe(y = density(value, from = 0, to = max_t, n = 101, bw = 0.05)$y) |>
  as.data.frame()
d$x <- seq(0, 10, l = 101)

d |>
  ggplot() + theme_bw() +
  geom_line(aes(x = x, y = y, col = infection)) +
  facet_wrap(~individual)
