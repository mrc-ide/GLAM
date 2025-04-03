
library(tidyverse)

dat <- read.csv("/Users/rverity/Desktop/PBC_cid_526_for_GLAM.csv")
names(dat) <- c("ind", "haplo", "time", "positive")
dat <- dat |>
  filter(time < 10) |>
  mutate(ind = sprintf("ind_%s", ind),
         haplo = match(haplo, unique(haplo)),
         time = as.numeric(time),
         positive = as.numeric(positive))

#dat <- sim2$df_data

n_haplos <- length(unique(dat$haplo))
haplo_freqs <- rep(1 / n_haplos, n_haplos)

g <- glam_mcmc$new(df_data = dat)

g$init(start_time = 0, end_time = 10, chains = 1, rungs = 1, haplo_freqs = haplo_freqs,
       lambda = 1, theta = 0.1, decay_rate = 1, sens = 0.9,
       n_infections = 2, infection_times = NULL,
       max_infections = 10, w_list = NULL)

#private2 <- g$.__enclos_env__$private
private <- g$.__enclos_env__$private
class(private$obs_time_list[[1]])

g$burn(iterations = 1e2)

g

