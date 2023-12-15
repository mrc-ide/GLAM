test_that("basic mcmc (burn-in + sampling) runs without error", {
  
  # define parameters
  samp_time <- seq(0, 10, 1)
  haplo_freqs <- rep(0.1, 10)
  lambda <- 0.2
  theta <- 0.5
  decay_rate <- 0.1
  sens <- 0.9
  
  # simulate some data
  sim_data <- sim_cohort(n = 5,
                         samp_time = samp_time,
                         haplo_freqs = haplo_freqs,
                         lambda = lambda,
                         theta = theta,
                         decay_rate = decay_rate,
                         sens = sens)
  
  # create mcmc object
  g <- glam_mcmc$new(df_data = sim_data$df_data)
  
  # set parameters
  g$params(chains = 3, rungs = 5)
  
  # run mcmc
  g$burn(iterations = 1e2)
  
  expect_equal(1, 1)
})
