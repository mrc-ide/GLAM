
#------------------------------------------------
# log-normal distribution parametrized in terms of mean and SD
#' @importFrom stats dlnorm
#' @noRd
dlnorm_reparam <- function(x, mean, sd, return_log = FALSE) {
  sigma_sq <- log(sd^2 / mean^2 + 1)
  dlnorm(x, meanlog = log(mean) - sigma_sq / 2, sdlog = sqrt(sigma_sq), log = return_log)
}

#------------------------------------------------
# gamma distribution parametrized in terms of mean and SD
#' @importFrom stats dgamma
#' @noRd
dgamma_reparam <- function(x, mean, sd, return_log = FALSE) {
  dgamma(x, shape = mean^2 / sd^2, rate = mean / sd^2, log = return_log)
}

#------------------------------------------------
# zero-truncated Poisson with rate parameter lambda (expectation
# lambda/(1-exp(-lambda)))
#' @importFrom stats qpois
#' @noRd
rztpois <- function(n, lambda) {
  if (lambda <= 0) {
    stop("lambda must be positive")
  }
  qpois(runif(n, exp(-lambda), 1), lambda)
}

#------------------------------------------------
# create nested list of phases in chains
#' @noRd
create_chain_phase_list <- function(chains, base) {
  lapply(1:chains, function(x) {
    list(
      tune = base,
      burn = base,
      sample = base
    )
  })
}
