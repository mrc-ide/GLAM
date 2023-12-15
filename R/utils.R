
#------------------------------------------------
# draw from zero-truncated Poisson with rate parameter lambda (expectation
# lambda/(1-exp(-lambda)))
#' @importFrom stats qpois rmultinom
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
