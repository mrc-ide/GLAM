
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

