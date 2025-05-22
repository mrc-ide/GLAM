# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("stats", "end", "start", "time", "haplo", "positive", "t_inf",
                           "t_start", "t_clear", "state_true", "state_obs", "state_combined"))
}

#------------------------------------------------
#' @title Compute time-varying protection from repeated treatments
#'
#' @description
#' This function calculates protective efficacy at a series of observation times, based on a set of treatment times. Protection is modelled as a sigmoidal (logistic) waning curve that resets to `max_protection` at each treatment time, and declines from that point forward until the next treatment.
#'
#' @param observation_times A numeric vector of observation times (e.g., in days) at which protection is evaluated.
#' @param treatment_times A numeric vector of treatment times (must be on the same time scale as `observation_times`).
#' @param treatment_duration_90 Time at which protection is still at 90% of `max_protection` after a treatment. Must be less than `treatment_duration_50`.
#' @param treatment_duration_50 Time (in the same units as `observation_times`) at which protection wanes to 50% of `max_protection` after a treatment. Must be greater than `treatment_duration_90`.
#' @param max_protection The maximum protection conferred immediately after treatment (default is 1.0).
#'
#' @return A numeric vector of the same length as `observation_times`, containing the estimated protection at each time point.
#'
#' @details
#' The protection at any observation time is computed using a logistic function:
#' \deqn{Protection(t) = \frac{P_{\text{max}}}{1 + \exp(k \cdot (t - t_{50}))}}
#' where \eqn{k = \log(9)/(t_{50} - t_{90})} ensures that protection equals 90% at `treatment_duration_90` and 50% at `treatment_duration_50`
#' after each treatment. If no treatment has occurred by an observation time, protection is zero.
#'
#' @examples
#' obs_times <- seq(0, 100, by = 1)
#' tx_times <- c(10, 40, 70)
#' protection <- get_treatment_protection(
#'   observation_times = obs_times,
#'   treatment_times = tx_times,
#'   treatment_duration_90 = 10,
#'   treatment_duration_50 = 20
#' )
#' plot(obs_times, protection, type = "l", main = "Waning Protection", ylab = "Protection")
#' abline(v = tx_times, col = "blue", lty = 2)
#'
#' @export
get_treatment_protection <- function(observation_times,
                                     treatment_times,
                                     treatment_duration_90,
                                     treatment_duration_50,
                                     max_protection = 1.0) {
  
  # checks
  if (treatment_duration_90 >= treatment_duration_50) {
    stop("treatment_duration_90 must be less than treatment_duration_50")
  }
  
  # steepness of waning
  k <- log(9) / (treatment_duration_50 - treatment_duration_90)
  
  # sort treatment times just in case
  treatment_times <- sort(treatment_times)
  
  # initialize output vector
  protection <- numeric(length(observation_times))
  
  for (i in seq_along(observation_times)) {
    t_obs <- observation_times[i]
    
    # find most recent treatment before or at observation time
    past_treatments <- treatment_times[treatment_times <= t_obs]
    
    if (length(past_treatments) == 0) {
      protection[i] <- 0  # no treatment yet
    } else {
      t_last <- max(past_treatments)
      time_since_treatment <- t_obs - t_last
      protection[i] <- max_protection / (1 + exp(k * (time_since_treatment - treatment_duration_50)))
    }
  }
  
  return(protection)
}
