#------------------------------------------------
#' @title Extract the MCMC draws from a glam_mcmc object
#'
#' @description Raw MCMC draws are stored within an object of class
#'   \code{glam_mcmc} in a series of lists. This makes it easier to extend the
#'   MCMC as needed, but is not the most convenient format for exploring or
#'   plotting. This function extracts these draws into a convenient long-form
#'   data.frame for use in downstream functions.
#'
#' @param x an object of class \code{glam_mcmc}, as produced by
#'   \code{glam_mcmc$new()}.
#'
#' @export

get_mcmc_draws <- function(x) {
  
  # check inputs
  assert_class(x, "glam_mcmc")
  
  
  
  # return
  return(NULL)
}
