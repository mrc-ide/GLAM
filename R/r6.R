#' @title R6 class representing a GLAM MCMC object
#'
#' @description 
#' Create a new glam_mcmc object to run an MCMC.
#' 
#' @export

glam_mcmc <- R6::R6Class(
  classname = "glam_mcmc",
  
  private = list(
    
    ### Private variables ###
    
    # data
    df_data = NULL,
    n_samp = NULL,
    
    # program flow
    sampling_called = FALSE,
    
    ### Private functions ###
    private_func = function(x) {
      print(x)
    }
  ),
  
  public = list(
    
    ### Public variables ###
    example_public_var = NULL,
    
    ### Public functions ###
    
    #--------------------
    #' @title Load data
    #' 
    #' @description
    #' Load data into the glam_mcmc object.
    #' 
    #' @param df_data a data.frame of data.
    initialize = function(df_data) {
      
      # input checks
      assert_dataframe(df_data)
      
      private$df_data <- df_data
      private$n_samp <- nrow(df_data)
    },
    
    #--------------------
    #' @title Run burn-in MCMC
    #' 
    #' @description
    #' Run burn-in MCMC
    #' 
    #' @param iterations number of burn-in iterations
    #' @param target_acceptance all Metropolis-Hastings proposals will be
    #'   adaptively tuned to aim for this acceptance rate.
    #' @param silent if \code{TRUE} then suppress all console output. Default =
    #'   \code{FALSE}.
    burnin = function(iterations, target_acceptance = 0.44, silent = FALSE) {
      
      # check inputs
      assert_single_pos_int(iterations)
      assert_greq(iterations, 10)
      assert_single_bounded(target_acceptance, inclusive_left = FALSE, inclusive_right = FALSE)
      assert_single_logical(silent)
      
      # cannot run burnin after sampling
      if (private$sampling_called) {
        stop("Cannot run burnin after sampling called")
      }
      
      mcmc_cpp()
      
    },
    
    #--------------------
    #' @title Print
    #' 
    #' @description
    #' Print MCMC object summary
    print = function() {
      message("GLAM MCMC object")
      
      # return invisibly
      invisible(self)
    }
  )
)
