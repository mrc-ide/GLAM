#' @title R6 class representing a GLAM MCMC object
#'
#' @description 
#' Create a new glam_mcmc object to run an MCMC.
#' 
#' @import R6
#' @import dust
#' @export

glam_mcmc <- R6::R6Class(
  classname = "glam_mcmc",
  
  private = list(
    
    ### Private variables ###
    
    # data
    df_data = NULL,
    n_samp = NULL,
    n_haplos = NULL,
    haplo_freqs = NULL,
    obs_time_list = NULL,
    
    # program flow
    sample_called = FALSE,
    
    # MCMC parameters
    chains = NULL,
    rungs = NULL,
    max_infections = NULL,
    param_list = NULL,
    proposal_sd = NULL,
    n_proposal_sd = NULL,
    rng_list = NULL,
    iteration_counter = NULL,
    acceptance_counter = NULL,
    swap_acceptance_counter = NULL,
    duration = NULL,
    
    ### Private functions ###
    private_func = function(x) {
      print(x)
    }
  ),
  
  public = list(
    
    ### Public variables ###
    
    #' @field example_public_var example.
    example_public_var = NULL,
    
    ### Public functions ###
    
    #--------------------
    ## Load data
     
    #' @description
    #' Load data into the glam_mcmc object.
    #' 
    #' @param df_data a data.frame of data.
    initialize = function(df_data) {
      
      # input checks
      assert_dataframe(df_data)
      
      message("Loading data")
      
      # get sample info
      private$df_data <- df_data
      private$n_samp <- length(unique(df_data$ind))
      
      # get haplotype info
      private$n_haplos <- length(unique(df_data$haplo))
      private$haplo_freqs <- df_data |>
        group_by(haplo) |>
        summarise(haplo_freq = mean(positive)) |>
        arrange(haplo) |>
        pull(haplo_freq)
      
      # get observation time info
      df_data_time <- df_data |>
        group_by(ind) |>
        reframe(time = unique(time)) |>
        arrange(ind, time)
      
      obs_time_list <- split(df_data_time$time, f = df_data_time$ind)
      private$obs_time_list <- obs_time_list
      
      # print summary of data import
      message(sprintf("%s samples", private$n_samp))
      message(sprintf("%s haplotypes", private$n_haplos))
      if (length(unique(obs_time_list)) == 1) {
        message(sprintf("%s observation times, identical for all samples", length(obs_time_list[[1]])))
      } else {
        message("variable observation times between samples")
      }
      
    },
    
    #--------------------
    ## Define parameters
     
    #' @description
    #' Define model parameters.
    #' 
    #' @param chains the number of independent chains.
    #' @param rungs TODO
    #' @param lambda TODO
    #' @param theta TODO
    #' @param decay_rate TODO
    #' @param sens TODO
    #' @param n_infections TODO
    #' @param haplo_freqs TODO
    #' @param max_infections TODO
    params = function(chains = 5,
                      rungs = 1,
                      lambda = NULL,
                      theta = NULL,
                      decay_rate = NULL,
                      sens = NULL,
                      n_infections = NULL,
                      haplo_freqs = NULL,
                      max_infections = 5) {
      
      # input checks
      assert_single_pos_int(chains, zero_allowed = FALSE)
      
      # which parameters need updating
      lambda_fixed <- !is.null(lambda)
      theta_fixed <- !is.null(theta)
      decay_rate_fixed <- !is.null(decay_rate)
      sens_fixed <- !is.null(sens)
      n_infections_fixed <- !is.null(n_infections)
      
      # initialise parameters in nested list over chains and then rungs
      param_list <- list()
      for (i in 1:chains) {
        param_list[[i]] <- list()
        for (j in 1:rungs) {
          param_list[[i]][[j]] <- list(
            lambda = ifelse(lambda_fixed, lambda, runif(1)),
            theta = ifelse(theta_fixed, theta, runif(1)),
            decay_rate = ifelse(decay_rate_fixed, decay_rate, runif(1)),
            sens = ifelse(sens_fixed, sens, runif(1))
            )
          if (!n_infections_fixed) {
            param_list[[i]][[j]]$n_infections <- sample(max_infections, private$n_samp, replace = TRUE)
          }
        }
      }
      
      # initialise proposal standard deviations for all parameters
      proposal_sd <- list()
      n_proposal_sd <- 5
      for (i in 1:chains) {
        proposal_sd[[i]] <- list()
        for (j in 1:rungs) {
          proposal_sd[[i]][[j]] <- list(lambda = 1,
                                        theta = 1,
                                        decay_rate = 1,
                                        sens = 1,
                                        infection_time = 1)
        }
      }
      
      # initialise random number pointer for each chain
      rng_list <- dust::dust_rng_distributed_pointer(n_nodes = chains)
      
      # initialise counters
      iteration_counter <- create_chain_phase_list(chains = chains, base = 0)
      duration <- create_chain_phase_list(chains = chains, base = 0)
      acceptance_counter <- create_chain_phase_list(chains = chains,
                                                    base = matrix(
                                                      data = 0,
                                                      nrow = rungs, 
                                                      ncol = n_proposal_sd
                                                    ))
      swap_acceptance_counter <- create_chain_phase_list(chains = chains,
                                                    base = rep(0, rungs - 1))
      
      # overwrite haplo_freqs if defined manually
      if (!is.null(haplo_freqs)) {
        private$haplo_freqs <- haplo_freqs
      }
      
      # store values
      private$chains <- chains
      private$rungs <- rungs
      private$max_infections <- max_infections
      private$param_list <- param_list
      private$proposal_sd <- proposal_sd
      private$n_proposal_sd <- n_proposal_sd
      private$rng_list <- rng_list
      private$iteration_counter <- iteration_counter
      private$acceptance_counter <- acceptance_counter
      private$swap_acceptance_counter <- swap_acceptance_counter
      private$duration <- duration
    },
    
    #--------------------
    ## Run burn-in MCMC
     
    #' @description
    #' Run burn-in MCMC
    #' 
    #' @param iterations number of burn-in iterations
    #' @param target_acceptance all Metropolis-Hastings proposals will be
    #'   adaptively tuned to aim for this acceptance rate.
    #' @param silent if \code{TRUE} then suppress all console output. Default =
    #'   \code{FALSE}.
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE) {
      
      # check inputs
      assert_single_pos_int(iterations)
      assert_greq(iterations, 10)
      assert_single_bounded(target_acceptance, inclusive_left = FALSE, inclusive_right = FALSE)
      assert_single_logical(silent)
      
      # cannot run burnin after sampling
      if (private$sample_called) {
        stop("Cannot run burnin after sampling called")
      }
      
      # loop over chains
      chains <- private$chains
      for (chain in 1:chains) {
        
        message(sprintf("Running chain %s", chain))
        
        # run this chain
        output_raw <- mcmc_cpp(iterations,                                        # iterations
                               TRUE,                                              # burnin
                               private$iteration_counter[[chain]][["burn"]],      # iteration_counter_init
                               private$proposal_sd[[chain]],                      # proposal_sd
                               rep(1, private$rungs),                             # beta
                               private$rng_list[[chain]]                          # rng_ptr
        )
        
        # sync RNG
        private$rng_list[[chain]]$sync()
        
        # update counters
        #private$iteration_counter[[chain]][["burn"]] <- private$iteration_counter[[chain]][["burn"]] + iterations
        #private$acceptance_counter[[chain]][["burn"]] <- private$acceptance_counter[[chain]][["burn"]] + output_raw$acceptance_out
        #private$swap_acceptance_counter[[chain]][["burn"]] <- private$swap_acceptance_counter[[chain]][["burn"]] + output_raw$swap_acceptance_out
        
      } # end loop over chains
      
      #print(private$iteration_counter)
      #print(private$acceptance_counter)
      #print(private$swap_acceptance_counter)
    },
    
    #--------------------
    ## Print
     
    #' @description
    #' Print MCMC object summary
    print = function() {
      message("GLAM MCMC object")
      
      # print summary of data
      message(sprintf("%s samples", private$n_samp))
      message(sprintf("%s haplotypes", private$n_haplos))
      if (length(unique(private$obs_time_list)) == 1) {
        message(sprintf("%s observation times, identical for all samples", length(private$obs_time_list[[1]])))
      } else {
        message("variable observation times between samples")
      }
      
      # return invisibly
      invisible(self)
    }
  )
)
