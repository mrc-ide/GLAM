
#------------------------------------------------
#' GLAM MCMC Class
#'
#' An R6 generator for running Markov Chain Monte Carlo (MCMC) in GLAM.
#'
#' @name glam_mcmc
#' @docType class
#' @format An \code{\link[R6]{R6Class}} generator object with methods:
#' @import R6
#' @import dust
#' @importFrom tidyr pivot_wider
#' @import progress
#' @export
glam_mcmc <- R6::R6Class(
  classname = "glam_mcmc",
  
  public = list(
    
    ### Public variables ###
    # (none)
    
    ### Public functions ###
    
    #--------------------
    #' @title Initialize the object with haplotype‐time series data
    #'
    #' @description Performs input validation on `df_data`.
    #'
    #' @param df_data A data.frame or tibble. Must contain the columns:
    #'   - `ind` (identifier for each sample/individual),
    #'   - `haplo` (haplotype label),
    #'   - `time` (observation time point),
    #'   - `positive` (binary indicator: 0 or 1).
    #' @param silent Logical scalar. If `FALSE` (the default), informational
    #'   messages about loading progress and data summary are printed to the
    #'   console. If `TRUE`, message printing is suppressed.
    #'
    #' @return Invisibly returns `NULL`. This is an R6 “initialize” method and
    #'   is used for its side effects on the object’s private fields.
    initialize = function(df_data, silent = FALSE) {
      
      # input checks
      assert_dataframe(df_data)
      assert_in(c("ind", "haplo", "time", "positive"), names(df_data),
                message = "input data must have columns: ind, haplo, time, positive")
      assert_in(df_data$positive, c(0,1), message = "`positive` column must contain 0 or 1")
      
      if (!silent) {
        message("Loading data")
      }
      
      # get sample info
      private$df_data <- df_data
      private$n_samp <- length(unique(df_data$ind))
      
      # convert to list (this is how we will pass to C++)
      private$list_data <- df_data |>
        dplyr::arrange(time) |>
        tidyr::pivot_wider(names_from = haplo, values_from = positive) |>
        dplyr::select(-time) |>
        dplyr::group_split(ind, .keep = FALSE)
      
      # get haplotype info
      private$n_haplos <- length(unique(df_data$haplo))
      private$haplo_freqs <- df_data |>
        dplyr::group_by(haplo) |>
        dplyr::summarise(haplo_freq = mean(positive)) |>
        dplyr::arrange(haplo) |>
        dplyr::pull(haplo_freq)
      
      # get observation time info
      df_data_time <- df_data |>
        dplyr::group_by(ind) |>
        dplyr::reframe(time = unique(time)) |>
        dplyr::arrange(ind, time)
      
      obs_time_list <- split(df_data_time$time, f = df_data_time$ind)
      private$obs_time_list <- obs_time_list
      
      # print summary of data import
      if (!silent) {
        message(sprintf("%s samples", private$n_samp))
        message(sprintf("%s haplotypes", private$n_haplos))
        if (length(unique(obs_time_list)) == 1) {
          message(sprintf("%s observation times, identical for all samples", length(obs_time_list[[1]])))
        } else {
          message("variable observation times between samples")
        }
      }
      
    },
    
    #--------------------
    #' @description
    #' Define model parameters.
    #' 
    #' @param start_time TODO
    #' @param end_time TODO
    #' @param chains the number of independent chains.
    #' @param rungs TODO
    #' @param max_infections TODO
    #' @param lambda TODO
    #' @param theta TODO
    #' @param decay_rate TODO
    #' @param sens TODO
    #' @param n_infections TODO
    #' @param infection_times TODO
    #' @param haplo_freqs TODO
    #' @param w_list TODO
    init = function(start_time,
                    end_time,
                    chains = 5,
                    rungs = 1,
                    max_infections = 5,
                    lambda = NULL,
                    theta = NULL,
                    decay_rate = NULL,
                    sens = NULL,
                    n_infections = NULL,
                    infection_times = NULL,
                    haplo_freqs = NULL,
                    w_list = NULL
                    ) {
      
      # check inputs
      assert_single_numeric(start_time)
      obs_time_range <- range(unlist(private$obs_time_list))
      assert_leq(start_time, obs_time_range[1], message = sprintf("start_time must be less than or equal to the first observation time (%s)", obs_time_range[1]))
      assert_single_numeric(end_time)
      assert_gr(end_time, start_time)
      assert_greq(end_time, obs_time_range[2], message = sprintf("end_time must be greater than or equal to the last observation time (%s)", obs_time_range[2]))
      assert_single_pos_int(chains, zero_allowed = FALSE)
      assert_single_pos_int(rungs, zero_allowed = FALSE)
      assert_single_pos_int(max_infections, zero_allowed = FALSE)
      if (!is.null(lambda)) {
        assert_single_pos(lambda, zero_allowed = TRUE)
      }
      if (!is.null(theta)) {
        assert_single_pos(theta, zero_allowed = FALSE)
      }
      if (!is.null(decay_rate)) {
        assert_single_pos(decay_rate, zero_allowed = TRUE)
      }
      if (!is.null(sens)) {
        assert_single_bounded(sens)
      }
      if (!is.null(n_infections)) {
        assert_vector_int(n_infections)
        assert_length(n_infections, private$n_samp)
        assert_bounded(n_infections, left = 0, right = max_infections)
      }
      if (!is.null(infection_times)) {
        assert_list(infection_times)
        assert_length(infection_times, private$n_samp, message = sprintf("infection_times must be of length %s to match the number of samples found in the data", private$n_samp))
        if (is.null(n_infections)) {
          message("Defining n_infections from infection_times")
          n_infections <- mapply(length, infection_times) # define n_infections from infection_times
        } else {
          assert_eq(mapply(length, infection_times), n_infections, message = "If both n_infections and infection_times are defined then the lengths of infection_times must match the values in n_infections")
        }
        mapply(function(x) {
          assert_vector_bounded(x, left = start_time, right = end_time, message = "Infection times must be within the window defined by start_time and end_time")
        }, infection_times)
      }
      if (!is.null(haplo_freqs)) {
        assert_vector_bounded(haplo_freqs)
        assert_length(haplo_freqs, private$n_haplos, message = sprintf("Must define %s haplotype frequencies to match the number found in the data", private$n_haplos))
        assert_eq(sum(haplo_freqs), 1)
      }
      if (!is.null(w_list)) {
        assert_list(w_list)
        assert_length(w_list, private$n_samp, message = sprintf("w_list must be of length %s to match the number of samples found in the data", private$n_samp))
        if (is.null(n_infections)) {
          message("Defining n_infections from w_list")
          n_infections <- mapply(nrow, w_list) # define n_infections from w_list
        } else {
          assert_eq(mapply(function(x) ifelse(is.null(x), 0, nrow(x)), w_list), n_infections, message = "If both n_infections and w_list are defined then the number of rows in w_list must match the values in n_infections")
        }
        mapply(function(x) {
          if (!is.null(x)) {
            assert_matrix(x)
            assert_logical(x)
          }
        }, w_list)
      }
      
      # reformat w_list to list of lists, not list of matrices
      if (!is.null(w_list)) {
        w_list <- mapply(function(x) {
          if (is.null(x)) {
            list()
          } else {
            apply(x, 1, function(y) y, simplify = FALSE)
          }
        }, w_list, SIMPLIFY = FALSE)
      }
      
      # which parameters need updating
      lambda_fixed <- !is.null(lambda)
      theta_fixed <- !is.null(theta)
      decay_rate_fixed <- !is.null(decay_rate)
      sens_fixed <- !is.null(sens)
      n_infections_fixed <- !is.null(n_infections)
      infection_times_fixed <- !is.null(infection_times)
      w_list_fixed <- !is.null(w_list)
      
      param_update_list <- list(lambda_fixed = lambda_fixed,
                                theta_fixed = theta_fixed,
                                decay_rate_fixed = decay_rate_fixed,
                                sens_fixed = sens_fixed,
                                n_infections_fixed = n_infections_fixed,
                                infection_times_fixed = infection_times_fixed,
                                w_list_fixed = w_list_fixed)
      
      # initialise parameters in nested list over chains and then rungs
      param_list <- list()
      for (i in 1:chains) {
        param_list[[i]] <- list()
        for (j in 1:rungs) {
          lambda <- ifelse(lambda_fixed, lambda, runif(1))
          theta <- ifelse(theta_fixed, theta, runif(1))
          decay_rate <- ifelse(decay_rate_fixed, decay_rate, runif(1))
          sens <- ifelse(sens_fixed, sens, runif(1))
          if (n_infections_fixed) {
            n_infections <- as.integer(n_infections)
          } else {
            n_infections <- sample(max_infections, private$n_samp, replace = TRUE)
          }
          if (!infection_times_fixed) {
            infection_times <- list()
            for (k in 1:private$n_samp) {
              infection_times[[k]] <- sort(runif(n_infections[k], min = start_time, max = end_time))
            }
          }
          if (!w_list_fixed) {
            w_list <- mapply(function(n) {
              replicate(n, sample(c(TRUE, FALSE), private$n_haplos, replace = TRUE), simplify = FALSE)
            }, n_infections, SIMPLIFY = FALSE)
          }
          
          param_list[[i]][[j]] <- list(
            lambda = lambda,
            theta = theta,
            decay_rate = decay_rate,
            sens = sens,
            n_infections = n_infections,
            infection_times = infection_times,
            w_list = w_list
            )
        }
      }
      
      # initialise proposal standard deviations for all parameters
      proposal_sd <- list()
      n_proposal_sd <- 5
      for (i in 1:chains) {
        proposal_sd[[i]] <- replicate(rungs, c(lambda = 1,
                                               theta = 1,
                                               decay_rate = 1,
                                               sens = 1,
                                               infection_time = 1),
                                      simplify = FALSE)
      }
      
      # initialise objects for storing results
      private$lambda_store <- replicate(chains, NULL)
      private$theta_store <- replicate(chains, NULL)
      private$decay_rate_store <- replicate(chains, NULL)
      private$sens_store <- replicate(chains, NULL)
      private$n_infections_store <- replicate(chains, NULL)
      private$infection_times_store <- replicate(chains, NULL)
      
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
      
      # haplo_freqs will have already been calculated when loading data, but
      # overwrite if defined manually here
      if (!is.null(haplo_freqs)) {
        private$haplo_freqs <- haplo_freqs
      }
      
      # store values
      private$start_time <- start_time
      private$end_time <- end_time
      private$chains <- chains
      private$rungs <- rungs
      private$max_infections <- max_infections
      private$param_list <- param_list
      private$param_update_list <- param_update_list
      private$proposal_sd <- proposal_sd
      private$n_proposal_sd <- n_proposal_sd
      private$rng_list <- rng_list
      private$iteration_counter <- iteration_counter
      private$acceptance_counter <- acceptance_counter
      private$swap_acceptance_counter <- swap_acceptance_counter
      private$duration <- duration
      
      # update program flow
      private$init_called <- TRUE
      private$burn_called <- FALSE
      private$sample_called <- FALSE
    },
    
    #--------------------
    #' @description
    #' Run burn-in MCMC
    #' 
    #' @param iterations number of burn-in iterations.
    #' @param target_acceptance all Metropolis-Hastings proposals will be
    #'   adaptively tuned to aim for this acceptance rate.
    #' @param silent if \code{TRUE} then suppress all console output. Default =
    #'   \code{FALSE}.
    burn = function(iterations, target_acceptance = 0.44, silent = FALSE) {
      
      # program flow
      assert_eq(private$init_called, TRUE, message = "MCMC must be initialised via init() before running burn-in phase")
      assert_eq(private$sample_called, FALSE, message = "Cannot return to burn-in phase after sample() has been called")
      
      # check inputs
      assert_single_pos_int(iterations)
      assert_greq(iterations, 10)
      assert_single_bounded(target_acceptance, inclusive_left = FALSE, inclusive_right = FALSE)
      assert_single_logical(silent)
      
      # update program flow
      private$burn_called <- TRUE
      
      # run main loop
      private$run_burn_sample(burnin = TRUE,
                              iterations = iterations,
                              target_acceptance = target_acceptance,
                              silent = silent)
    },
    
    #--------------------
    #' @description
    #' Run sampling MCMC
    #' 
    #' @param iterations number of sampling iterations.
    #' @param silent if \code{TRUE} then suppress all console output. Default =
    #'   \code{FALSE}.
    sample = function(iterations, silent = FALSE) {
      
      # program flow
      assert_eq(private$burn_called, TRUE, message = "MCMC must be burned in via burn() before running sampling phase")
      
      # check inputs
      assert_single_pos_int(iterations)
      assert_greq(iterations, 10)
      assert_single_logical(silent)
      
      # update program flow
      private$sample_called <- TRUE
      
      # run main loop
      private$run_burn_sample(burnin = FALSE,
                              iterations = iterations,
                              target_acceptance = 0.1, # this value is ignored
                              silent = silent)
    },
    
    #--------------------
    #' @description
    #' Get global parameter output
    get_output_global = function() {
      
      ret <- mapply(function(i) {
        iteration_counter <- unlist(private$iteration_counter[[i]][c("burn", "sample")])
        
        v_chain <- rep(i, sum(iteration_counter))
        v_phase <- rep(c("burnin", "sampling"), times = iteration_counter)
        v_iteration <- 1:sum(iteration_counter)
        
        data.frame(chain = v_chain,
                   phase = v_phase,
                   iteration = v_iteration,
                   lambda = private$lambda_store[[i]],
                   theta = private$theta_store[[i]],
                   decay_rate = private$decay_rate_store[[i]],
                   sens = private$sens_store[[i]])
      }, seq_len(private$chains), SIMPLIFY = FALSE) |>
        bind_rows() |>
        mutate(phase = factor(phase, levels = c("burnin", "sampling")),
               chain = factor(chain, levels = 1:private$chains))
      
      return(ret)
    },
    
    #--------------------
    #' @description
    #' Get number of infections output.
    get_output_n_infections = function() {
      
      ret <- mapply(function(i) {
        iteration_counter <- unlist(private$iteration_counter[[i]][c("burn", "sample")])
        
        v_chain <- rep(i, length(private$n_infections_store[[i]]))
        v_phase <- rep(c("burnin", "sampling"), times = iteration_counter * private$n_samp)
        v_iteration <- rep(1:sum(iteration_counter), each = private$n_samp)
        v_ind <- rep(1:private$n_samp, times = sum(iteration_counter))
        
        data.frame(chain = v_chain,
                   phase = v_phase,
                   iteration = v_iteration,
                   ind = v_ind,
                   value = as.vector(t(private$n_infections_store[[i]])))
      }, seq_len(private$chains), SIMPLIFY = FALSE) |>
        bind_rows() |>
        mutate(phase = factor(phase, levels = c("burnin", "sampling")),
               chain = factor(chain, levels = 1:private$chains))
      
      return(ret)
    },
    
    #--------------------
    #' @description
    #' Get infection times output.
    get_output_infection_times = function() {
      
      seq_len_V <- Vectorize(function(n) seq_len(n), SIMPLIFY = FALSE)
      
      ret <- mapply(function(i) {
        iteration_counter <- unlist(private$iteration_counter[[i]][c("burn", "sample")])
        
        n_infections <- as.vector(t(private$n_infections_store[[i]]))
        
        v_iter <- rep(1:sum(iteration_counter), times = rowSums(private$n_infections_store[[i]]))
        v_phase <- ifelse(v_iter <= iteration_counter[1], "burnin", "sampling")
        v_ind <- rep(rep(1:private$n_samp, times = sum(iteration_counter)), times = n_infections)
        v_inf <- unlist(seq_len_V(n_infections))
        
        # sort by infection time before indexing infections
        data.frame(chain = i,
                   phase = v_phase,
                   iteration = v_iter,
                   individual = v_ind,
                   value = unlist(private$infection_times_store[[i]])) |>
          arrange(chain, phase, iteration, individual, value) |>
          add_column(infection = v_inf, .after = "individual")
        
      }, seq_len(private$chains), SIMPLIFY = FALSE) |>
        bind_rows() |>
        mutate(phase = factor(phase, levels = c("burnin", "sampling")),
               chain = factor(chain, levels = 1:private$chains),
               infection = factor(infection, levels = 1:private$max_infections))
      
      return(ret)
    },
    
    #--------------------
    #' @description
    #' Print MCMC object summary
    print = function() {
      message("Data")
      message("--------")
      
      # print summary of data
      message(sprintf("%s samples", private$n_samp))
      message(sprintf("%s haplotypes", private$n_haplos))
      if (length(unique(private$obs_time_list)) == 1) {
        message(sprintf("%s observation times, identical for all samples", length(private$obs_time_list[[1]])))
      } else {
        message("variable observation times between samples")
      }
      message("")
      
      if (private$burn_called) {
        message("MCMC")
        message("--------")
        
        burnin_iterations <- sum(mapply(function(x) x$burn, private$iteration_counter))
        message(sprintf("Burn-in: %s iterations", burnin_iterations))
        
        sampling_iterations <- sum(mapply(function(x) x$sample, private$iteration_counter))
        message(sprintf("Sampling: %s iterations", sampling_iterations))
      }
      
      # return invisibly
      invisible(self)
    }
    
  ),
  
  private = list(
    
    ### Private variables ###
    
    # data
    df_data = NULL,
    list_data = NULL,
    n_samp = NULL,
    n_haplos = NULL,
    haplo_freqs = NULL,
    obs_time_list = NULL,
    
    # program flow
    init_called = FALSE,
    burn_called = FALSE,
    sample_called = FALSE,
    
    # MCMC parameters
    chains = NULL,
    rungs = NULL,
    max_infections = NULL,
    param_list = NULL,
    param_update_list = NULL,
    proposal_sd = NULL,
    n_proposal_sd = NULL,
    rng_list = NULL,
    iteration_counter = NULL,
    acceptance_counter = NULL,
    swap_acceptance_counter = NULL,
    duration = NULL,
    
    # misc parameters
    start_time = NULL,
    end_time = NULL,
    
    # parameter draws
    lambda_store = NULL,
    theta_store = NULL,
    decay_rate_store = NULL,
    sens_store = NULL,
    n_infections_store = NULL,
    infection_times_store = NULL,
    
    ### Private functions ###
    
    #--------------------
    run_burn_sample = function(burnin, iterations, target_acceptance = 0.44, silent = FALSE) {
      
      # check inputs
      assert_single_logical(burnin)
      assert_single_pos_int(iterations)
      assert_greq(iterations, 10)
      assert_single_bounded(target_acceptance, inclusive_left = FALSE, inclusive_right = FALSE)
      assert_single_logical(silent)
      
      phase_name <- ifelse(burnin, "burn", "sample")
      
      # loop over chains
      chains <- private$chains
      for (chain in 1:chains) {
        
        if (!silent) {
          message(sprintf("\nRunning chain %s", chain))
        }
        
        # run this chain
        output_raw <- mcmc_cpp(private$list_data,                                 # data in list format
                               private$obs_time_list,                             # observation times
                               private$haplo_freqs,                               # haplo freqs
                               iterations,                                        # iterations
                               burnin,                                            # burnin
                               private$param_list[[chain]],                       # params
                               private$param_update_list,                         # which params to update
                               private$proposal_sd[[chain]],                      # proposal_sd
                               private$iteration_counter[[chain]][[phase_name]],  # iteration_counter_init
                               rep(1, private$rungs),                             # beta
                               private$start_time,                                # start_time
                               private$end_time,                                  # end_time
                               private$max_infections,                            # max infections
                               private$rng_list[[chain]],                         # rng_ptr
                               interactive()                                      # logical; if running interactively. Useful for supressing progress bars when running Rmarkdown
        )
        
        # sync RNG
        private$rng_list[[chain]]$sync()
        
        # convert raw output objects if needed
        acceptance_out <- matrix(unlist(output_raw$acceptance_out),
                                 nrow = length(output_raw$acceptance_out),
                                 byrow = TRUE)
        n_infections <- matrix(unlist(output_raw$n_infections),
                               nrow = length(output_raw$n_infections),
                               byrow = TRUE)
        
        # update params_list to final values of chain
        private$param_list[[chain]] <- output_raw$param_list_out
        
        # append parameter draws
        private$lambda_store[[chain]] <- c(private$lambda_store[[chain]], output_raw$lambda)
        private$theta_store[[chain]] <- c(private$theta_store[[chain]], output_raw$theta)
        private$decay_rate_store[[chain]] <- c(private$decay_rate_store[[chain]], output_raw$decay_rate)
        private$sens_store[[chain]] <- c(private$sens_store[[chain]], output_raw$sens)
        private$n_infections_store[[chain]] <- rbind(private$n_infections_store[[chain]], n_infections)
        private$infection_times_store[[chain]] <- c(private$infection_times_store[[chain]], output_raw$infection_times)
        
        # update counters
        private$iteration_counter[[chain]][[phase_name]] <- private$iteration_counter[[chain]][[phase_name]] + iterations
        private$acceptance_counter[[chain]][[phase_name]] <- private$acceptance_counter[[chain]][[phase_name]] + acceptance_out
        private$swap_acceptance_counter[[chain]][[phase_name]] <- private$swap_acceptance_counter[[chain]][[phase_name]] + output_raw$swap_acceptance_out
        private$duration[[chain]][[phase_name]] <- private$duration[[chain]][[phase_name]] + output_raw$dur
        
      } # end loop over chains
      
    },
    
    #--------------------
    debug_algo1 = function() {
      chain <- 1
      debug_algo1_cpp(private$list_data,                                 # data in list format
                      private$obs_time_list,                             # observation times
                      private$haplo_freqs,                               # haplo freqs
                      private$param_list[[chain]],                       # params
                      private$param_update_list,                         # which params to update
                      private$proposal_sd[[chain]],                      # proposal_sd
                      rep(1, private$rungs),                             # beta
                      private$start_time,                                # start_time
                      private$end_time,                                  # end_time
                      private$max_infections,                            # max infections
                      private$rng_list[[chain]])                         # rng_ptr
    }
  )
)
