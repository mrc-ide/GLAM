#------------------------------------------------
#' @title Plot simulated infection trajectories
#'
#' @description
#' Given a \code{sim_list} produced by your simulation routines, \code{plot_sim()}
#' will draw, for each individual:
#' * Vertical dashed lines at each true infection time.
#' * Pink segments showing any initial infections.
#' * Black segments for new infections occurring after the start.
#' * Points coloured by observed vs.\ unobserved haplotype states.
#'
#' @param sim_list A list with at least two elements:
#'   \describe{
#'     \item{\code{raw_list}}{A list of simulation output objects, each with
#'       components \code{t_inf}, \code{clear_init}, \code{w_init},
#'       \code{clear_inf}, \code{w_inf}, \code{samp_time},
#'       \code{state_true}, and \code{state_obs}.}
#'     \item{\code{df_data}}{A \code{data.frame} with columns
#'       \code{haplo} and \code{time}, used to determine the number of
#'       haplotypes (\code{n_haplo}) and observation times (\code{n_obs}).}
#'   }
#' @param ind Integer vector or \code{NULL}. Which indices of
#'   \code{sim_list$raw_list} to plot. Defaults to all individuals.
#' @param nrow passed through to \code{ggplot2::facet_wrap()} to set the number
#'   of rows when plotting multiple individuals.
#'
#' @return A \pkg{ggplot2} object with:
#'   \itemize{
#'     \item Infection times as vertical dashed lines.
#'     \item Initial infections (pink segments).
#'     \item Subsequent infections (black segments).
#'     \item Observed haplotype states as coloured points.
#'   }
#'   
#' @import ggplot2
#' @importFrom dplyr mutate
#' @export

plot_sim <- function(sim_list, ind = NULL, nrow = NULL) {
  
  # subset list
  if (is.null(ind)) {
    ind <- 1:length(sim_list$raw_list)
  }
  sim_list$raw_list <- sim_list$raw_list[ind]
  
  # get basic properties
  n_ind <- length(sim_list$raw_list)
  n_haplo <- length(unique(sim_list$df_data$haplo))
  n_obs <- length(unique(sim_list$df_data$time))
  start_time <- mapply(function(x) {
    min(x$samp_time)
  }, sim_list$raw_list) |>
    min()
  end_time <- mapply(function(x) {
    max(x$samp_time)
  }, sim_list$raw_list) |>
    max()
  
  # plot infection times
  df_t_inf <- mapply(function(x, i) {
    if (length(x$t_inf) == 0) {
      return(NULL)
    }
    data.frame(ind = i,
               t_inf = x$t_inf)
  }, sim_list$raw_list, 1:n_ind, SIMPLIFY = FALSE) |>
    bind_rows()
  
  plot1 <- ggplot(df_t_inf) + theme_bw() +
    geom_vline(aes(xintercept = t_inf), linetype = "dashed") +
    facet_wrap(~ind, nrow = nrow) +
    scale_x_continuous(limits = c(start_time, end_time))
  
  
  # overlay starting infection segments
  df_init_inf <- mapply(function(x, i) {
    if (any(x$w_init)) {
      w <- which(x$w_init)
      data.frame(ind = i,
                 haplo = w,
                 t_start = start_time,
                 t_clear = x$clear_init[w])
    } else {
      NULL
    }
  }, sim_list$raw_list, 1:n_ind, SIMPLIFY = FALSE) |>
    bind_rows()
  
  plot1 <- plot1 +
    geom_segment(aes(x = t_start, xend = t_clear, y = haplo), col = "pink", data = df_init_inf)
  
  
  # overlay new infection segments
  df_new_inf <- mapply(function(x, i) {
    if (any(x$w_inf)) {
      w <- which(x$w_inf, arr.ind = TRUE)
      data.frame(ind = i,
                 haplo = w[,2],
                 t_start = x$t_inf[w[,1]],
                 t_clear = x$clear_inf[w])
    } else {
      NULL
    }
  }, sim_list$raw_list, 1:n_ind, SIMPLIFY = FALSE) |>
    bind_rows()
  
  plot1 <- plot1 +
    geom_segment(aes(x = t_start, xend = t_clear, y = haplo), data = df_new_inf)
  
  
  # overlay observed haplos
  df_haplo <- mapply(function(x, i) {
    data.frame(ind = i,
               time = x$samp_time,
               haplo = rep(1:n_haplo, each = n_obs),
               state_true = c("Absent", "Present")[x$state_true + 1],
               state_obs = c("Unobserved", "Observed")[x$state_obs + 1]) |>
      mutate(state_combined = paste(state_true, state_obs, sep = "_"),
             state_combined = factor(state_combined, levels = c("Absent_Unobserved",
                                                                "Present_Unobserved",
                                                                "Present_Observed")))
  }, sim_list$raw_list, 1:n_ind, SIMPLIFY = FALSE) |>
    bind_rows()
  
  plot1 <- plot1 + 
    geom_point(aes(x = time, y = haplo, col = state_combined), size = 2, data = df_haplo) +
    scale_color_manual(values = c("#00000010", "firebrick1", 'dodgerblue'),
                       labels = c("Absent, unobserved",
                                  "Present, unobserved",
                                  "Present, observed"), 
                       name = "Observed status")
  
  plot1
}

#------------------------------------------------
#' @title Plot observed haplotypes
#'
#' @description Produces a plot showing the haplotypes observed at each
#'   observation time. If the input data spans multiple individuals then these
#'   will be faceted.
#'
#' @param df_data a data.frame with columns for \code{ind} (factor),
#'   \code{haplo} (factor), \code{time} (numeric), and \code{positive}
#'   (logical).
#' @param nrow passed through to \code{ggplot2::facet_wrap()} to set the number
#'   of rows when plotting multiple individuals.
#'
#' @import ggplot2
#' @export

plot_data <- function(df_data, nrow = NULL) {
  
  # check inputs
  assert_dataframe(df_data)
  assert_in(c("ind", "haplo", "time", "positive"), names(df_data))
  assert_int(df_data$time)
  assert_in(df_data$positive, c(0, 1))
  
  df_data |>
    ggplot() + theme_bw() +
    geom_point(aes(x = time, y = haplo, alpha = positive), size = 1) +
    facet_wrap(~ind, nrow = nrow) +
    theme(panel.grid = element_blank()) +
    guides(alpha = "none") +
    xlab("Time") + ylab("Haplotype")
}
