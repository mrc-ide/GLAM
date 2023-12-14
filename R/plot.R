#------------------------------------------------
#' @title Plot a single simulated individual
#'
#' @description Produces a plot showing the infections in a single individual,
#'   and which haplotypes were introduced in each infection. Requires that the
#'   individual was simulated with hidden states returned (\code{return_full =
#'   TRUE}).
#'
#' @param ind a single individual simulated using \code{sim_ind()} with
#'   \code{return_full = TRUE}.
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @export

plot_ind <- function(ind) {
  
  # avoid "no visible binding" notes
  start <- haplo <- end <- time <- NULL
  
  # check inputs
  assert_in(c("state_obs", "t_inf", "w_init", "clear_init", "w_inf", "clear_inf",
              "state_true", "samp_time"), names(ind))
  
  # extract values
  state_obs <- ind$state_obs
  t_inf <- ind$t_inf
  w_init <- ind$w_init
  clear_init <- ind$clear_init
  w_inf <- ind$w_inf
  clear_inf <- ind$clear_inf
  state_true <- ind$state_true
  samp_time <- ind$samp_time
  
  n_haplo <- length(w_init)
  n_samp <- length(samp_time)
  start_time <- samp_time[1]
  end_time <- samp_time[n_samp]
  n_inf <- length(t_inf)
  
  # plot infection times
  df_plot <- data.frame(start = rep(start_time, sum(w_init)),
                        end = clear_init[w_init == 1],
                        haplo = which(w_init == 1))
  if (n_inf > 0) {
    df_plot <- df_plot |>
      bind_rows(data.frame(start = (w_inf*t_inf)[w_inf == 1],
                           end = clear_inf[w_inf == 1],
                           haplo = which(w_inf == 1, arr.ind = TRUE)[,2]))
  }
  
  plot1 <- ggplot(df_plot) + theme_bw() +
    geom_vline(xintercept = t_inf, linetype = "dashed") +
    geom_point(aes(x = start, y = haplo), size = 3, pch = 4, stroke = 1.3) +
    geom_segment(aes(x = start, y = haplo, xend = end, yend = haplo)) +
    xlab("Time") + ylab("Haplotype") +
    xlim(c(start_time, end_time)) + ylim(c(1, n_haplo)) +
    theme(panel.grid = element_blank())
  
  # add haplotypes to plot
  df_plot <- data.frame(time = samp_time,
                        haplo = rep(1:n_haplo, each = n_samp),
                        state_true = c("Absent", "Present")[state_true + 1],
                        state_obs = c("Unobserved", "Observed")[state_obs + 1]) |>
    mutate(state_true = factor(state_true, levels = c("Absent", "Present")),
           state_obs = factor(state_obs, levels = c("Unobserved", "Observed")))
  
  plot2 <- plot1 + 
    geom_point(aes(x = time, y = haplo, alpha = state_true, col = state_obs),
               size = 3, data = df_plot) +
    scale_color_discrete(name = "Observed status") +
    scale_alpha_manual(values = c(0.1, 1), name = "True status")
  
  # tweak legends etc
  plot2
}
