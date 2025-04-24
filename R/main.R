# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("stats", "end", "start", "time", "haplo", "positive", "t_inf",
                           "t_start", "t_clear", "state_true", "state_obs", "state_combined"))
}
