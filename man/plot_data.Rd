% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_data}
\alias{plot_data}
\title{Plot observed haplotypes}
\usage{
plot_data(df_data, nrow = NULL)
}
\arguments{
\item{df_data}{a data.frame with columns for \code{ind} (factor),
\code{haplo} (factor), \code{time} (numeric), and \code{positive}
(logical).}

\item{nrow}{passed through to \code{ggplot2::facet_wrap()} to set the number
of rows when plotting multiple individuals.}
}
\description{
Produces a plot showing the haplotypes observed at each
  observation time. If the input data spans multiple individuals then these
  will be faceted.
}
