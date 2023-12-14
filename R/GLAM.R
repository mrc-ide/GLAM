#------------------------------------------------
#' @title GLAM
#'
#' @description Malaria infections often contain multiple haplotypes. When
#'   sequenced at several points in time, these haplotypes may be present for a
#'   period beforing cleaing as parasite densities decrease. We may also fail to
#'   detect a haplotype at a given point in time due to imperfect sensitivity of
#'   sequencing. GLAM takes these factors into account in a probabilistic model
#'   and produces estimates of the time(s) at which each individual in a
#'   longitudinal cohort became infected, for how long this infection lasted,
#'   and various other parameters.
#'
#' @docType package
#' @name GLAM
NULL

#------------------------------------------------
# required for cpp11 package
#' @useDynLib GLAM, .registration = TRUE
NULL
