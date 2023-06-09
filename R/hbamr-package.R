#' Hierarchical Bayesian Aldrich-McKelvey Scaling via Stan
#'
#' @description Fit hierarchical Bayesian Aldrich-McKelvey (HBAM) models using a form of Hamiltonian Monte Carlo via Stan. Aldrich-McKelvey (AM) scaling is a method for estimating the ideological positions of survey respondents and political actors on a common scale using positional survey data. The hierarchical versions of the Bayesian AM model included in this package outperform other versions by a considerable margin both in terms of yielding meaningful posterior distributions for respondent positions and in terms of recovering true respondent positions in simulations. The package contains functions for preparing data, fitting models, extracting estimates, plotting key results, and comparing models using cross-validation.
#'
#' @author Jørgen Bølstad
#'
#' @seealso \itemize{\item \url{https://github.com/jbolstad/hbamr/}}
#'
#' @docType package
#' @name hbamr-package
#' @aliases hbamr
#' @useDynLib hbamr, .registration = TRUE
#' @import methods stats
#' @import Rcpp ggplot2 RColorBrewer parallel tidyr loo
#' @importFrom rstan sampling sflist2stanfit
#' @importFrom matrixStats colLogSumExps
#' @importFrom plyr llply
#' @importFrom dplyr bind_cols
#' @importFrom pbmcapply pbmclapply
#' @importFrom rlang .data
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstantools rstan_config
#'
#' @references
#' - Bølstad, Jørgen (2023). Hierarchical Bayesian Aldrich-McKelvey Scaling. \emph{Political Analysis}. \url{https://doi.org/10.1017/pan.2023.18}.
#' - Stan Development Team (2023). RStan: the R interface to Stan. R package version 2.21.8. \url{https://mc-stan.org}.
#'
NULL
