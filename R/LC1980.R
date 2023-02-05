#' 1980 Liberal-Conservative Scales
#'
#' Liberal-Conservative 7-point scales from the 1980 National Election Study.
#' Includes (in order) self-placement, and rankings of Carter, Reagan, Kennedy,
#' Anderson, Republican party, Democratic Party. Stored as a matrix of integers.
#' The numbers 0, 8, and 9 are considered to be missing values.
#'
#' @docType data
#'
#' @usage data(LC1980)
#'
#' @keywords datasets
#'
#' @source American National Election Studies: \href{https://www.electionstudies.org/}{www.electionstudies.org/}.
#'
#' This dataset was originally part of the `basicspace` package under the same name ("LC1980").
#'
#' @examples
#' data(LC1980)
#' dat <- LC1980
#' dat[dat == 0 | dat == 8 | dat == 9] <- NA
#' head(dat)
"LC1980"
