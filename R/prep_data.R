#' Prepare data to fit an HBAM model
#'
#' This function prepares data to fit a hierarchical Bayesian Aldrich-McKelvey (HBAM) model. It can be run ahead of fitting the models, or it can be run implicitly as part of a single function call to fit the models using `hbam`. It applies a set of inclusion criteria, performs any necessary data transformation, and returns a list of data suited for sampling in `rstan`. The data provided to `prep_data` can be centered, but they do not have to be: The function will detect un-centered data and attempt to center these automatically, assuming that the highest and lowest observed values in the data mark the extremes of the scale.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA.
#' @param prefs An N × J matrix of numerical stimulus ratings or preference scores. These data are only required by the `"HBAM_R"` and `"HBAM_R_MINI"` models and will be ignored when fitting other models.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`.
#' @param B Integer specifying the upper bound of the survey scale after centering. If not supplied, this information will be inferred from the data.
#' @return A list of data to be used by `hbam`. The returned list includes the logical vector `keep`, which identifies the rows in the original data that have been kept for further analysis. The stimuli data are stored in a vector as a long-form sparse matrix. If the stimuli data include column-names, these will be preserved for later use.
#' @examples
#' # Loading and re-coding ANES 1980 data:
#' data(LC1980)
#' dat <- LC1980
#' dat[dat == 0 | dat == 8 | dat == 9] <- NA
#' self <- dat[, 1]
#' stimuli <- dat[, -1]
#'
#' # Prepare data for model fitting, using defaults:
#' dat <- prep_data(self, stimuli)
#'
#' # Prepare data for model fitting, using using alternative settings:
#' dat2 <- prep_data(self, stimuli, allow_miss = 0, req_unique = 3)
#'
#' # Obtain the data that are included in the analysis:
#' self2 <- self[dat2$keep]
#' stimuli2 <- stimuli[dat2$keep, ]

prep_data <- function(self, stimuli,
                      prefs = NULL,
                      allow_miss = 2,
                      req_valid = NA, req_unique = 2, B = NULL) {

  #dimnames(stimuli) <- NULL
  #self <- as_numeric(self)

  # req_valid takes precedence over allow_miss they are not consistent:
  if (!is.na(req_valid)) { allow_miss <- ncol(stimuli) - req_valid }

  all_dat <- cbind(stimuli, self)

  # Counting number of unique stimuli positions reported per respondent:
  n_unique <- apply(stimuli, 1, function(x) {length(unique(x[!is.na(x)]))} )

  # Center position-data if necessary:
  minval <- min(all_dat, na.rm = T)
  if (minval >= 0) {
    maxval <- max(all_dat, na.rm = T)
    stimuli <- (stimuli - minval) - ((maxval - minval) / 2)
    self <- (self - minval) - ((maxval - minval) / 2)
  }

  if (!is.null(prefs)) {
    prefs <- prefs / max(prefs, na.rm = T)

    # Rescale preferences to range from 0 to 1 per for each individual:
    suppressWarnings(pref_min <- apply(prefs, 1, min, na.rm = T))
    suppressWarnings(pref_max <- apply(prefs, 1, max, na.rm = T))
    sel <- pref_min != pref_max
    prefs[sel,] <- (prefs[sel,] - pref_min[sel]) / (pref_max[sel] - pref_min[sel])

    # Indicators of non-constant preference data for each respondent:
    has_var <- apply(prefs, 1, var, na.rm = T) > 0

    # Which observations are missing either prefs or positions:
    is_any_na <- apply((is.na(prefs) | is.na(stimuli)), 2, as.numeric)

    keep <- !is.na(self) & apply(is_any_na, 1, sum) <= allow_miss & n_unique >= req_unique & has_var
    prefs <- prefs[keep,]
  } else {
    keep <- !is.na(self) & apply(is.na(stimuli), 1, sum) <= allow_miss & n_unique >= req_unique
  }

  stimuli <- stimuli[keep, ]
  self <- self[keep]
  mean_spos <- apply(stimuli, 2, mean, na.rm = T)

  # Coding Y as a long-form sparse matrix:
  stimuli_vec <- as.numeric(as.matrix(stimuli))
  if (!is.null(prefs)) {
    prefs_vec <- as.numeric(as.matrix(prefs))
    drop <- is.na(stimuli_vec) | is.na(prefs_vec)
    prefs_vec <- prefs_vec[!drop]
  } else{
    drop <- is.na(stimuli_vec)
    prefs_vec <- 0
  }
  ii <- rep(1:nrow(stimuli), ncol(stimuli))[!drop]
  jj <- rep(1:ncol(stimuli), each = nrow(stimuli))[!drop]
  stimuli_vec <- stimuli_vec[!drop]

  if(is.null(B)) { B <- max(abs(c(stimuli_vec, self)), na.rm = T) }

  return(list(J = ncol(stimuli), N = nrow(stimuli), B = B, N_obs = length(stimuli_vec),
              V = self, Y = stimuli_vec, U = prefs_vec, L = which.min(mean_spos), R = which.max(mean_spos),
              ii = ii, jj = jj, mean_spos = mean_spos, keep = keep, names = colnames(stimuli),
              CV = 0, holdout = rep(0, length(stimuli_vec))))
}
