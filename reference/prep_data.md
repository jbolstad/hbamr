# Prepare data to fit an HBAM or FBAM model

This function prepares data to fit a hierarchical Bayesian
Aldrich-McKelvey (HBAM) model. It can be run ahead of fitting the
models, or it can be run implicitly as part of a single function call to
fit the models using
[`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) or
[`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md). It
applies a set of inclusion criteria, performs any necessary data
transformation, and returns a list of data suited for sampling in
`rstan`. The data provided to `prep_data()` can be centered, but they do
not have to be: The function will detect uncentered data and attempt to
center these automatically, assuming that the highest and lowest
observed values in the data mark the extremes of the scale.

## Usage

``` r
prep_data(
  self = NULL,
  stimuli,
  prefs = NULL,
  allow_miss = 2,
  req_valid = NA,
  req_unique = 2,
  B = NULL,
  group_id = NULL,
  quiet = FALSE
)
```

## Arguments

- self:

  An optional numerical vector of N ideological self-placements. Any
  missing data must be coded as NA. If this argument is not supplied,
  respondent positions will not be estimated.

- stimuli:

  An N × J matrix of numerical stimulus placements, where J is the
  number of stimuli. Any missing data must be coded as NA.

- prefs:

  An N × J matrix of numerical stimulus ratings or preference scores.
  These data are only required by the `"HBAM_R_MINI"` model and will be
  ignored when fitting other models.

- allow_miss:

  Integer specifying how many missing stimulus positions should be
  accepted for an individual still to be included in the analysis.
  Defaults to 2.

- req_valid:

  Integer specifying how many valid observations to require for a
  respondent to be included in the analysis. The default is
  `req_valid = J - allow_miss`, but if specified, `req_valid` takes
  precedence.

- req_unique:

  Integer specifying how many unique positions on the ideological scale
  each respondent is required to have used when placing the stimuli in
  order to be included in the analysis. The default is `req_unique = 2`.

- B:

  Scalar specifying the upper bound of the survey scale after centering.
  If not supplied, this information will be inferred from the data.

- group_id:

  Vector of length N identifying which group each respondent belongs to.
  The format can be factor, character, integer, or numeric. Respondents
  with NAs on `group_id` will be dropped when `group_id` is supplied.
  These data are only required by models with `"MULTI"` in their name
  and will be ignored when fitting other models.

- quiet:

  Logical: Should information about the data be printed to the console?
  Defaults to `FALSE`.

## Value

A list of data to be used by
[`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) or
[`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md). The
returned list includes the logical vector `keep`, which identifies the
rows in the original data that have been kept for further analysis. The
stimuli data are stored in a vector as a long-form sparse matrix. If the
stimuli data include column-names, these will be preserved for later
use.

## Examples

``` r
# Loading and re-coding ANES 1980 data:
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA
self <- LC1980[, 1]
stimuli <- LC1980[, -1]

# Prepare data for model fitting, using defaults:
dat <- prep_data(self, stimuli)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 881 (888)
#> - Number of stimuli: 6 (6)
#> - Number of stimuli obs.: 4937 (4973)
#> - Range of observations: [-3, 3] ([1, 7])

# Prepare data for model fitting, using alternative settings:
dat2 <- prep_data(self, stimuli, allow_miss = 0, req_unique = 3)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 603 (888)
#> - Number of stimuli: 6 (6)
#> - Number of stimuli obs.: 3618 (4973)
#> - Range of observations: [-3, 3] ([1, 7])

# Obtain the data that are included in the analysis:
self2 <- self[dat2$keep]
stimuli2 <- stimuli[dat2$keep, ]
```
