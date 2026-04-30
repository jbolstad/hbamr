# Fit an FBAM model using optimization

Fit a simplified Bayesian Aldrich-McKelvey model with fixed
hyperparameters using optimization via `rstan`. Users may replace the
default priors by supplying their own values for the hyperparameters.

## Usage

``` r
fbam(
  self = NULL,
  stimuli = NULL,
  model = "FBAM",
  allow_miss = 2,
  req_valid = NA,
  req_unique = 2,
  group_id = NULL,
  data = NULL,
  seed = sample.int(.Machine$integer.max, 1),
  sigma_alpha = NULL,
  sigma_beta = 0.3,
  sigma_mu_alpha = NULL,
  sigma_mu_beta = 0.2,
  ...
)
```

## Arguments

- self:

  An optional numerical vector of N ideological self-placements. Any
  missing data must be coded as NA. If this argument is not supplied
  (either here or in a previous call to
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)),
  respondent positions will not be estimated. If the data have been
  prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function, the argument supplied here will be ignored.

- stimuli:

  An N × J matrix of numerical stimulus placements, where J is the
  number of stimuli. Any missing data must be coded as NA. This argument
  will not be used if the data have been prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function.

- model:

  Character: Name of the model to be used. Defaults to FBAM. The
  available options are the three models with "FBAM" in their name. See
  the documentation for the
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md)
  function for descriptions of the models.

- allow_miss:

  Integer specifying how many missing stimulus positions should be
  accepted for an individual still to be included in the analysis. This
  argument will not be used if the data have been prepared in advance
  via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function. Defaults to 2.

- req_valid:

  Integer specifying how many valid observations to require for a
  respondent to be included in the analysis. The default is
  `req_valid = J - allow_miss`, but if specified, `req_valid` takes
  precedence. This argument will not be used if the data have been
  prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function.

- req_unique:

  Integer specifying how many unique positions on the ideological scale
  each respondent is required to have used when placing the stimuli in
  order to be included in the analysis. The default is `req_unique = 2`.
  This argument will not be used if the data have been prepared in
  advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function.

- group_id:

  Vector of length N identifying which group each respondent belongs to.
  The format can be factor, character, integer, or numeric. Respondents
  with NAs on `group_id` will be dropped when `group_id` is supplied.
  These data are only required by models with `"MULTI"` in their name
  and will be ignored when fitting other models.

- data:

  List of data that have been prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function. Not required if the arguments `self` and `stimuli` are
  provided.

- seed:

  A positive integer specifying an optional seed for reproducibility. If
  this argument is not supplied, a random seed will be generated and the
  function will produce slightly different results on each run.

- sigma_alpha:

  A positive numeric value specifying the standard deviation of the
  prior on the shift parameters in the FBAM model, or the standard
  deviation of the parameters' deviation from the group-means in
  FBAM_MULTI models. (This argument will be ignored by HBAM models.)
  Defaults to B / 5, where B measures the length of the survey scale as
  the number of possible placements on one side of the center.

- sigma_beta:

  A positive numeric value specifying the standard deviation of the
  prior on the logged stretch parameters in the FBAM model, or the
  standard deviation of the logged parameters' deviation from the
  group-means in FBAM_MULTI models. (This argument will be ignored by
  HBAM models.) Defaults to .3.

- sigma_mu_alpha:

  A positive numeric value specifying the standard deviation of the
  prior on the group-means of the shift parameters in MULTI-type models.
  Defaults to B / 10.

- sigma_mu_beta:

  A positive numeric value specifying the standard deviation of the
  prior on the group-means of the logged stretch parameters in
  MULTI-type models. Defaults to .2.

- ...:

  Arguments passed to
  [`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html).

## Value

A list produced by
[`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html).

## Examples

``` r
# \donttest{
# Loading ANES 2012 data:
data(LC2012)

# Making a small subset of the data for illustration:
self <- LC2012[1:1000, 2]
stimuli <- LC2012[1:1000, -c(1:2)]

# Fitting the FBAM model:
fit_fbam <- fbam(self, stimuli)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 654 (1000)
#> - Number of stimuli: 4 (4)
#> - Number of stimuli obs.: 2559 (3564)
#> - Range of observations: [-3, 3] ([1, 7])

# Obtaining point estimates for the latent stimulus positions:
theta_est <- get_est(fit_fbam, par = "theta")
# }
```
