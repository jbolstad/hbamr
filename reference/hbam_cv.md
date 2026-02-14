# Perform K-fold cross-validation

This function performs K-fold cross-validation for an HBAM or FBAM model
in order to estimate the expected log pointwise predictive density for a
new dataset (ELPD). Multiple chains for one or more folds can be run in
parallel using the `future` package.

## Usage

``` r
hbam_cv(
  self = NULL,
  stimuli = NULL,
  model = "HBAM",
  allow_miss = 0,
  req_valid = NA,
  req_unique = 2,
  prefs = NULL,
  group_id = NULL,
  prep_data = TRUE,
  data = NULL,
  K = 10,
  chains = 2,
  warmup = 1000,
  iter = 3000,
  seed = 1,
  control = list(max_treedepth = 7),
  sigma_alpha = NULL,
  sigma_beta = 0.35,
  sigma_mu_alpha = NULL,
  sigma_mu_beta = 0.3,
  ...
)
```

## Arguments

- self:

  A numerical vector of N ideological self-placements. Any missing data
  must be coded as NA. This argument will not be used if the data have
  been prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function.

- stimuli:

  An N × J matrix of numerical stimulus placements, where J is the
  number of stimuli. Any missing data must be coded as NA. This argument
  will not be used if the data have been prepared in advance via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function.

- model:

  Character: Name of the model to be used. Defaults to HBAM.

- allow_miss:

  Integer specifying how many missing stimulus positions should be
  accepted for an individual still to be included in the analysis. This
  argument will not be used if the data have been prepared in advance
  via the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  function. Defaults to 0.

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

- prefs:

  An N × J matrix of numerical stimulus ratings or preference scores.
  These data are only required by the HBAM_R_MINI model and will be
  ignored when fitting other models.

- group_id:

  Integer vector of length N identifying which group each respondent
  belongs to. The supplied vector should range from 1 to the total
  number of groups in the data, and all integers between these numbers
  should be represented in the supplied data. These data are only
  required by models with "MULTI" in their name and will be ignored when
  fitting other models.

- prep_data:

  Logical: Should the data be prepared before fitting the model? (Or
  have the data been prepared in advance by first running the
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  and
  [`prep_data_cv()`](https://jbolstad.github.io/hbamr/reference/prep_data_cv.md)
  functions)? If so, set `prep_data = FALSE`.) Defaults to
  `prep_data = TRUE`.

- data:

  A list of data produced by
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  followed by
  [`prep_data_cv()`](https://jbolstad.github.io/hbamr/reference/prep_data_cv.md).

- K:

  An integer above 2, specifying the number of folds to use in the
  analysis. Defaults to 10.

- chains:

  A positive integer specifying the number of Markov chains to use per
  fold. Defaults to 2.

- warmup:

  A positive integer specifying the number of warmup (aka burn-in)
  iterations per chain. It defaults to 1000. The number of warmup
  iterations should be smaller than `iter`.

- iter:

  A positive integer specifying the number of iterations for each chain
  (including warmup). It defaults to 3000 as running fewer chains for
  longer is a more efficient way to obtain a certain number of draws
  (and cross-validation can be computationally expensive).

- seed:

  An integer passed on to `set.seed` before creating the folds to
  increase reproducibility and comparability. Defaults to 1 and only
  applies to fold-creation when the argument `prep_data` is `TRUE`. The
  supplied `seed` argument is also used to generate seeds for the
  sampling algorithm.

- control:

  A named list of parameters to control the sampler's behavior. See the
  documentation for
  [`rstan::stan`](https://mc-stan.org/rstan/reference/stan.html) for
  more details.

- sigma_alpha:

  A positive numeric value specifying the standard deviation of the
  prior on the shift parameters in the FBAM model, or the standard
  deviation of the parameters' deviation from the group-means in
  FBAM_MULTI models. (This argument will be ignored by HBAM models.)
  Defaults to B / 4, where B measures the length of the survey scale as
  the number of possible placements on one side of the center.

- sigma_beta:

  A positive numeric value specifying the standard deviation of the
  prior on the logged stretch parameters in the FBAM model, or the
  standard deviation of the logged parameters' deviation from the
  group-means in FBAM_MULTI models. (This argument will be ignored by
  HBAM models.) Defaults to .35.

- sigma_mu_alpha:

  A positive numeric value specifying the standard deviation of the
  prior on the group-means of the shift parameters in MULTI-type models.
  Defaults to B / 5.

- sigma_mu_beta:

  A positive numeric value specifying the standard deviation of the
  prior on the group-means of the logged stretch parameters in
  MULTI-type models. Defaults to .3.

- ...:

  Arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

A list of classes `kfold` and `loo`, which contains the following named
elements:

- `"estimates"`: A `1x2` matrix containing the ELPD estimate and its
  standard error. The columns have names `"Estimate"` and `"SE"`.

- `"pointwise"`: A `Nx1` matrix with column name `"elpd_kfold"`
  containing the pointwise contributions for each data point.

## Examples

``` r
# \donttest{
# Loading and re-coding ANES 1980 data:
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA

# Making a small subset of the data for illustration:
self <- LC1980[1:50, 1]
stimuli <- LC1980[1:50, -1]

# Preparing to run chains in parallel using 2 cores via the future package:
  # Note: You would normally want to use all physical cores for this.
future::plan(future::multisession, workers = 2)

# Performing 4-fold cross-validation for the HBAM_MINI model:
  # Note: You would typically want to run the chains for more iterations.
cv_hbam_mini <- hbam_cv(self, stimuli, model = "HBAM_MINI", K = 4,
                        chains = 1, warmup = 500, iter = 1000)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 42 (50)
#> - Number of stimuli: 6 (6)
#> - Number of stimuli obs.: 252 (289)
#> - Range of observations: [-3, 3] ([1, 7])
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# Performing 4-fold cross-validation for the FBAM model:
cv_fbam <- hbam_cv(self, stimuli, model = "FBAM", K = 4,
                        chains = 1, warmup = 500, iter = 1000)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 42 (50)
#> - Number of stimuli: 6 (6)
#> - Number of stimuli obs.: 252 (289)
#> - Range of observations: [-3, 3] ([1, 7])

# Comparing the results using the loo package:
loo::loo_compare(list(HBAM_MINI = cv_hbam_mini,
                 FBAM = cv_fbam))
#>           elpd_diff se_diff
#> HBAM_MINI  0.0       0.0   
#> FBAM      -8.3       7.9   

# Stop the cluster of parallel sessions:
future::plan(future::sequential)
# }
```
