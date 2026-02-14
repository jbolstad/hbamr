# Extract point estimates or other summaries of marginal posterior distributions

For objects produced by
[`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md), this
function is a wrapper for `rstan::summary()`. For objects produced by
[`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md) it offers
a way to extract point estimates.

## Usage

``` r
get_est(
  object,
  par = "theta",
  format_orig = FALSE,
  probs = c(0.025, 0.5, 0.975),
  simplify = TRUE,
  ...
)
```

## Arguments

- object:

  An instance of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md), or a
  list produced by
  [`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md).

- par:

  Character: Name of the parameter type to be extracted. Typically
  `"theta"` (stimuli positions) or `"chi"` (respondent positions).

- format_orig:

  Logical: Should individual-level parameters be mapped to the original
  dataset by returning rows of NAs for respondents who were not included
  in the analysis? Defaults to `FALSE`.

- probs:

  A numeric vector of quantiles of interest for summarizing `stanfit`
  objects.

- simplify:

  Logical: Should the returned object be simplified by dropping the
  Monte Carlo standard error and the posterior standard deviation?
  Defaults to `TRUE`.

- ...:

  Other arguments are passed on to `rstan::summary()` when summarizing
  `stanfit` objects.

## Value

A tibble containing summaries of marginal posterior distributions. For
objects produced by
[`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md), only
maximum a posteriori estimates are returned.
