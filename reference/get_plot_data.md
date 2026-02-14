# Extract data for plotting results from an HBAM model

Extract data for plotting results from an HBAM model.

## Usage

``` r
get_plot_data(object, n_draws = 15, seed = 1)
```

## Arguments

- object:

  An instance of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) or a
  list produced by
  [`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md).

- n_draws:

  Integer specifying the number of posterior draws to use when
  illustrating the uncertainty of the population distribution. This only
  applies for `stanfit` objects.

- seed:

  A positive integer specifying an optional seed for reproducibility.
  The seed is used to select respondent position draws for illustrating
  uncertainty. This only applies for `stanfit` objects.

## Value

A list of three `tibble`s: The first element contains the posterior mean
stimulus positions, as well as the x- and y-values of the posterior
modes (which can be useful for labeling the distributions). The second
element contains the posterior draws for the stimulus positions (which
can be used to calculate marginal posterior densities). The third
element contains the selected number of posterior draws for each
respondent (which form the key ingredient for `plot_respondents`).
