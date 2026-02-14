# Plot estimated respondent positions

Plot the distribution of estimated respondent positions from an HBAM
model.

## Usage

``` r
plot_respondents(
  object,
  inc_stimuli = TRUE,
  n_draws = 15,
  color = "#053061",
  fill = "#2166AC",
  alpha_color = 0.6,
  alpha_fill = 0.7/n_draws,
  seed = 1
)
```

## Arguments

- object:

  An instance of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md).

- inc_stimuli:

  Logical: Should estimated stimulus positions also be shown?

- n_draws:

  Integer specifying the number of posterior draws to use when
  illustrating the uncertainty of the population distribution. Defaults
  to 15.

- color:

  Color of lines illustrating uncertainty.

- fill:

  Fill color for density plots.

- alpha_color:

  Number in \[0,1\]: Inverse level of transparency for line color.

- alpha_fill:

  Number in \[0,1\]: Inverse level of transparency for fill color.

- seed:

  A positive integer specifying an optional seed for reproducibility.
  The seed is used to select respondent position draws for illustrating
  uncertainty.

## Value

A `ggplot` object.
