# Plot posterior densities of parameter averages by group

Plot posterior densities of group summaries of individual parameters.
The respondents can be grouped by any categorical variable and the
function works whether the fitted model is of "MULTI"-type or not.

## Usage

``` r
plot_by_group(
  object,
  par = "abs_beta",
  group_id = NULL,
  ascending_means = TRUE,
  fill = "#2166AC",
  color = "#053061",
  alpha = 0.5,
  ncol = max(1, round(length(unique(group_id))/10))
)
```

## Arguments

- object:

  An instance of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md).

- par:

  Character: Name of the parameter to be plotted. One of the following:
  `"alpha"`, `"beta"`, `"abs_beta"`, `"lambda"`, or `"chi"`. Defaults to
  `"abs_beta"`, which means the absolute values of the draws for beta
  will be used. Further individual-level parameters like `"eta"` can be
  specified if these have been passed to
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) via the
  argument `extra_pars` when fitting the model. (Note that homoskedastic
  models have no `"eta"` parameters and "NF"-type models have no
  `"lambda"` or `"kappa"` parameters.)

- group_id:

  An optional vector that will be used to split the respondents into
  groups. The vector must either be as long as the number of rows in the
  original dataset, or as long as the number of respondents included in
  the analysis. If a `group_id` was previously supplied to
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md)
  or [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) and
  if no `group_id` is supplied here, the default is to use the existing
  `group_id`. If a `group_id` is supplied here, it will be used instead
  of any previously supplied vector. The `group_id` supplied here does
  not have to coincide with the `group_id` used to fit a "MULTI"-type
  model: Any vector that can be used to group the respondents is
  allowed.

- ascending_means:

  Logical: Should the groups be placed in ascending order based on their
  posterior means (`TRUE`) or should they be ordered based on their
  names (`FALSE`)? Defaults to `TRUE`.

- fill:

  Fill color. Passed on to
  [`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html).

- color:

  Color of outer lines. Passed on to
  [`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html).

- alpha:

  Number in \[0,1\]: Inverse level of transparency.

- ncol:

  Number of columns. The default uses a formula to have approximately
  ten subplots per column.

## Value

A `ggplot` object.
