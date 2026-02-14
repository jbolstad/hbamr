# Plot individual parameter estimates over self-placements

Create a boxplot of individual parameter point estimates from an HBAM
model over self-placements

## Usage

``` r
plot_over_self(
  object,
  par = "chi",
  estimate = "median",
  names = NULL,
  parlabel = NULL,
  fill = "#2166AC",
  color = "#053061",
  width = 0.7,
  alpha = 0.5,
  outlier.size = 0.3,
  median_color = "black",
  median_lwd = 0.7
)
```

## Arguments

- object:

  An object of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md), a list
  produced by
  [`fbam()`](https://jbolstad.github.io/hbamr/reference/fbam.md), or a
  `list` of such objects, which will produce a faceted plot.

- par:

  Character: Name of the parameter to be plotted. One of the following:
  `"alpha"`, `"beta"`, `"abs_beta"`, `"lambda"`, or `"chi"`. Defaults to
  `"chi"`. Further individual-level parameters like `"eta"` can be
  specified if these have been passed to
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) via the
  argument `extra_pars` when fitting the model. (Note that homoskedastic
  models have no `"eta"` parameters and "NF"-type models have no
  `"lambda"` or `"kappa"` parameters.)

- estimate:

  Character: Specifying which type of posterior point estimate to use.
  One of `"median"` and `"mean"`. Defaults to `"median"`. This only
  applies for `stanfit` objects.

- names:

  An optional character vector of model names of same length as the
  supplied list of models.

- parlabel:

  An optional character containing an alternative label for the
  parameter (will be parsed if passed as an expression).

- fill:

  Fill color of boxes. Passed on to
  [`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- color:

  Color of outer lines. Passed on to
  [`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- width:

  Width of boxes. Passed on to
  [`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- alpha:

  Number in \[0,1\]: Inverse level of transparency for fill color.

- outlier.size:

  Size of dots representing outliers. Passed on to
  [`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- median_color:

  Color of solid line representing the median.

- median_lwd:

  Thickness of solid line representing the median.

## Value

A `ggplot` object.
