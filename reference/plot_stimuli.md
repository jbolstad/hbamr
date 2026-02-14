# Plot estimated stimulus positions

Plot marginal posterior distributions of stimulus positions from an HBAM
model

## Usage

``` r
plot_stimuli(object, rev_color = FALSE, alpha = 0.55)
```

## Arguments

- object:

  An instance of class `stanfit` produced by
  [`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md).

- rev_color:

  Logical: Display low positions as red and high positions as blue.

- alpha:

  Number in \[0,1\]: Inverse level of transparency for fill color.

## Value

A `ggplot` object.
