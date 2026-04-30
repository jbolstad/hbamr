# Prepare data for a K-fold cross-validation of an HBAM model

This function turns data prepared for
[`hbam()`](https://jbolstad.github.io/hbamr/reference/hbam.md) into a
list of K versions, where each version includes a different vector
identifying holdout data.

## Usage

``` r
prep_data_cv(data, K = 10, seed = 1)
```

## Arguments

- data:

  A list of data produced by
  [`prep_data()`](https://jbolstad.github.io/hbamr/reference/prep_data.md).

- K:

  An integer above 2, specifying the number of folds to use in the
  analysis. Defaults to 10.

- seed:

  An integer passed on to `set.seed` before creating the folds to
  increase reproducibility. Defaults to 1.

## Value

A list of K data objects where each version includes a different vector
identifying holdout data.

## Examples

``` r
# Loading and re-coding ANES 1980 data:
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA
self <- LC1980[, 1]
stimuli <- LC1980[, -1]
dat <- prep_data(self, stimuli)
#> Summary of prepared data (values for supplied data in parentheses)
#> - Number of respondents: 881 (888)
#> - Number of stimuli: 6 (6)
#> - Number of stimuli obs.: 4937 (4973)
#> - Range of observations: [-3, 3] ([1, 7])

# Prepare data for cross-validation:
dat_cv <- prep_data_cv(dat, K = 10)
```
