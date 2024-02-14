# hbamr <a href="https://jbolstad.github.io/hbamr/"><img src="man/figures/logo.png" align="right" height="139" alt="hbamr website" /></a>

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/hbamr)](https://CRAN.R-project.org/package=hbamr) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/hbamr)](https://cran.r-project.org/package=hbamr) [![R-CMD-check](https://github.com/jbolstad/hbamr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jbolstad/hbamr/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

### Hierarchical Bayesian Aldrich-McKelvey Scaling

This is an R package for performing hierarchical Bayesian Aldrich-McKelvey (HBAM) scaling using Hamiltonian Monte Carlo simulations via Stan. Aldrich-McKelvey (AM) scaling is a method for estimating the ideological positions of survey respondents and political actors on a common scale using positional survey data (Aldrich & McKelvey 1977). The hierarchical versions of the AM model included in this package outperform other versions both in terms of yielding meaningful posterior distributions for respondent positions and in terms of recovering true respondent positions in simulations ([Bølstad 2024](https://doi.org/10.1017/pan.2023.18)). The package contains functions for preparing data, fitting models, extracting estimates, plotting key results, and comparing models using cross-validation.

### Important Updates

**Version 2.1.0**:

-   All models now simulate errors in respondents' self-placements to yield realistic levels of uncertainty in estimated respondent positions while offering faster sampling.

**Version 2.0.1**:

-   New MULTI-type models explicitly model group-differences.
-   Models of FBAM-type now allow users to specify key priors.
-   Most models have been revised to offer faster and better sampling.

See the [changelog](https://jbolstad.github.io/hbamr/news/index.html) for a more comprehensive discussion of the updates.

### Installation

The package is available from [CRAN](https://CRAN.R-project.org/package=hbamr) and can be installed using the standard method:

``` r
install.packages("hbamr")
```

This is the easiest and fastest way to install the package, as the binaries on CRAN include pre-compiled models that are ready for use.

### Vignette

A vignette showing how to use all key functions in the package is available [here](https://jbolstad.github.io/hbamr/articles/hbamr.html). It can also be viewed locally, after installing the package:

``` r
vignette("hbamr")
```

### Usage

Load the package:

``` r
library("hbamr")
```

Load and re-code example data:

``` r
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA 
self <- LC1980[, 1]
stimuli <- LC1980[, -1]
```

Fit the standard HBAM model:

``` r
fit_hbam <- hbam(self, stimuli)
```

Fit the HBAM_MINI model:

``` r
fit_hbam_mini <- hbam(self, stimuli, model = "HBAM_MINI")
```

Plot the estimated stimuli positions:

``` r
plot_stimuli(fit_hbam)
```

<img src="https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_stim.svg" width="850px"/>

Plot the distribution of estimated respondent positions:

``` r
plot_respondents(fit_hbam)
```

<img src="https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_resp.svg" width="850px"/>

Plot the estimated scale-stretching parameters over respondents' self-placements:

``` r
plot_over_self(list(fit_hbam, fit_hbam_mini), "abs_beta")
```

<img src="https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_abs_beta.svg" width="850px"/>

### References

-   Aldrich, John H, and Richard D McKelvey. 1977. "A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections." *American Political Science Review* 71(1): 111-130.
-   Bølstad, Jørgen. 2020. "Capturing Rationalization Bias and Differential Item Functioning: A Unified Bayesian Scaling Approach." *Political Analysis* 28(3): 340-355.
-   Bølstad, Jørgen. 2024. "[Hierarchical Bayesian Aldrich-McKelvey Scaling](https://doi.org/10.1017/pan.2023.18)." *Political Analysis* 32(1): 50-64.
-   Hare, Christopher et al. 2015. "Using Bayesian Aldrich-McKelvey Scaling to Study Citizens' Ideological Preferences and Perceptions." *American Journal of Political Science* 59(3): 759-774.
