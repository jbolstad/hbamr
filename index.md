# hbamr

### Hierarchical Bayesian Aldrich-McKelvey Scaling

This R package implements hierarchical Bayesian Aldrich-McKelvey (HBAM)
scaling using Hamiltonian Monte Carlo via Stan. Aldrich-McKelvey (AM)
scaling is a method for estimating the latent positions of survey
respondents and external objects on a common scale using positional
survey data (Aldrich & McKelvey 1977). The hierarchical versions of the
AM model included in this package outperform other versions both in
terms of yielding meaningful posterior distributions for respondent
positions and in terms of recovering true respondent positions in
simulations ([Bølstad 2024](https://doi.org/10.1017/pan.2023.18)). The
package provides functions for preparing data, fitting models,
extracting estimates, plotting key results, and comparing models using
cross-validation.

### Installation

The package is available from
[CRAN](https://CRAN.R-project.org/package=hbamr) and can be installed
using the standard method:

``` r

install.packages("hbamr")
```

This is the easiest way to install the package, as the binaries on CRAN
include pre-compiled models that are ready for use.

### Usage

Load the package:

``` r

library("hbamr")
```

Load and recode example data:

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

Plot estimated stimuli positions:

``` r

plot_stimuli(fit_hbam)
```

![](https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_stim.svg)

Plot the distribution of estimated respondent positions:

``` r

plot_respondents(fit_hbam)
```

![](https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_resp.svg)

Plot estimated scale-stretching over respondents’ self-placements:

``` r

plot_over_self(list(fit_hbam, fit_hbam_mini), "abs_beta")
```

![](https://raw.githubusercontent.com/jbolstad/hbamr/main/vignettes/p_abs_beta.svg)

### References

- Aldrich, John H, and Richard D McKelvey. 1977. “A Method of Scaling
  with Applications to the 1968 and 1972 Presidential Elections.”
  *American Political Science Review* 71(1): 111-130.
- Bølstad, Jørgen. 2020. “Capturing Rationalization Bias and
  Differential Item Functioning: A Unified Bayesian Scaling Approach.”
  *Political Analysis* 28(3): 340-355.
- Bølstad, Jørgen. 2024. “[Hierarchical Bayesian Aldrich-McKelvey
  Scaling](https://doi.org/10.1017/pan.2023.18).” *Political Analysis*
  32(1): 50-64.
- Hare, Christopher et al. 2015. “Using Bayesian Aldrich-McKelvey
  Scaling to Study Citizens’ Ideological Preferences and Perceptions.”
  *American Journal of Political Science* 59(3): 759-774.
