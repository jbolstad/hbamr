# hbamr: Hierarchical Bayesian Aldrich-McKelvey Scaling in R via Stan

An R package for performing hierarchical Bayesian Aldrich-McKelvey (HBAM) scaling using Hamiltonian Monte Carlo simulations via Stan. Aldrich-McKelvey (AM) scaling is a method for estimating the ideological positions of survey respondents and political actors on a common scale using positional survey data (Aldrich & McKelvey 1977). The hierarchical versions of the AM model included in this package outperform other versions by a considerable margin both in terms of yielding meaningful posterior distributions for respondent positions and in terms of recovering true respondent positions in simulations (Bølstad 2023). The package contains functions for preparing data, fitting models, extracting estimates, plotting key results, and comparing models using cross-validation.

### Usage

```{r}
library("hbamr")
```

Load and re-code example data:

```{r}
data(LC1980)
LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA 
self <- LC1980[, 1]
stimuli <- LC1980[, -1]
```

Fit standard HBAM model:

```{r}
fit_hbam <- hbam(self, stimuli)
```

Fit HBAM_MINI model:

```{r}
fit_hbam_mini <- hbam(self, stimuli, model = "HBAM_MINI")
```

Plot estimated stimuli positions:
```{r}
plot_stimuli(fit_hbam)
```

<img src="https://github.com/jbolstad/hbamr/blob/f2ddbcac26b56e59b8fad22a898a6fee145f06c5/vignettes/p_stim.svg?raw=true" width="700px">

Plot distribution of estimated respondent positions:

```{r}
plot_respondents(fit_hbam)
```

<img src="https://github.com/jbolstad/hbamr/blob/f2ddbcac26b56e59b8fad22a898a6fee145f06c5/vignettes/p_resp.svg?raw=true" width="700px">

Plot estimated scale stretching by self-placement:

```{r}
dat <- prep_data(self, stimuli)
plot_over_self(list(fit_hbam, fit_hbam_mini), dat, "abs_beta")
```

<img src="https://github.com/jbolstad/hbamr/blob/f2ddbcac26b56e59b8fad22a898a6fee145f06c5/vignettes/p_abs_beta.svg?raw=true" width="700px">

### References

-   Aldrich, John H, and Richard D McKelvey. 1977. "A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections." *American Political Science Review* 71(1): 111--130.
-   Bølstad, Jørgen. 2020. "Capturing Rationalization Bias and Differential Item Functioning: A Unified Bayesian Scaling Approach." *Political Analysis* 28(3): 340--355.
-   Bølstad, Jørgen. 2023. "Hierarchical Bayesian Aldrich-McKelvey Scaling." *Political Analysis*.
-   Hare, Christopher et al. 2015. "Using Bayesian Aldrich-McKelvey scaling to study citizens' ideological preferences and perceptions." *American Journal of Political Science* 59(3): 759--774.
