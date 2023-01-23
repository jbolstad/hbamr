# hbamr: Hierarchical Bayesian Aldrich-McKelvey Scaling in R via Stan

An R package for performing hierarchical Bayesian Aldrich-McKelvey (HBAM) scaling using Hamiltonian Monte Carlo simulations via Stan. Aldrich-McKelvey (AM) scaling is a method for estimating the ideological positions of survey respondents and political actors on a common scale using positional survey data (Aldrich & McKelvey 1977). The hierarchical versions of the AM model included in this package outperform other versions by a considerable margin both in terms of yielding meaningful posterior distributions for respondent positions and in terms of recovering true respondent positions in simulations (Bølstad 2023). The package contains functions for preparing data, fitting models, extracting estimates, plotting key results, and comparing models using cross-validation.

### Installing and loading the package

The **hbamr** package requires a working installation of **rstan**. Instructions for installing **rstan** are available [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). (As explained in these instructions, you need to configure your R installation to be able to compile C++ code and you can drastically improve the estimation speed of the models by optimizing compiler usage.)

Users of the **hbamr** package should install the latest version from github -- it is not available on CRAN. Those who do not have the **devtools** package installed, should first install this package by running `install.packages("devtools")`. **hbamr** can then be installed by running `devtools::install_github("jbolstad/hbamr")`. Once this has been done, the package can be loaded by running `library("hbamr")`.

### References

Aldrich, John H, and Richard D McKelvey. 1977. "A Method of Scaling with Applications to the 1968 and 1972 Presidential Elections." *American Political Science Review* 71(1): 111--130.

Bølstad, Jørgen. 2020. "Capturing Rationalization Bias and Differential Item Functioning: A Unified Bayesian Scaling Approach." *Political Analysis* 28(3): 340--355.

Bølstad, Jørgen. 2023. "Hierarchical Bayesian Aldrich-McKelvey Scaling." *Political Analysis*.

Hare, Christopher et al. 2015. "Using Bayesian Aldrich-McKelvey scaling to study citizens' ideological preferences and perceptions." *American Journal of Political Science* 59(3): 759--774.
