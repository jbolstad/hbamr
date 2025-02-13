# hbamr 2.4.2

### Revisions to existing models

-   The HBAM_R_MINI model now allows analyzing data with an even number of answering categories. 

# hbamr 2.4.1

### Revisions to existing models

-   All models except HBAM_R_MINI now allow analyzing data with an even number of answering categories. 

# hbamr 2.4.0

### Refactored code

-   All models have been collected in a single Stan file to reduce the footprint of the package. It now requires significantly less disk space. 

### Discontinued function

-   The `show_code()` function has been removed, as there is no longer easily readable code to show for any particular model. 

# hbamr 2.3.2

### Revisions to existing models

-  The following models now generate posterior predictions if the user runs `hbam()` with the argument `extra_pars = "Y_pred"`: HBAM, HBAM_MULTI, HBAM_NF, HBAM_MULTI_NF, BAM.

# hbamr 2.3.1

### Improvements of existing functions

-   The `prep_data()` function now prints a summary of the input and output data to the console. It also throws an error if the selection criteria for inclusion in the analysis are too strict to retain any respondents.
-   The appearance of the plot from `plot_stimuli()` has been slightly improved.

### Revisions to existing models

-   The priors on `theta` and `rho` have been made narrower to make the models more robust to situations with extremely scarce data. 

# hbamr 2.3.0

### Revisions to existing models

-   All FBAM-type models have been revised to allow for heteroskedastic errors. 
    -   The **FBAM_MINI** model has been renamed **FBAM** to reflect this change.  

# hbamr 2.2.1

### Revisions to existing models

-   The prior on `rho` (which appears in all heteroskedastic models) has been made narrower to reduce the risk of divergent transitions when there are few observations per stimuli (which could be the case for e.g. expert surveys). 
-   The default priors on `mu_beta` and `mu_alpha` in MULTI-type models and the priors on `sigma_beta` and `sigma_alpha` in HBAM-type models have also been made narrower to reduce the risk of sampling issues.

# hbamr 2.2.0

### Improvements of existing functions

-   The `get_est()` function now takes the logical argument `format_orig`, which if `TRUE` makes the function return posterior summaries for individual-level parameters in a format that matches the rows in the original dataset.
-   The `hbam()` and `fbam()` functions now store the input data within the returned objects, which allows simplifying the interfaces for other functions like `get_est()` and `plot_over_self()`. As a result, the `plot_over_self()` function no longer requires a data argument.
-   The `prep_data()`, `hbam()`, and `fbam()` functions now allow users to not supply self-placements. In this case, no meaningful respondent positions will be estimated, but all other parameters are unaffected.
-   The `prep_data()` function now allows the `group_id` argument to take various forms, such as factor or character. It also allows missing values in the `group_id` vector and will drop respondents who do not have a valid `group_id`. Missing values would previously generate an uninformative error message.
-   The models now accept even numbers of answering categories, but `prep_data()` will throw a warning.
-   The `prep_data()` function now identifies the left and right poles for the BAM model as the stimuli with the most non-NA observations on each side of the center. This can be advantageous when analyzing datasets where some stimuli have a much higher number of valid observations than others. 

# hbamr 2.1.2

### Revisions to existing models

-   The behavior of FBAM models now depends on whether they are fit via MCMC or optimization. When fit via MCMC, the models simulate errors in respondents' self-placements and generate draws for the latent discrete flipping parameters, kappa, using a Bernoulli function -- just like other models. When fit via optimization, there is no simulation of self-placement errors and lambda is rounded to yield a MAP estimate of kappa.

# hbamr 2.1.1

### Improved cross-validation function

-   The `hbam_cv()` function no longer uses `parallel::mclapply()` for parallel computation as the latter relies on forking, which is not available on Windows. `hbam_cv()` has been revised to work with the `future` package, where the user decides the computational strategy and options are available for parallel computation on all systems. 
-   The return value of `hbam_cv()` has been changed to comply with the standards of the `loo` package. The function now returns a list with classes `kfold` and `loo`. This allows the user to compare estimated ELPDs and obtain standard errors for their differences via `loo::loo_compare()`. 

# hbamr 2.1.0

### Revisions to existing models

-   All models in the package now simulate errors in respondents' self-placements. In the original HBAM model, respondents' latent positions were treated as parameters, which had several consequences: It yielded a realistic level of uncertainty in the estimated positions, but also led to slower sampling. Furthermore, the hierarchical prior on the latent respondent positions shrunk them toward zero, which was undesirable in most applications as it affected the distances between respondents and stimuli. Finally, the original approach led to a more nuanced posterior distribution where combinations of individual-level parameter values that would lead to implausible values for the respondent positions were weighted down. However, the model's log-normal prior on the stretch parameters is by itself sufficient to yield meaningful respondent positions estimates, so there is no significant cost to not treating the latent respondent positions as parameters. The current versions of the model instead simulate errors in the self-placements to yield the same level of uncertainty as the original model, while sampling considerably faster.

### Discontinued model

-   **HBAM_MAX** has been removed as it offered little extra after the other models were revised to simulate errors in respondents' self-placements.

# hbamr 2.0.1

### New function

-   `show_code()` shows Stan code for any model in the package.

### New models

-   **HBAM_MULTI** is a version that models differences between groups defined by the user. It requires a vector identifying the groups to be supplied as the argument `group_id`. The model gives each group separate hyperparameters for the locations of the prior distributions for the shift and stretch parameters. Rather than shrinking the estimates toward the mode for the whole dataset, this model shrinks the estimates toward the mode for the group. The vectors of hyperparameters are called `mu_alpha` and `mu_beta` and are constructed to have means of 0. The scales of the priors on these hyperparameters can be set by the user via the arguments `sigma_mu_alpha` and `sigma_mu_beta`. One potential use for this model is to supply self-placements as `group_id`, and thus give each self-placement group its own prior distribution for the shift and stretch parameters.
-   **HBAM_MULTI_NF** is a version of the HBAM_MULTI model that does not allow for scale flipping.
-   **FBAM_MULTI** is a version of the FBAM_MINI model that shares the group-modeling features of the HBAM_MULTI model. It allows the user to set the scales of the priors for the shift and stretch parameters via the arguments `sigma_alpha` and `sigma_beta`, and set the scales of the priors on `mu_alpha` and `mu_beta` via the arguments `sigma_mu_alpha` and `sigma_mu_beta`.
-   **FBAM_MULTI_NF** is a version of the FBAM_MULTI model that does not allow for scale flipping.

### Renamed models

-   The original HBAM model has been renamed **HBAM_MAX**.
-   The HBAM_NE model has been renamed **HBAM** as it makes a better default model.
-   HBAM_0 has been renamed **HBAM_NF**. To be consistent with the other models and offer faster sampling, the model no longer models the latent respondent positions as parameters.

### Revisions to existing models

-   **FBAM_MINI** (as well as the new FBAM models) now allows for user-defined prior distributions.
-   All models (except HBAM_MAX) have been revised to not treat the respondent positions as parameters. This results in considerably faster sampling.
-   All HBAM models that allow for scale flipping have been given a logit-normal prior on the mixing proportions, lambda (i.e. the expectations of the latent discrete flipping parameters, kappa). This replaces the original beta prior, which could trigger divergent transitions.

### Discontinued models

-   **HBAM_2** has been replaced by the more general HBAM_MULTI model.
-   **HBAM_HM** has been removed as it offered little extra and the number of models should be kept somewhat limited.
-   **HBAM_R** has been removed as users would typically be better off with the simpler HBAM_R_MINI.
