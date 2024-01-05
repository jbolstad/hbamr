# hbamr 2.0.1

## New features

-   New function: `show_code()` shows Stan code for any model in the package.
-   New models:
    -   **HBAM_MULTI** is a version that models differences between groups defined by the user. It requires an integer vector identifying the groups to be supplied as the argument `group_id`. The model gives each group separate hyperparameters for the locations of the prior distributions for the shift and stretch parameters. Rather than shrinking the estimates toward the mode for the whole dataset, this model shrinks the estimates toward the mode for the group. The vectors of hyperparameters are called `mu_alpha` and `mu_beta` and are constructed to have means of 0. The scales of the priors on these hyperparameters can be set by the user via the arguments `sigma_mu_alpha` and `sigma_mu_beta`. One potential use for this model is to supply self-placements (appropriately transformed) as `group_id`, and thus give each self-placement group its own prior distribution for the shift and stretch parameters.
    -   **HBAM_MULTI_NF** is a version of the HBAM_MULTI model that does not allow for scale flipping.
    -   **FBAM_MULTI** is a version of the FBAM_MINI model that shares the group-modeling features of the HBAM_MULTI model. It allows the user to set the scales of the priors for the shift and stretch parameters via the arguments `sigma_alpha` and `sigma_beta`, and set the scales of the priors on `mu_alpha` and `mu_beta` via the arguments `sigma_mu_alpha` and `sigma_mu_beta`.
    -   **FBAM_MULTI_NF** is a version of the FBAM_MULTI model that does not allow for scale flipping.
-   Renamed models:
    -   The original HBAM model has been renamed **HBAM_MAX**.
    -   The HBAM_NE model has been renamed **HBAM** as it makes a better default model. 
    -   HBAM_0 has been renamed **HBAM_NF**. To be consistent with the other models and offer faster sampling, the model no longer models the latent respondent positions as parameters.
- Revisions to existing models:
    -   **FBAM_MINI** (as well as the new FBAM models) now allows for user-defined prior distributions.
    -   All models (except HBAM_MAX) have been revised to not treat the respondent positions as parameters. This results in considerably faster sampling. 
    -   All HBAM models that allow for scale flipping have been given a logit-normal prior on the mixing proportions, lambda (i.e. the expectations of the latent discrete flipping parameters, kappa). This replaces the original beta prior, which could trigger divergent transitions. 
-   Deprecated models:
    -   **HBAM_2** has been replaced by the more general HBAM_MULTI model.
    -   **HBAM_HM** has been removed as it offered little extra and the number of models should be kept somewhat limited.
    -   **HBAM_R** has been removed as users would typically be better off with the simpler HBAM_R_MINI.


