# hbamr 2.0.1

## New features

-   New function: `show_code()` shows Stan code for any model in the package.
-   New models:
    -   **HBAM_MULTI** is a version of the HBAM_MINI model that models differences between groups defined by the user. It requires an integer vector identifying the groups to be supplied as the argument `group_id`. The model gives each group separate hyperparameters for the locations of the prior distributions for the shift and stretch parameters. Rather than shrinking the estimates toward the mode for the whole dataset, this model shrinks the estimates toward the mode for the group. The vectors of hyperparameters are called `mu_alpha` and `mu_beta` and are constructed to have means of 0. One potential use for this model is to supply self-placements (appropriately transformed) as `group_id`, and thus give each self-placement group its own prior distribution for the shift and stretch parameters.
    -   **HBAM_MULTI_0** is a version of the HBAM_NE model that does not allow for scale flipping, but shares the group-modeling features of the HBAM_MULTI model. If differs from HBAM_MULTI by modeling heteroskedastic prediction errors and not allowing for scale flipping.
    -   **FBAM_MULTI** is a version of the FBAM_MINI model that shares the group-modeling features of the HBAM_MULTI model.
    -   **FBAM_MULTI_0** is a version of the FBAM_MULTI model that does not allow for scale flipping.
