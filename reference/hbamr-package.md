# Hierarchical Bayesian Aldrich-McKelvey Scaling via Stan

Fit hierarchical Bayesian Aldrich-McKelvey (HBAM) models using a form of
Hamiltonian Monte Carlo via Stan. Aldrich-McKelvey (AM) scaling is a
method for estimating the ideological positions of survey respondents
and political actors on a common scale using positional survey data. The
hierarchical versions of the Bayesian AM model included in this package
outperform other versions both in terms of yielding meaningful posterior
distributions for respondent positions and in terms of recovering true
respondent positions in simulations. The package contains functions for
preparing data, fitting models, extracting estimates, plotting key
results, and comparing models using cross-validation.

## References

- Bølstad, Jørgen (2024). Hierarchical Bayesian Aldrich-McKelvey
  Scaling. *Political Analysis*. 32(1): 50–64.
  [doi:10.1017/pan.2023.18](https://doi.org/10.1017/pan.2023.18) .

- Stan Development Team (2024). RStan: the R interface to Stan.
  <https://mc-stan.org>.

## See also

- <https://jbolstad.github.io/hbamr/>

## Author

Jørgen Bølstad
