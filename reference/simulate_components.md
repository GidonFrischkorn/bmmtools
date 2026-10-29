# Simulate several models for the same simulated people

Draws the subject values of every component, and any covariates, from
one multivariate normal distribution, then simulates each component's
data with its own model and generator. Parameters of different models
can therefore be correlated across subjects, which is what a study of
individual differences across tasks or models needs.

## Usage

``` r
simulate_components(
  components,
  n_subjects,
  cors = NULL,
  covariates = NULL,
  seed = NULL
)
```

## Arguments

- components:

  A list of
  [`recovery_component()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_component.md)s.
  If the list has names, they must be the component names.

- n_subjects:

  The number of subjects, shared by every component.

- cors:

  A correlation matrix over prefixed terms (`<name>_<term>`, such as
  `m3_c_task1`) and covariates; terms it does not name are uncorrelated.
  `NULL` means none are correlated. May be a function with no arguments,
  as in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).

- covariates:

  Observed, error-free person variables drawn with the parameters, as in
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md).
  A covariate name must differ from every component name and from every
  model's parameters and columns.

- seed:

  A seed applied with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  around everything that is drawn; `NULL` leaves the random number
  generator alone.

## Value

A list of class `bmmtools_simulation_set` with

- `components`: a named list of `bmmtools_simulation`s, one per
  component, with its own unprefixed terms, no covariates and `seed`
  `NA`;

- `truth`: the tables of
  [`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
  (`population`, `subjects`, `sd`, `cor`, `covariates`) over every
  component, with prefixed terms; `cor` holds every pair of varying
  terms and covariates;

- `pars`, `sds`: the realised values, prefixed and expanded over tasks;

- `cors`: the realised correlation matrix over the varying terms, then
  the covariates (`NULL` with fewer than two);

- `covariates` and `covariate_data`, a tibble with `id` and one column
  per covariate;

- `links`: every model's link table under prefixed names;

- `specs`: the components as given;

- `n_subjects` and `seed` (`NA` when none).

## Details

**Order of the draws.** Inside the seed, each component's `pars` and
then `sds` are evaluated, in component order, then `cors`. The subject
values of all components are one draw, `Z %*% chol(cors)`, over the
varying terms, component by component. The generators then run in
component order, and the covariates are drawn last, so adding a
covariate never changes the responses. A set with one component and no
covariates gives exactly the simulation
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
gives with the same arguments and seed. Adding a component after another
leaves the first one's subject values unchanged, but not its responses,
because the draw for the added terms comes before the generators.

**Covariates** are not added to the components' data; they are in
`covariate_data` and in `truth$covariates`.

**What separate fits can recover.** Each component is fitted on its own
by
[`fit_components()`](https://www.gfrischkorn.org/bmmtools/reference/fit_components.md),
so no fit knows about the correlation between its parameters and another
model's. Correlating the estimates across fits attenuates the
correlation by the reliabilities `rel_a` and `rel_b` of the two subject
estimates: the `point` estimator of
[`extract_correlations()`](https://www.gfrischkorn.org/bmmtools/reference/extract_correlations.md)
is roughly `rho * sqrt(rel_a * rel_b)`, and the `draws` estimator
roughly `rho * rel_a * rel_b`, so neither recovers `rho`. A structural
equation model of true against estimated values, from
[`subject_table()`](https://www.gfrischkorn.org/bmmtools/reference/subject_table.md),
does; so would a joint multivariate fit of all components, which
bmmtools does not support yet.

## Examples

``` r
if (FALSE) { # \dontrun{
components <- list(
  recovery_component(
    bmm::mixture2p(resp_error = "y"),
    pars = c(kappa = log(8), thetat = qlogis(0.75)),
    n_trials = 60, sds = c(kappa = 0.3, thetat = 0.5), name = "m2p"
  ),
  recovery_component(
    bmm::sdm(resp_error = "y"),
    pars = c(c = log(3), kappa = log(5)),
    n_trials = 60, sds = c(c = 0.4), name = "sdm"
  )
)
cors <- diag(2)
dimnames(cors) <- rep(list(c("m2p_thetat", "sdm_c")), 2)
cors[1, 2] <- cors[2, 1] <- 0.6
set <- simulate_components(components,
  n_subjects = 100, cors = cors,
  seed = 1
)
fits <- fit_components(set, dir = "fits")
recover_correlations(fits, set, estimator = c("draws", "point"))
} # }
```
