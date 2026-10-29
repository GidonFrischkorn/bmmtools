# Build a correlation matrix from factor loadings

The model-implied correlation matrix of a factor model, for the `cors`
argument of
[`simulate_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/simulate_recovery.md)
and
[`recovery_grid()`](https://www.gfrischkorn.org/bmmtools/reference/recovery_grid.md):
tasks that load on a common ability, parameters that share a latent
source, or a covariate that measures the factor. The result is `L Φ Lᵀ`
with its diagonal set to 1, so each term's unique variance is one minus
its communality.

## Usage

``` r
cors_from_factors(loadings, factor_cors = NULL)
```

## Arguments

- loadings:

  A numeric matrix, terms by factors, with the terms as row names; or a
  named numeric vector for a single factor.

- factor_cors:

  The factor correlation matrix. `NULL` means uncorrelated factors. When
  both it and `loadings` carry factor names, they must agree.

## Value

A correlation matrix with the terms as dimnames.

## Details

This is the **population** matrix. A simulation that instead takes the
empirical correlation of a finite sample of factor scores (as the indSim
scripts do with 1000 lavaan draws) adds sampling noise to the generating
values.

## Examples

``` r
# two tasks per ability, the abilities correlated .5
loadings <- matrix(
  c(0.8, 0.7, 0, 0, 0, 0, 0.75, 0.85),
  ncol = 2,
  dimnames = list(
    c("kappa_task1", "kappa_task2", "thetat_task1", "thetat_task2"),
    c("precision", "memory")
  )
)
phi <- matrix(c(1, 0.5, 0.5, 1), 2,
  dimnames = rep(list(c("precision", "memory")), 2)
)
cors_from_factors(loadings, phi)
#>              kappa_task1 kappa_task2 thetat_task1 thetat_task2
#> kappa_task1         1.00      0.5600       0.3000       0.3400
#> kappa_task2         0.56      1.0000       0.2625       0.2975
#> thetat_task1        0.30      0.2625       1.0000       0.6375
#> thetat_task2        0.34      0.2975       0.6375       1.0000

# one factor
cors_from_factors(c(a = 0.8, b = 0.6, G = 0.5))
#>      a    b   G
#> a 1.00 0.48 0.4
#> b 0.48 1.00 0.3
#> G 0.40 0.30 1.0
```
