# Compare validation and target nearest-neighbour distances

Calculates exact integrals for weighted empirical distributions, without
a numerical integration grid. `A` is mean target distance minus mean
validation distance; `W` is their one-dimensional Wasserstein distance.
`J1` and `J2` are validation minus target averages of correlation and
squared correlation. `U1` and `U2` are Wasserstein distances between
those transformed values.

## Usage

``` r
distance_contrasts(
  validation,
  target,
  correlation = NULL,
  validation_weights = NULL,
  target_weights = NULL
)
```

## Arguments

- validation, target:

  Nonempty numeric vectors of finite, nonnegative nearest-training
  distances.

- correlation:

  A function mapping distances to correlations in `[-1, 1]`, or a
  [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md)
  object. If `NULL`, correlation diagnostics are `NA`.

- validation_weights, target_weights:

  Nonnegative weights, normalised separately. Defaults are equal weights
  within each distribution.

## Value

A named numeric vector: `A`, `W`, `J1`, `J2`, `U1`, `U2`, `Aplus`,
`Aminus`, and `crosses` (one if both positive and negative CDF areas
exceed a relative floating-point tolerance, otherwise zero).

## Details

A positive `J2` means higher average squared nearest-observation
correlation during validation. Under Gaussian prediction using one
neighbour with known mean and covariance, target MSE minus expected CV
MSE equals total variance times `J2`. This is not a general bias
correction for fitted prediction models. The integral interpretation of
`U1` and `U2` as unsigned weighted distance-CDF contrasts assumes
nonnegative correlation decreasing with distance.

## Examples

``` r
model <- spatial_covariance(range = 1, signal = 0.8)
distance_contrasts(c(0.2, 0.5, 1), c(0.5, 1, 2), model)
#>         A         W        J1        J2        U1        U2     Aplus    Aminus 
#> 0.6000000 0.6000000 0.1822388 0.1390943 0.1822388 0.1390943 0.6000000 0.0000000 
#>   crosses 
#> 0.0000000 
```
