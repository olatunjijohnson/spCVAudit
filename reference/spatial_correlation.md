# Evaluate correlation between distinct spatial observations

Evaluate correlation between distinct spatial observations

## Usage

``` r
spatial_correlation(model, distance)
```

## Arguments

- model:

  A
  [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md)
  or
  [`estimate_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/estimate_covariance.md)
  object.

- distance:

  Finite, nonnegative distances.

## Value

A numeric vector with one correlation per distance. At zero distance
this is the signal fraction, because the observations have independent
noise.

## Examples

``` r
spatial_correlation(spatial_covariance(1, 0.8), c(0, 1, 2))
#> [1] 0.8000000 0.2943036 0.1082682
```
