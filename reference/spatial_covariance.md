# Specify a stationary isotropic spatial covariance

`range` is the scale parameter, not a common practical correlation range
across families. For distinct observations, correlation is
`signal * exp(-d/range)` for the exponential family and
`signal * (1 + sqrt(3)*d/range) * exp(-sqrt(3)*d/range)` for the
Matern-3/2 family. Independent observation noise has variance
`variance * (1 - signal)`, including at coincident locations.

## Usage

``` r
spatial_covariance(
  range,
  signal = 1,
  variance = 1,
  family = c("exponential", "matern32")
)
```

## Arguments

- range:

  Positive correlation scale in coordinate-distance units.

- signal:

  Spatial signal fraction in `[0, 1]`.

- variance:

  Positive total response variance.

- family:

  Correlation family, `"exponential"` or `"matern32"`.

## Value

A `spcv_covariance` object with the supplied parameters, spatial
variance `psill`, independent noise variance `nugget`, and
`estimated = FALSE`.

## Examples

``` r
model <- spatial_covariance(range = 0.2, signal = 0.8)
spatial_correlation(model, c(0, 0.2, 1))
#> [1] 0.800000000 0.294303553 0.005390358
```
