# Evaluate conditional target MSE under a specified Gaussian distribution

This reference calculation uses the supplied mean and covariance as the
generating distribution. Conditional on observed training responses,
target MSE is conditional variance plus squared difference between the
supplied predictions and the true conditional mean. It is available in a
simulation with known generating parameters; estimated parameters do not
make it an exact evaluation of real-world prediction error.

## Usage

``` r
gaussian_target_mse(
  locations,
  response,
  target,
  covariance,
  predictions = NULL,
  mean = 0
)
```

## Arguments

- locations:

  Numeric two-column training coordinates or projected points.

- response:

  Finite training responses.

- target:

  A
  [`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md)
  or target coordinates.

- covariance:

  The specified generating
  [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md).

- predictions:

  Finite predictions at each target location. If `NULL`, uses the
  conditional mean under the supplied Gaussian distribution.

- mean:

  Known constant mean of the generating distribution.

## Value

A list with weighted `mse`, vectors `pointwise_mse`, `conditional_mean`,
`conditional_variance`, and `predictions`.

## Examples

``` r
xy <- cbind(c(0, 1, 2), 0)
gaussian_target_mse(xy, c(0.2, -0.1, 0.5), matrix(c(1.5, 0), 1),
                    spatial_covariance(1, 0.8))
#> $mse
#> [1] 0.6353526
#> 
#> $pointwise_mse
#> [1] 0.6353526
#> 
#> $conditional_mean
#> [1] 0.1565101
#> 
#> $conditional_variance
#> [1] 0.6353526
#> 
#> $predictions
#> [1] 0.1565101
#> 
```
