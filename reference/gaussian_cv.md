# Calculate Gaussian reference MSE, its decomposition, and a bound

Assumes jointly Gaussian responses with known mean and the supplied
stationary covariance, fixed locations, folds and target weights.
Prediction uses all available observations. Target observations have
fresh independent noise even when their coordinates coincide with
training locations. The calculation averages over repeated response
datasets; it is not an observed CV score.

## Usage

``` r
gaussian_cv(locations, folds, target, covariance)
```

## Arguments

- locations:

  Numeric two-column coordinates or projected spatial points.

- folds:

  A partition accepted by
  [`cv_folds()`](https://olatunjijohnson.github.io/spCVAudit/reference/cv_folds.md).

- target:

  A
  [`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md)
  or target coordinates.

- covariance:

  A
  [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md)
  object. A correlation function alone is insufficient because the
  calculation requires total variance.

## Value

A `spcv_gaussian` list with scalar `target_mse`, `cv_mse`, `optimism`,
`J2`, `one_neighbour`, `additional_neighbours`, `lower`, `upper`, and
`estimated_covariance`; and location-level `validation` and `target`
data frames of MSE reductions and bounds. Costs are cubic in
training-set size.

## Details

`optimism` is target MSE minus expected CV MSE. It equals
`one_neighbour` (total variance times `J2`) plus
`additional_neighbours`. The latter is the validation-minus-target
difference in the extra MSE reduction supplied by observations beyond
the nearest. With strictly positive nugget, `lower` and `upper` bound
optimism using conditional covariance residuals. These are deterministic
limits under the supplied model, not confidence intervals. Bounds are
`NA` when the nugget is zero. Supplying estimated covariance gives a
plug-in model calculation, not a guarantee about actual CV bias.

## Examples

``` r
xy <- cbind(c(-1, 0, 1), 0)
q <- prediction_target(matrix(c(2, 0), nrow = 1))
result <- gaussian_cv(xy, 1:3, q, spatial_covariance(1, 0.8))
result
#> Gaussian reference calculation
#>            target_mse                cv_mse              optimism 
#>           0.912868743           0.889812660           0.023056082 
#>         one_neighbour additional_neighbours                 lower 
#>           0.000000000           0.023056082          -0.002661682 
#>                 upper 
#>           0.116354041 
# Identical nearest distances, but full-neighbour optimism is about 0.02306.
```
