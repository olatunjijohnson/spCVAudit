# Audit a CV partition against a prediction target

Each validation observation is compared with its fold's available
training observations. Each target location is compared with the full
training set. Validation observations receive weight `1/n`; target
weights come from
[`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md).
This function needs no response values.

## Usage

``` r
audit_cv(locations, folds, target, covariance = NULL)
```

## Arguments

- locations:

  Numeric two-column coordinates or projected spatial points.

- folds:

  Fold labels or an explicit list accepted by
  [`cv_folds()`](https://olatunjijohnson.github.io/spCVAudit/reference/cv_folds.md).

- target:

  A
  [`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md)
  object, or target coordinates.

- covariance:

  A
  [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md)
  or
  [`estimate_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/estimate_covariance.md)
  object, a correlation function, or `NULL` for distance diagnostics
  only.

## Value

An object of class `spcv_audit` containing `diagnostics`, validation and
target nearest distances, normalised target weights, validated folds,
and the covariance specification. `diagnostics` is described in
[`distance_contrasts()`](https://olatunjijohnson.github.io/spCVAudit/reference/distance_contrasts.md).

## Examples

``` r
xy <- cbind(seq_len(12) / 12, rep(c(0, 1), 6))
q <- prediction_target(rbind(c(0.2, 0.5), c(0.8, 0.5)))
folds <- spatial_folds(xy, k = 3)
audit_cv(xy, folds, q, spatial_covariance(0.3, 0.8))
#> Spatial CV audit
#>         A         W        J1        J2        U1        U2 
#> 0.2649988 0.2649988 0.2338896 0.1370341 0.2338896 0.1370341 
#> Diagnostics compare prediction conditions; they do not certify CV accuracy.
```
