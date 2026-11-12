# Compare candidate CV schemes using a common diagnostic

`A`, `J1` and `J2` are ranked by absolute value; `W`, `U1` and `U2` by
value. All candidates use the same supplied covariance and target.
Schemes within `tolerance` of the minimum are tied and receive equal
selection probabilities. No observed target errors are used to select or
break ties. When evaluating such a rule, average the individual tied
schemes' losses with these probabilities; do not average their signed
errors before taking absolute values.

## Usage

``` r
compare_cv_schemes(
  locations,
  schemes,
  target,
  covariance = NULL,
  criterion = c("J2", "W", "J1", "A", "U1", "U2"),
  tolerance = 1e-12
)
```

## Arguments

- locations:

  Numeric two-column coordinates or projected spatial points.

- schemes:

  A nonempty named list of candidate fold partitions.

- target:

  A
  [`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md)
  or target coordinates.

- covariance:

  A covariance object or correlation function; required for
  correlation-based selection.

- criterion:

  One of `"J2"`, `"W"`, `"J1"`, `"A"`, `"U1"`, or `"U2"`.

- tolerance:

  Absolute numerical tolerance for ties, in criterion units.

## Value

A `spcv_comparison` list with `table` (one row per scheme, diagnostics,
`criterion_value`, `selected`, and `selection_probability`), individual
`audits`, the selection criterion and tolerance.

## Examples

``` r
xy <- cbind(seq_len(20) / 20, rep(c(0, 1), 10))
q <- prediction_target(rbind(c(0.2, 0.5), c(0.8, 0.5)))
schemes <- list(random = spatial_folds(xy, 4),
                block = spatial_folds(xy, 4, "block"))
compare_cv_schemes(xy, schemes, q, spatial_covariance(0.3, 0.8))
#> CV comparison by |J2| 
#>  scheme     W        J2 criterion_value selection_probability
#>  random 0.395 0.2977614       0.2977614                     0
#>   block 0.200 0.1074095       0.1074095                     1
```
