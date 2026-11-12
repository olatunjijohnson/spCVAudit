# Construct random, spatial-block, buffered, or kNNDM folds

Spatial blocks use coordinate k-means. Buffered folds remove training
observations within `buffer` of any test observation. The kNNDM method
calls
[`CAST::knndm()`](https://hannameyer.github.io/CAST/reference/knndm.html)
and requires the optional packages 'CAST' and 'sf'. Its target must have
equal weights. All methods return explicit indices suitable for
[`audit_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/audit_cv.md).
These constructors do not use response values.

## Usage

``` r
spatial_folds(
  locations,
  k = 5L,
  method = c("random", "block", "buffered", "knndm"),
  buffer = NULL,
  target = NULL,
  seed = 1L
)
```

## Arguments

- locations:

  Numeric two-column coordinates or projected spatial points.

- k:

  Number of folds, between two and the number of observations.

- method:

  One of `"random"`, `"block"`, `"buffered"`, or `"knndm"`.

- buffer:

  Nonnegative exclusion distance, required for buffered CV and expressed
  in the same units as the coordinates.

- target:

  Prediction target, required for kNNDM.

- seed:

  An integer seed. The default is reproducible and restores the caller's
  random-number state, including after calls to external constructors.
  Use `NULL` to use and advance the caller's random-number stream.

## Value

Validated
[`cv_folds()`](https://olatunjijohnson.github.io/spCVAudit/reference/cv_folds.md)
with explicit training and test indices.

## Examples

``` r
xy <- cbind(seq_len(12), rep(c(0, 1), 6))
spatial_folds(xy, k = 3, method = "random", seed = 10)
#> [[1]]
#> [[1]]$train
#> [1]  2  3  4  6  7  9 11 12
#> 
#> [[1]]$test
#> [1]  1  5  8 10
#> 
#> 
#> [[2]]
#> [[2]]$train
#> [1]  1  3  4  5  8 10 11 12
#> 
#> [[2]]$test
#> [1] 2 6 7 9
#> 
#> 
#> [[3]]
#> [[3]]$train
#> [1]  1  2  5  6  7  8  9 10
#> 
#> [[3]]$test
#> [1]  3  4 11 12
#> 
#> 
#> attr(,"class")
#> [1] "spcv_folds" "list"      
#> attr(,"n")
#> [1] 12
```
