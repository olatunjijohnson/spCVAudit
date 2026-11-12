# Plot validation and target distance distributions

The solid curve is validation and the dashed curve is the weighted
target.

## Usage

``` r
# S3 method for class 'spcv_audit'
plot(x, ...)
```

## Arguments

- x:

  A
  [`audit_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/audit_cv.md)
  result.

- ...:

  Additional arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, the audit object. Called for its plotting side effect.

## Examples

``` r
xy <- cbind(1:8, rep(c(0, 1), 4))
a <- audit_cv(xy, rep(1:2, 4), rbind(c(2, 2), c(6, 2)))
plot(a)
```
