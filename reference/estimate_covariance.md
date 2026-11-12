# Estimate exponential covariance by profiling a binned variogram

Fits `nugget + psill * (1 - exp(-distance / range))` to an empirical
semivariogram, with bin weights proportional to pair counts. For each
range, constrained nugget and partial sill are obtained by evaluating
the quadratic optimum and all parameter boundaries. An 81-point
log-range grid and local refinements reduce dependence on optimiser
starting values. This does not account for covariance estimation
uncertainty or establish model adequacy.

## Usage

``` r
estimate_covariance(
  locations,
  response,
  n_bins = 15L,
  cutoff = NULL,
  max_pairs = 100000L,
  seed = 1L
)
```

## Arguments

- locations:

  Numeric two-column coordinates or projected spatial points.

- response:

  Finite numeric observations or residuals, one per location.

- n_bins:

  Number of equal-width distance bins.

- cutoff:

  Maximum pair distance, default half the maximum sampled distance.

- max_pairs:

  Maximum number of unordered pairs; larger datasets use a sample
  without replacement, before allocating pair coordinates.

- seed:

  Seed for pair sampling; restores the caller's random-number state.

## Value

A `spcv_covariance` object with `estimated = TRUE`, its empirical
`variogram`, and `fit` diagnostics including objective, boundary flags
and scaled parameter limits. Boundary fits remain visible for
assessment.

## Details

The fit is isotropic. Supply residuals as `response` to describe
residual rather than marginal dependence. To avoid leakage in CV
predictions, estimate covariance separately using each fold's training
responses. A common estimate from all training responses may be used to
compare diagnostics across schemes; that selection then depends on
responses as well as locations.

## Examples

``` r
xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
                           y = seq(0, 1, length.out = 6)))
y <- sin(3 * xy[, 1]) + cos(5 * xy[, 2])
fit <- estimate_covariance(xy, y)
fit
#> Estimated exponential covariance
#>   scale: 1.237599  signal: 1  total variance: 1.955322 
#>   fit boundary: zero nugget, spatial variance upper limit 
```
