# Package index

## Define your prediction task

Declare the target and construct or import a complete CV partition.

- [`prediction_target()`](https://olatunjijohnson.github.io/spCVAudit/reference/prediction_target.md)
  : Declare the locations and weights of a prediction target
- [`cv_folds()`](https://olatunjijohnson.github.io/spCVAudit/reference/cv_folds.md)
  : Validate an explicit cross-validation partition
- [`spatial_folds()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_folds.md)
  : Construct random, spatial-block, buffered, or kNNDM folds

## Describe spatial dependence

- [`spatial_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_covariance.md)
  : Specify a stationary isotropic spatial covariance
- [`spatial_correlation()`](https://olatunjijohnson.github.io/spCVAudit/reference/spatial_correlation.md)
  : Evaluate correlation between distinct spatial observations
- [`estimate_covariance()`](https://olatunjijohnson.github.io/spCVAudit/reference/estimate_covariance.md)
  : Estimate exponential covariance by profiling a binned variogram

## Audit and compare schemes

- [`audit_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/audit_cv.md)
  : Audit a CV partition against a prediction target
- [`distance_contrasts()`](https://olatunjijohnson.github.io/spCVAudit/reference/distance_contrasts.md)
  : Compare validation and target nearest-neighbour distances
- [`compare_cv_schemes()`](https://olatunjijohnson.github.io/spCVAudit/reference/compare_cv_schemes.md)
  : Compare candidate CV schemes using a common diagnostic
- [`plot(`*`<spcv_audit>`*`)`](https://olatunjijohnson.github.io/spCVAudit/reference/plot.spcv_audit.md)
  : Plot validation and target distance distributions

## Understand the Gaussian reference model

- [`gaussian_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/gaussian_cv.md)
  : Calculate Gaussian reference MSE, its decomposition, and a bound
- [`gaussian_target_mse()`](https://olatunjijohnson.github.io/spCVAudit/reference/gaussian_target_mse.md)
  : Evaluate conditional target MSE under a specified Gaussian
  distribution
