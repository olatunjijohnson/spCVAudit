# spCVAudit

**Compare spatial cross-validation with the places where you need
predictions.**

[![R package
checks](https://github.com/olatunjijohnson/spCVAudit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/olatunjijohnson/spCVAudit/actions/workflows/R-CMD-check.yaml)

spCVAudit compares validation-to-training distances with
target-to-training distances, adds information about spatial
correlation, and helps you assess candidate CV schemes. It implements
the nearest-neighbour contrasts and Gaussian reference calculations in
*Nearest-neighbour distances and prediction error in spatial
cross-validation*.

**[Website and
tutorials](https://olatunjijohnson.github.io/spCVAudit/)** · **[Paper
and reproducibility
code](https://github.com/olatunjijohnson/spLeakage-paper)** · **[Report
a problem](https://github.com/olatunjijohnson/spCVAudit/issues)**

## Install

``` r

install.packages("remotes")
remotes::install_github("olatunjijohnson/spCVAudit")
```

The tutorials are available on the website. To also install local copies
for `vignette(package = "spCVAudit")`, use:

``` r

install.packages(c("knitr", "rmarkdown"))
remotes::install_github("olatunjijohnson/spCVAudit", build_vignettes = TRUE)
```

Version 0.1.0 is available from GitHub. It is prepared for CRAN
submission but has not been submitted to, or accepted by, CRAN. Core
calculations use only base and recommended R packages. Install `sf` for
projected spatial objects and `CAST` for the optional kNNDM constructor.

## A first comparison

``` r

library(spCVAudit)

# Numeric coordinates must share a projected coordinate system and units.
set.seed(42)
training <- matrix(runif(160), ncol = 2)
target <- prediction_target(
  as.matrix(expand.grid(x = seq(0.05, 0.95, length.out = 10),
                       y = seq(0.05, 0.95, length.out = 10)))
)

schemes <- list(
  random = spatial_folds(training, k = 5, method = "random", seed = 1),
  block = spatial_folds(training, k = 5, method = "block", seed = 1)
)

covariance <- spatial_covariance(range = 0.2, signal = 0.8)
comparison <- compare_cv_schemes(training, schemes, target, covariance,
                                 criterion = "J2")
comparison
plot(comparison$audits$random)
```

This comparison chooses the smallest absolute `J2`. It does not use
target responses or choose the smallest observed prediction error. Every
tied scheme receives an explicit selection probability.

## What the diagnostics tell you

| Quantity | Meaning |
|----|----|
| `A` | Mean target nearest distance minus mean validation nearest distance |
| `W` | Overall difference between the two nearest-distance distributions |
| `J1` | Validation minus target mean nearest-observation correlation |
| `J2` | Validation minus target mean squared nearest-observation correlation |
| `U1`, `U2` | Unsigned discrepancies between correlation and squared-correlation distributions |

`W = 0` means the nearest-distance distributions match. `J2 = 0` means
their mean squared correlations match. Neither guarantees unbiased CV
for a model that uses all training observations.

For Gaussian prediction using **one neighbour**, with known mean and
covariance, target MSE minus expected CV MSE equals total response
variance times `J2`.
[`gaussian_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/gaussian_cv.md)
calculates the full-neighbour MSE difference, separates the
additional-neighbour contribution, and supplies its bound under positive
independent observation noise. These are model-based reference
calculations, not a correction for an arbitrary fitted prediction model.

## Learn with worked examples

- [Start with a spatial CV
  audit](https://olatunjijohnson.github.io/spCVAudit/articles/getting-started.html)
- [Compare schemes, weights and
  ties](https://olatunjijohnson.github.io/spCVAudit/articles/choosing-cv.html)
- [Understand the Gaussian identity and its
  limits](https://olatunjijohnson.github.io/spCVAudit/articles/theory-and-limits.html)
- [Audit log-cadmium predictions in the Swiss
  Jura](https://olatunjijohnson.github.io/spCVAudit/articles/jura.html)

## Relationship to spLeakage and the paper

spCVAudit is a focused successor to the spatial CV diagnostics developed
in spLeakage. It uses the paper’s names `J1` and `J2`, supports weighted
targets and explicit ties, and exposes the Gaussian decomposition and
bound. Older leakage checks and experimental error-correction functions
are not part of this API. The original spLeakage package remains
separate.

The paper repository preserves the exact spLeakage source snapshot and
fitting procedure used for its retained results. Installing spCVAudit is
not a replacement for running that fixed reproduction environment. Its
tutorials identify which calculations reproduce the theory and which
illustrate a new analysis.

## Citation

Use `citation("spCVAudit")` for the software citation. The associated
manuscript, supplementary material, data attribution and complete
reproduction instructions are in the linked paper repository. No
publication DOI is claimed before one is assigned.

Developed by Olatunji Johnson. Released under the MIT licence.
