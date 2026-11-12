## Submission

This is a new submission, version 0.1.0.

spCVAudit provides weighted spatial cross-validation diagnostics and Gaussian
reference calculations. The name has not previously been used on CRAN. It is a
separate package from the development package spLeakage; the associated research
repository preserves the earlier implementation used for the paper's results.

The package has no compiled code. Core functions need only base and recommended
R packages. Optional spatial and kNNDM dependencies are used conditionally. The
Jura vignette reads data distributed with gstat; no dataset is downloaded or
bundled in this package.

## R CMD check results

Local Ubuntu 24.04, R 4.6.1, `R CMD check --as-cran`, including PDF manual,
HTML validation, examples, tests and rebuilt vignettes:

0 errors | 0 warnings | 1 note

The note is:

    New submission

The public URLs pass the incoming checks. GitHub Actions passed on Linux,
Windows and macOS with R release, and on Linux with R devel (run 43322419093).
See the release validation record for the completed run and exact source version.

## Reverse dependencies

There are no known reverse dependencies because this is the initial release.
