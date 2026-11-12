# Release preparation

Version **0.1.0** is available from GitHub. It has **not** been submitted to CRAN.

The local `R CMD check --as-cran` passed on Ubuntu 24.04 with R 4.6.1:
**0 errors, 0 warnings, 1 note**. The sole note identifies a new submission.
The check includes the PDF manual, HTML validation, examples, tests and all four
rebuilt vignettes. Public URLs were checked after publication.

Additional verification:

- Installation from the public GitHub repository into a fresh temporary library,
  without GitHub credentials, followed by an installed-package calculation.
- Agreement with 984 retained paper and Gaussian reference cases, with maximum
  absolute difference 1.72e-14.
- 28 desktop/mobile browser checks and 656 resolved internal documentation links.
- Automated Linux, Windows, macOS and R-development checks in
  [GitHub Actions](https://github.com/olatunjijohnson/spCVAudit/actions/workflows/R-CMD-check.yaml).

All four platform jobs passed in
[run 43322419093](https://github.com/olatunjijohnson/spCVAudit/actions/runs/43322419093)
for the implementation at commit `7365863`. Subsequent release preparation
adds installation guidance and check records without changing the R implementation.

The release assets include the source archive and a check log. The package's
maintainer should submit that archive through the
[CRAN submission form](https://cran.r-project.org/submit.html) and complete the
maintainer confirmation. Passing these checks does not guarantee CRAN acceptance.

The source under `spLeakage-paper/package/` remains a separate, frozen research
dependency. It is not the spCVAudit package submitted for review.
