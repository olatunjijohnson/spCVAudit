# Contributing

Please report reproducible problems through GitHub Issues. Include the package
version, `sessionInfo()`, a small example, and the behaviour you expected. Do not
include confidential coordinates or responses.

For changes, work on a branch and run:

```r
devtools::document()
devtools::test()
devtools::check()
pkgdown::build_site()
```

Mathematical changes should include an independently calculated example or a
relevant invariant. Keep projected-distance units and the distinction between
observed error, conditional MSE and expected CV bias explicit. New diagnostics
must document their assumptions and limitations.
