# Validate an explicit cross-validation partition

Every observation must be tested exactly once. Training and test indices
must be nonempty, unique, disjoint, and within `1:n`. Training sets may
exclude additional observations, as in buffered CV. Repeated CV should
be audited one complete partition at a time.

## Usage

``` r
cv_folds(folds, n)
```

## Arguments

- folds:

  A list whose elements each contain integer vectors `train` and `test`;
  or a vector of fold labels of length `n`.

- n:

  Number of training observations.

## Value

A list of class `spcv_folds`, with one `train`/`test` list per fold and
an `n` attribute recording the number of observations.

## Examples

``` r
cv_folds(c(1, 2, 1, 2), n = 4)
#> [[1]]
#> [[1]]$train
#> [1] 2 4
#> 
#> [[1]]$test
#> [1] 1 3
#> 
#> 
#> [[2]]
#> [[2]]$train
#> [1] 1 3
#> 
#> [[2]]$test
#> [1] 2 4
#> 
#> 
#> attr(,"class")
#> [1] "spcv_folds" "list"      
#> attr(,"n")
#> [1] 4
```
