# Compare schemes, weights and ties

A selection rule compares candidate CV schemes for the **same target and
final prediction model**. It seeks a useful estimate of that model’s
target MSE. It does not seek the smallest CV score.

## Specify candidates and a weighted target

``` r

set.seed(101)
training <- matrix(runif(120), ncol = 2)
q <- as.matrix(expand.grid(x = seq(0.1, 1.3, length.out = 8),
                          y = seq(0.1, 0.9, length.out = 8)))
target <- prediction_target(q, weights = ifelse(q[, 1] > 1, 2, 1),
                            label = "Target with extra weight beyond the sampled region")
schemes <- list(
  random = spatial_folds(training, 5, "random", seed = 4),
  block = spatial_folds(training, 5, "block", seed = 4),
  buffered = spatial_folds(training, 5, "buffered", buffer = 0.04, seed = 4)
)
model <- spatial_covariance(range = 0.2, signal = 0.8)
```

The buffer removes training observations within 0.04 coordinate units of
a test observation. If this empties a training set, the constructor
reports an error; reduce the buffer or reconsider the partition.
Matching seeds gives block and buffered schemes the same test groups.

## Compare two selection criteria

``` r

by_W <- compare_cv_schemes(training, schemes, target, model, criterion = "W")
by_J2 <- compare_cv_schemes(training, schemes, target, model, criterion = "J2")
by_W
#> CV comparison by W 
#>    scheme         W         J2 criterion_value selection_probability
#>    random 0.1005988  0.1171102       0.1005988                   0.0
#>     block 0.0730827 -0.1249033       0.0730827                   0.5
#>  buffered 0.0730827 -0.1249033       0.0730827                   0.5
by_J2
#> CV comparison by |J2| 
#>    scheme         W         J2 criterion_value selection_probability
#>    random 0.1005988  0.1171102       0.1171102                     1
#>     block 0.0730827 -0.1249033       0.1249033                     0
#>  buffered 0.0730827 -0.1249033       0.1249033                     0
```

`W` measures a distance-distribution discrepancy. `abs(J2)` compares
average squared correlation and can allow cancellation. They need not
select the same scheme. A single common covariance specification is used
for all candidates.

If covariance is estimated from observed training responses, pass an
`estimate_covariance(training, response)` result in place of `model`.
Selection then depends on those responses. The covariance used to make
each held-out prediction should still be estimated from that fold’s
training data only.

## Add kNNDM when appropriate

kNNDM is a candidate CV construction, whereas `W` and `J2` are criteria
that can evaluate it alongside the other candidates. The optional
constructor uses the implementation in **CAST**. It currently constructs
folds for equally weighted target points. For a weighted target you may
audit externally supplied folds, but the CAST construction itself does
not optimise those weights.

``` r

install.packages(c("sf", "CAST"))
equal_target <- prediction_target(q)
schemes$knndm <- spatial_folds(training, 5, "knndm",
                              target = equal_target, seed = 4)
compare_cv_schemes(training, schemes, equal_target, model, criterion = "J2")
```

## Evaluate a rule when the generating distribution is known

In a simulation,
[`gaussian_cv()`](https://olatunjijohnson.github.io/spCVAudit/reference/gaussian_cv.md)
can calculate expected CV MSE and target MSE without simulating
responses. This example assumes the known-mean Gaussian model with the
covariance specified above. All schemes share the same target MSE.

``` r

reference <- lapply(schemes, function(f) gaussian_cv(training, f, target, model))
evaluation <- data.frame(
  scheme = names(schemes),
  target_mse = vapply(reference, `[[`, numeric(1), "target_mse"),
  cv_mse = vapply(reference, `[[`, numeric(1), "cv_mse"),
  absolute_bias = abs(vapply(reference, `[[`, numeric(1), "optimism"))
)
evaluation
#>            scheme target_mse    cv_mse absolute_bias
#> random     random  0.7183768 0.5756394     0.1427374
#> block       block  0.7183768 0.8796564     0.1612796
#> buffered buffered  0.7183768 0.8796564     0.1612796
```

The best scheme by `absolute_bias` is the one whose **expected CV MSE is
closest to target MSE**. Target MSE is used here only for evaluation,
after the diagnostic choices were made. In an application it is not
usually known.

## Handle ties without using target errors

``` r

same <- list(first = schemes$random, second = schemes$random)
tied <- compare_cv_schemes(training, same, target, model)
tied$table[, c("scheme", "criterion_value", "selection_probability")]
#>   scheme criterion_value selection_probability
#> 1  first       0.1171102                   0.5
#> 2 second       0.1171102                   0.5
```

If two schemes tie, a rule that chooses uniformly between them has
expected evaluation loss equal to the average of their **individual
losses**. For example, signed estimation errors $`-0.1`$ and $`0.1`$
give average absolute error 0.1. Averaging the signed errors first would
give zero and would assess a different estimator.

``` r

sum(by_J2$table$selection_probability * evaluation$absolute_bias)
#> [1] 0.1427374
```

The default absolute tie tolerance is $`10^{-12}`$ in criterion units.
Report and check that tolerance when very weak correlation makes all
contrasts tiny. With zero spatial signal, all `J2` values are zero and
all candidates tie.
