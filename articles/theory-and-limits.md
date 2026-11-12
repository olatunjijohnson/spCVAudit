# The Gaussian identity and what other neighbours change

The reference calculations assume jointly Gaussian responses, known mean
and covariance, and fixed training locations, target weights and CV
folds chosen independently of responses. They average over repeated
response datasets at those locations. They do not equate an individual
observed CV score with target MSE.

## Prediction using one neighbour

Let the total response variance be $`V`$ and correlation with the
nearest available training response be $`\rho(d)`$. With known mean
$`\mu`$, the Gaussian conditional mean based on that response $`Y_1`$ is

``` math
\widehat Y(t)=\mu+\rho(d)(Y_1-\mu).
```

Its MSE is the conditional variance, $`V\{1-\rho(d)^2\}`$. Averaging
separately over target and validation locations and subtracting gives

``` math
O_1=R_Q-C=V\left[\mathrm{E}_C\rho(D_C)^2-\mathrm{E}_Q\rho(D_Q)^2\right]=VJ_2.
```

Therefore minimising $`|J_2|`$ minimises **absolute expected CV bias**
in this model. It does not minimise prediction MSE itself.

## Prediction using all observations

Let $`e(t,T)`$ be the MSE reduction from conditioning on all training
responses in $`T`$. The other observations provide the additional
reduction

``` math
a(t,T)=e(t,T)-V\rho[D(t,T)]^2\geq0.
```

The full difference is

``` math
O_{\mathrm{full}}=VJ_2+\overline a_C-\mathrm{E}_Q a(t,S).
```

The final two terms compare the benefits of the other neighbours during
validation and at the target. They need not cancel, and can change which
CV scheme has the smallest absolute bias.

## Reproduce the paper’s counterexample

``` r

training <- cbind(c(-1, 0, 1), 0)
target <- prediction_target(matrix(c(2, 0), nrow = 1))
model <- spatial_covariance(range = 1, signal = 0.8, variance = 1)
folds <- cv_folds(1:3, n = 3)

audit_cv(training, folds, target, model)
#> Spatial CV audit
#>             A             W            J1            J2            U1 
#>  0.000000e+00  0.000000e+00  0.000000e+00 -1.387779e-17  0.000000e+00 
#>            U2 
#>  0.000000e+00 
#> Diagnostics compare prediction conditions; they do not certify CV accuracy.
result <- gaussian_cv(training, folds, target, model)
result
#> Gaussian reference calculation
#>            target_mse                cv_mse              optimism 
#>           0.912868743           0.889812660           0.023056082 
#>         one_neighbour additional_neighbours                 lower 
#>           0.000000000           0.023056082          -0.002661682 
#>                 upper 
#>           0.116354041
```

Every nearest-training distance is one, so `W` and `J2` are zero.
Nevertheless, target MSE exceeds expected CV MSE by approximately
0.02306. All of this difference comes from `additional_neighbours`.

When the observation at zero is held out, training observations lie one
unit away on both sides. At the target location two, the training
observations lie one, two and three units away, all on the same side.
Matching only the closest distance omits this difference in available
information.

``` r

result$validation[, c("nearest", "additional", "bound")]
#>      nearest   additional       bound
#> 1 0.08661458 0.0005133434 0.002344402
#> 2 0.08661458 0.0696915894 0.344373320
#> 3 0.08661458 0.0005133434 0.002344402
result$target[, c("nearest", "additional", "bound")]
#>      nearest   additional       bound
#> 1 0.08661458 0.0005166762 0.002661682
```

## Bound the contribution

With independent observation noise variance $`\tau^2>0`$, write
$`\widetilde b_R`$ for the remaining observations’ target covariances
after removing their association through the nearest observation. If the
spatial variance is $`\sigma^2`$, the additional reduction satisfies

``` math
0\leq a(t,T)\leq B(t,T)=\min\left\{
\frac{\|\widetilde b_R\|^2}{\tau^2},\ \sigma^2-\frac{b_1^2}{V}\right\}.
```

Consequently,

``` math
VJ_2-\mathrm{E}_Q B(t,S)\leq O_{\mathrm{full}}\leq VJ_2+\overline B_C.
```

``` r

c(lower = result$lower, optimism = result$optimism, upper = result$upper)
#>        lower     optimism        upper 
#> -0.002661682  0.023056082  0.116354041
```

A positive lower bound guarantees underestimation of target MSE on
average under the specified model. A negative upper bound guarantees
overestimation. An interval containing zero does not establish a
direction. These are deterministic bounds, not confidence intervals;
they can be wide and require more covariance information than nearest
distances alone. The complete proof is in the paper’s supplementary
material.

## Conditional MSE for an observed training dataset

Conditional on training responses, the MSE of an arbitrary supplied
prediction is

``` math
\operatorname{Var}(Y(t)\mid Y_S)
+\{\mathrm{E}(Y(t)\mid Y_S)-\widehat Y(t)\}^2.
```

Only when the prediction equals the conditional mean does MSE equal
conditional variance. The following reference calculation makes the
distinction explicit:

``` r

y <- c(0.2, -0.1, 0.5)
optimal <- gaussian_target_mse(training, y, target, model)
shifted <- gaussian_target_mse(training, y, target, model,
                              predictions = optimal$conditional_mean + 0.5)
c(conditional_mean_prediction = optimal$mse,
  shifted_prediction = shifted$mse,
  increase = shifted$mse - optimal$mse)
#> conditional_mean_prediction          shifted_prediction 
#>                   0.9128687                   1.1628687 
#>                    increase 
#>                   0.2500000
```

The increase is $`0.5^2=0.25`$. In simulations the generating
distribution is known. For real observations, replacing it with an
estimated covariance gives a model-based calculation, not an exact
measurement of prediction error.
