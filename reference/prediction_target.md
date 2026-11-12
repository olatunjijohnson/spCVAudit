# Declare the locations and weights of a prediction target

The target describes where prediction error is to be assessed. Target
response values are not required and never enter the diagnostics.
Coordinates use a common projected system; numeric matrices are assumed
to share their units.

## Usage

``` r
prediction_target(locations, weights = NULL, label = "Prediction target")
```

## Arguments

- locations:

  A numeric matrix or data frame with exactly two coordinate columns, or
  an `sf`/`sfc` object of points in a projected CRS.

- weights:

  Nonnegative target weights, normalised to sum to one. The default
  gives each target location equal weight.

- label:

  A single descriptive character string.

## Value

An object of class `spcv_target` containing locations, normalised
weights, label and, for spatial objects, CRS.

## Examples

``` r
target <- prediction_target(rbind(c(0, 0), c(1, 1)), weights = c(1, 3))
target
#> Prediction target
#> 2 locations; unequal weights
```
