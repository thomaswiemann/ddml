# Construct a Replicated RAL Inference Object

Creates a replicated RAL inference object from a list of `ral` objects.
Provides cross-resample aggregation for coefficients and covariance
matrices.

## Usage

``` r
ral_rep(fits, subclass = NULL, ...)
```

## Arguments

- fits:

  A list of at least 2 objects inheriting from class `"ral"`. All fits
  must share the same coefficient names and number of observations.

- subclass:

  Optional character string prepended to the class vector.

- ...:

  Additional named elements stored in the object.

## Value

An object of class `"ral_rep"` (or `c(subclass, "ral_rep")`).
