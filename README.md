| **Build Status** | **Coverage** |
|:--------------------:|:----------------:|
| [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]

# PoGO.jl: Piecewise-linear approximation for Global Optimization

The key functions provided by this package are:
- `approximate` this allows for the piecewise linear approximation of arbitrary functions. Including those that are discontinuous or otherwise non-smooth. Various methods are supported.
- `bilinear` this creates an approximation of the product of two variables. If either variable is integer or binary, then this will be modelled exactly.
- `power` this creates an approximation of one variable to the power of another.
- `interpolate` this creates a constraint forcing a set of variables to lie inside at least one set (from a group of sets), that may be defined in terms of other variables.
- `interpolate_fn` this creates one (or more) new variable(s) that is / are defined as a function of other variables at various points.
- `interpolate_points` this takes a set of multidimensional points and using multidimensional Delaunay triangulation over a subset of the dimensions, creates one (or more) new variable(s) that is / are defined by interpolation of the points over the space defined by the triangulation. If the triangulation is over all the dimensions, then rather than interpolating values for new variables, this function will define a feasible region for the variables.

[build-img]: https://github.com/adow031/PoGO.jl/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/EPOC-NZ/JuDGE.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/github/adow031/PoGO.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/adow031/PoGO.jl?branch=main

