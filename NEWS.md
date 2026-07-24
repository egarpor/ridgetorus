# ridgetorus 1.0.0

* Initial release.

# ridgetorus 1.0.1

* Fix UBSAN additional issues generated when `implicit_eq()` was called with an empty argument `theta2`.

# ridgetorus 1.0.2

* Update reference.
* Drop C++11 requirement to adhere to new CRAN policies.
* Drop `personList()` and `citEntry()`.

# ridgetorus 1.0.3

* Fix NOTEs on Rd cross-references and CRAN incoming feasibility.

# ridgetorus 1.1.0

* Add a package logo.
* Fix the moment estimator of the wrapped Cauchy dependence parameter in `fit_bwc_mm()`: use the correct Möbius tangent half-angle factor `(1 + xi) / (1 - xi)`.
* Speed up ridge estimation by vectorizing the eigenvalue check in `ridge_bvm()`/`ridge_bwc()`/`ridge_bwn()`, the wrapped normal density `d_bwn()`, and the node matching in `ridge_fourier_fit()`.
* `fit_bwn_mle()` now also returns the optimization object `opt`.
* Fix `ridge_pca()`'s validation message for `lrts`.
* Documentation fixes and new `@seealso` cross-references between the package functions.

# ridgetorus 1.1.1

* Cleaner equations in the documentation.