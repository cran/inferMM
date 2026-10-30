# inferMM 0.1.1

- Fixed the profile estimator's covariance inverse for singleton clusters by
  explicitly retaining one-by-one diagonal matrix dimensions.
- Clustered data cleaning now removes rows with a missing value in any
  clustering column, with a warning reporting the number removed. Both
  estimators, variance-model screening and reported cluster counts use the
  same retained data. Added regression tests for singleton clusters, missing
  and multi-column IDs, unused factor levels, NaN IDs and refitting.
- Strengthened projected-cluster numerical diagnostics: non-finite tail grid
  evaluations no longer discard an otherwise verified interior root, and
  refinement warnings are consolidated and recorded. Both the scalar and
  joint scaled scores are checked. Failures report their stage, reason and
  numerical diagnostics without returning estimates or intervals; there is
  no automatic switch to the profile method. Non-finite variance shapes,
  unrepresentable scales and
  non-positive joint covariance diagonals fail explicitly. Boundary variance
  moments and incomplete searches have distinct structured warning reasons.
  Legacy profile iteration-limit warnings now report the last parameter
  change and explicitly mark the fit/intervals as unconverged.
- Added opt-in `cluster_mm(method = "projected")`, an MM-specific
  Ma-Genton-style conditional projected-score construction with the identity
  instrument. It retains conditionally Gaussian residuals but does not specify
  the random-effect density. A scalar Km score is followed by explicit moment
  updates and a stacked independent-cluster sandwich covariance, not the
  legacy mean-information approximation or a general efficiency guarantee.
  Refits retain this method; pointwise mean-confidence intervals are available,
  but distribution-free prediction quantiles are not. The legacy profile
  method remains the default. Root, boundary and covariance checks stop invalid
  fits without returning estimates or intervals.
- Single-curve profile fits now verify standardized mean-score residuals,
  compare detected acceptable roots by weighted residual sum of squares,
  and reject boundary, non-positive, non-finite or numerically degenerate
  solutions. Score minimization supplies candidates, not automatic success.
  Single-curve covariance uses scaled unregularized information and errors on
  singular or poorly conditioned matrices instead of silently adding ridge.
  Diagnostics are stored in `convergence`; failures carry class
  `inferMM_profile_error`, also inheriting from `inferMM_convergence_error`.
  No estimates or intervals are returned after a validation error. Automatic
  variance selection records failed candidates and selects only successful
  fits. Screening/grouping and the higher-moment profile pilot share the new
  root checks; the separate clustered solver is unchanged.
- Expanded higher-moment convergence errors to identify the stage, stopping
  reason, iterations, scaled score norm, tolerance, and continuation progress.
  Failures can be caught as `inferMM_convergence_error`; no fit or confidence
  intervals are returned. Non-Gaussian continuation failures suggest a manual
  retry with normal working moments and the same variance function, without
  automatic fallback or a guarantee of convergence or non-Gaussian efficiency.
- Added the opt-in `fit_mm_moments()` estimator for independent curves. It
  jointly solves MM mean/variance estimating equations weighted using third
  and fourth conditional moments, without calling `nls()`.
- Working standardized moments can be estimated from pilot residuals, set
  using centered-gamma, Student-t or normal moments, or supplied as numbers
  or concentration-dependent functions. Fixed user-defined variance shapes
  are supported. This is not a parametric error-density likelihood fit.
- Added full joint working-information (default) and optional sandwich covariance,
  Wald intervals, pairs-bootstrap refits, and pointwise mean-confidence bands.
  The nuisance scale is included before extracting the mean covariance.
- Added continuation/root diagnostics and explicit failure for unconverged,
  boundary, or ill-conditioned solutions. Higher moments do not identify
  predictive quantiles, so this estimator does not report prediction bands.
- Added standardized Student-t simulation errors with degrees of freedom
  greater than four. The existing centered-gamma simulation remains available.
- Existing `fit_mm()`, variance screening/grouping, and clustered estimators
  retain their defaults. The new estimator does not claim automatic finite-
  sample coverage, robustness to an incorrect variance shape, or dominance
  over all weighted or transformation-based methods.

# inferMM 0.1.0

- `fit_mm()` accepts a fixed user-defined positive function of concentration;
  `variance_function()` can attach a descriptive label to it.
- `screen_mm()` and `group_mm()` accept named `variance_functions` lists for
  comparison with built-in candidates, including custom-only screening.
- Custom functions are validated without silently clipping invalid values and
  retained for prediction, plotting, and wild/pairs bootstrap refits.
- Fixed custom shapes retain the scalar Km profile-score calculation and
  closed-form Vmax and variance-scale updates. Shape parameters are not fitted;
  inference and information criteria condition on the supplied shape.
- The clustered fitting algorithm and built-in variance choices are unchanged.
  User-defined shapes do not implement higher-moment efficient scores.
- Corrected the help pages to call `set.seed()` outside package functions;
  random-number handling is unchanged.

# inferMM 0.0.3

- Removed the `nls` fitting path so all single-curve fits now use the manuscript's Michaelis-Menten profile-score estimator.
- Updated clustered initialization to use the pooled single-curve profile-score fit without an `nls` fallback.
- Refreshed package documentation to reflect the unified profile-score implementation.

# inferMM 0.0.1

- Initial development version.
- Added variance-aware Michaelis-Menten fitting and screening functions.
- Added grouped workflows, simulation helpers, tests, and demo data.
