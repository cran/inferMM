# inferMM

`inferMM` provides variance-aware Michaelis-Menten estimation and inference for
enzyme-kinetic data with concentration-dependent heteroscedasticity. The
package supports single-curve analyses, grouped workflows, and clustered
repeated-measurement settings through Michaelis-Menten-specific score-based
estimating routines and companion S3 summaries, plotting, and reporting tools.

The package is designed around a compact workflow:

- fit one curve with `fit_mm()`
- optionally fit independent curves with higher moments using `fit_mm_moments()`
- compare working variance models with `screen_mm()`
- repeat that workflow across many enzymes with `group_mm()`
- fit repeated or clustered assays with `cluster_mm()`
- summarize and plot fitted objects with standard S3 methods plus `report_mm()`

## Installation

```r
install.packages("inferMM")
```

## Bundled data

The package ships with two demo datasets.

- `sdl_demo`: self-driving laboratory Michaelis-Menten panel
- `alves_demo`: filtered soil exoenzyme kinetics panel from Alves et al. (2021)

```r
library(inferMM)
data(sdl_demo)
data(alves_demo)
```

## Minimal example

```r
library(inferMM)

one_curve <- subset(sdl_demo, enzyme == "1111")

fit <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "sqrt"
)

summary(fit)
confint(fit)
```

## Screening working variance models

Single-curve fits verify standardized mean-score residuals, positive interior
parameters, residual variance, and the condition of the scaled information
matrix before returning estimates and intervals. The profile search compares
all detected acceptable roots by weighted residual sum of squares rather than
selecting the first bracket. Score minimization supplies candidates only;
it is not automatically treated as convergence. Single-curve Wald covariance
uses an unregularized scaled information inverse, not a silent ridge correction.
No numerical search guarantees detection of every root or a global minimum.

Diagnostics are stored in `fit$convergence`. Errors carry class
`inferMM_profile_error` (also `inferMM_convergence_error`), with stage, reason,
search bounds, standardized score norm, tolerance and candidate count; a
failed fit returns no estimates or confidence intervals. Check `km_bounds`,
substrate coverage, outliers, the variance shape and identifiability.
Screening/grouping records failed curves without fabricating results.
Automatic variance selection retains failed candidates with `success = 0`
and their error messages in `variance_candidate_table`, selecting only
successful candidates and stopping if all fail. The higher-moment profile
pilot uses the same root checks; the clustered solver retains its own rules.

```r
screen <- screen_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  quiet = TRUE
)

screen$table[, c("model", "selected_model", "quasi_aic", "quasi_bic", "rmse")]
```

## User-defined fixed variance functions (0.1.0)

Pass a vectorized positive function directly as `variance`. A descriptive label
is optional and can be attached with `variance_function()`:

```r
h <- variance_function(
  function(S) 1 + S / 20,
  label = "Var = gamma (1 + S/20)"
)
custom_fit <- fit_mm(one_curve$s_uM, one_curve$v_uM_per_min, variance = h)
summary(custom_fit)
confint(custom_fit)
predict(custom_fit, newdata = c(0, 10, 40), interval = "prediction")

custom_shapes <- list(
  linear_floor = function(S) 1 + S / 20,
  saturation = function(S) 1 + 9 * S / (20 + S)
)
custom_screen <- screen_mm(
  one_curve$s_uM, one_curve$v_uM_per_min,
  variance_models = c("constant", "sqrt"),
  variance_functions = custom_shapes, quiet = TRUE
)
custom_screen$table[, c("model", "quasi_aic", "quasi_bic", "success")]
```

Use `variance_models = character(0)` to compare custom functions only. The same
`variance_functions = custom_shapes` argument works in `group_mm()`. Candidate
names identify fits and ranking rows; invalid functions are recorded as failed
candidates instead of preventing the other candidates from being fitted.

The function represents **variance shape**, not standard deviation or weights.
It must be deterministic and pointwise, returning a numeric vector of the same
length as its input, with finite strictly positive values at every fitting and
prediction concentration. Invalid custom values are rejected, not clipped.
Functions must also work on new grids and bootstrap samples; do not return a
cached vector or use the input's mean, range, or ordering to define the shape.
Use a local closure to capture fixed external constants if necessary.

For each fixed shape, the existing scalar Km profile-score fit and closed-form
Vmax and scale updates are retained. Custom functions are preserved for
prediction, plotting, and bootstrap refits. Internal normalization does not
change the supplied variance model; the reported scale refers to the original
shape. Numerical difficulties can still arise for ill-scaled shapes or weakly
identified curves.

This interface **does not estimate extra shape parameters**. Standard errors,
intervals, and the three-parameter working information criteria condition on the
supplied shape and omit external shape-estimation or selection uncertainty.
Prediction bands retain the working-normal approximation. This feature does not
estimate higher conditional moments, prove efficiency gains over same-weight
WNLS, or extend `cluster_mm()` to user-defined functions. The separate
`fit_mm_moments()` estimator below adds a higher-moment route.

User-defined variance functions require inferMM 0.1.0 or later. Check the
installed version with `packageVersion("inferMM")`.

## Higher-moment estimation for independent curves (0.1.1)

`fit_mm_moments()` uses conditional third and fourth moments to weight joint
mean/variance estimating equations, following Kim and Ma (2012), Section 3,
<https://doi.org/10.1007/s10463-011-0332-y>. It is an opt-in estimator, not a
replacement for `fit_mm()` and not a gamma or Student-t likelihood fit.

```r
set.seed(19)
dat <- simulate_mm_data(n = 200, error = "skewed", skew_shape = 4)
h <- function(S) 1 + 9 * S / (20 + S)

# Estimate constant standardized skewness/kurtosis from pilot residuals.
fit_hm <- fit_mm_moments(dat$x, dat$y, variance = h)
summary(fit_hm)
confint(fit_hm, level = 0.95)
fit_hm$vcov_joint  # Vmax, Km, and the nuisance variance scale together

# Alternatively, supply fixed working moments of centered gamma errors.
fit_gamma <- fit_mm_moments(dat$x, dat$y, variance = h,
                            moments = "gamma", shape = 4)
predict(fit_gamma, c(10, 20, 80), interval = "confidence")

# Concentration-dependent standardized moments are also possible.
fit_custom_moments <- fit_mm_moments(
  dat$x, dat$y, variance = h, moments = "custom",
  skewness = function(S) 0.8 + 0.2 * S / (20 + S),
  kurtosis = function(S) 3 + 1.5 * (0.8 + 0.2 * S / (20 + S))^2
)

set.seed(20)
confint(fit_hm, method = "bootstrap", B = 99)  # pairs; re-estimates moments
```

The default `covariance = "model"` uses joint inverse working information
and needs correct higher moments for its interpretation. Optional
`covariance = "sandwich"` differentiates the full three-parameter
estimating equations. The
leading mean-parameter block is extracted *after* accounting for the nuisance
variance scale. `moments = "student_t", df = 10` supplies symmetric heavy-tail
moments; `simulate_mm_data(error = "student_t", df = 10)` generates matching
standardized errors. Normal or constant symmetric working moments yield the
same mean point estimates as the profile equations, not an efficiency gain.

The new estimator uses a joint numerical solve with continuation from the
profile pilot; it does **not** retain the scalar/closed-form updates of the
default estimator. It reports a score norm and errors on failed/boundary roots
or singular covariance. It does not silently return an old fit as a new one.
Convergence errors identify the fitting stage and reason, total iterations,
scaled score norm, tolerance, and completed continuation fraction where relevant.
They suggest checks for the data, variance shape, and working moments, or a
larger iteration limit when appropriate. They can be handled with
`tryCatch(..., inferMM_convergence_error = function(e) ...)`; the condition
also carries `stage`, `reason`, `moments`, `iterations`, `score_norm`,
`tolerance`, and `progress` fields. Failed fits return no estimates or intervals.
Failures during non-Gaussian higher-moment continuation also suggest manually
retrying with `moments = "normal"` while keeping the same variance function.
This changes the working moments; it guarantees neither convergence nor
efficiency under the original non-Gaussian moments. There is no automatic fallback.
Only pointwise mean-confidence bands are available: four moments alone do not
specify predictive quantiles. Working-Gaussian AIC/BIC are not reported for
these non-likelihood roots. `screen_mm()`, `group_mm()`, and `cluster_mm()`
continue to use their original estimators, not this new route.

Efficiency in the conditional mean/variance model requires suitable moment
estimation and regularity. The default estimates *constant standardized*
moments, not arbitrary conditional moment functions. Correct mean/variance
remain important: a sandwich is not a remedy for bias caused by an incorrect
variance shape in these squared-residual equations. Small-sample performance,
especially high-moment estimation, must be checked; no universal efficiency
or 95% coverage guarantee is made. See `?fit_mm_moments` for the equations,
solver, assumptions, and covariance details.

## Grouped analyses

```r
grouped <- group_mm(
  data = sdl_demo,
  s = "s_uM",
  v = "v_uM_per_min",
  groups = "enzyme",
  variance_models = c("constant", "log", "sqrt", "cuberoot"),
  quiet = TRUE
)

grouped$comparison$best_by_group[
  , c("group_label", "model", "selected_model", "quasi_aic", "quasi_bic", "rmse")
]
```

## Clustered analyses

```r
cluster_fit <- cluster_mm(
  data = subset(alves_demo, enzyme == "BG"),
  s = "substrate_conc",
  v = "activity",
  cluster = "core",
  variance = "sqrt"
)

summary(cluster_fit)
confint(cluster_fit)
```

Both clustered estimators remove rows with a missing value in any clustering
column and issue a warning with the number removed. Cluster counts, cluster
sizes, fitted values and residuals then describe only the retained data.
Unused factor levels do not create extra clusters. The default profile method
supports singleton clusters alongside repeated-measurement clusters; the
projected method continues to require at least two readings per cluster.

For sparse clustered fits, default interval reporting is intentionally cautious:
printed summaries may suppress intervals, and bootstrap intervals should be read
as sensitivity analyses rather than routine default inference.

### Unspecified random-effect density

The opt-in `method = "projected"` adds an MM-specific adaptation of the
conditional-centering construction of Ma and Genton (2010), Section 3.2,
<https://doi.org/10.1111/j.1467-9868.2010.00741.x>. It leaves the distribution
of the cluster random effect unspecified while retaining conditionally
independent Gaussian residuals with the chosen fixed variance shape. This is
different from the non-Gaussian residual analysis in `fit_mm_moments()`.

```r
set.seed(351)
substrate <- c(2, 5, 10, 20, 50, 100, 200)
n_clusters <- 60L
b <- 12 * (rgamma(n_clusters, shape = 2) - 2) / sqrt(2)
dat <- data.frame(
  cluster = rep(seq_len(n_clusters), each = length(substrate)),
  S = rep(substrate, n_clusters)
)
dat$Y <- (100 + b[dat$cluster]) * dat$S / (20 + dat$S) +
  rnorm(nrow(dat), sd = sqrt(2 * sqrt(dat$S)))

fit_projected <- cluster_mm(
  dat, s = "S", v = "Y", cluster = "cluster",
  variance = "sqrt", method = "projected"
)
summary(fit_projected)
confint(fit_projected)
predict(fit_projected, c(10, 20, 80), interval = "confidence")
```

For each candidate Km, the new path projects observations onto the cluster's
MM loading, leaving a residual orthogonal to that loading. A scalar score
with the identity instrument estimates Km; Vmax, the residual scale and the
random-effect variance moment then have explicit updates. The identity
instrument belongs to a consistent estimating-equation class under the stated
model and identification conditions, but is **not generally the efficient
instrument**. This implementation does not claim universal efficiency or
better finite-sample performance than the legacy estimator.

The new path uses the **stacked independent-cluster sandwich**, including
mean and variance-moment equations, for standard errors and 95% Wald intervals.
It does not reuse the legacy mean-information covariance. At least five
independent clusters and two readings per cluster are required; larger cluster
counts are important for asymptotic inference. The variance shape must be
correct, and finite fourth random-effect moments are assumed for the joint
covariance. Intervals condition on the supplied or selected variance shape.
Quasi-AIC/BIC remain working-normal marginal screening scores, not likelihood
information criteria for the unspecified random-effect density.
Invalid roots, boundary roots, degenerate variance and singular derivatives
stop with `inferMM_cluster_projected_error`; no fallback or ridge is applied.
`km_bounds` can adjust the positive Km search interval for this path.

Failure conditions include `stage`, `reason`, and `diagnostics`, accessible
with `tryCatch(..., inferMM_cluster_projected_error = function(e) e)`.
Root failures report finite grid counts, candidate/accepted-root counts,
the smallest grid score residual and the associated Km; covariance failures
report available conditioning diagnostics. No fitted object or interval is
returned on failure, and no automatic switch to `method = "profile"` occurs.
A non-finite tail point does not by itself invalidate a verified interior
root. An incomplete search produces `inferMM_cluster_projected_warning`
with reason `"incomplete_search"`, stored in `fit$convergence$search_warning`.
The same warning class with reason `"variance_component_boundary"` reports
a truncated non-positive variance moment; it is not mean-score failure.
Inspect `fit$convergence$joint_score_norm`, `$search`, and `$bread_rcond`.
Passing these checks is not a guarantee of global optimality or adequate
information from a narrow concentration range.

Only pointwise **mean-confidence bands** are provided, not prediction bands:
the random-effect density has not been specified. A negative random-effect
variance moment is reported as zero with a warning; that truncation affects
marginal variance reporting, not mean estimation or its covariance. Bootstrap
refits preserve the projected method. Existing calls without `method` retain
the legacy profile/moment algorithm and its covariance.

## Reporting and plotting

```r
report_mm(fit, interval_type = "confidence")
plot(grouped, interval_type = "confidence")
predict(fit, newdata = seq(0, 80, length.out = 6), interval = "prediction")
```

## Repository contents

- `R/`: package source code
- `man/`: function documentation
- `data/`: bundled `.rda` demo data
- `inst/extdata/`: raw CSV mirrors of the demo data
- `vignettes/`: end-to-end workflow vignette
- `tests/`: unit tests

This package source is focused on the reusable software interface, bundled
example data, documentation, and tests. Manuscript-specific reproducibility
materials can be distributed separately from the package source when needed.
