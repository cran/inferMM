\name{fit_mm}
\alias{fit_mm}
\alias{coef.fit_mm}
\alias{confint.fit_mm}
\alias{plot.fit_mm}
\alias{predict.fit_mm}
\alias{print.fit_mm}
\alias{summary.fit_mm}
\alias{print.summary.fit_mm}
\alias{vcov.fit_mm}
\title{Fit a Variance-Aware Michaelis-Menten Model}
\usage{
fit_mm(x, y,
  variance = c("log", "sqrt", "cuberoot", "constant", "power", "auto"),
  power = NULL, power_selection = c("quasi_aic", "quasi_bic"),
  power_grid = default_power_grid(), interval_level = 0.95,
  km_bounds = NULL, allow_zero_substrate = TRUE)

\method{coef}{fit_mm}(object, ...)

\method{vcov}{fit_mm}(object, ...)

\method{confint}{fit_mm}(object, parm = c("Vmax", "Km"),
  level = object$interval_level, method = c("wald", "bootstrap"),
  B = 399L, bootstrap_type = c("wild", "pairs"),
  bootstrap_ci = c("studentized", "basic", "percentile"),
  wild_weights = c("mammen", "rademacher"), seed = NULL, ...)

\method{print}{fit_mm}(x, digits = max(3L, getOption("digits") - 2L), ...)

\method{summary}{fit_mm}(object, level = object$interval_level,
  method = c("wald", "bootstrap"), B = 399L,
  bootstrap_type = c("wild", "pairs"),
  bootstrap_ci = c("studentized", "basic", "percentile"),
  wild_weights = c("mammen", "rademacher"), seed = NULL, ...)

\method{print}{summary.fit_mm}(x,
  digits = max(3L, getOption("digits") - 2L), ...)

\method{predict}{fit_mm}(object, newdata = NULL, se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = object$interval_level, ...)

\method{plot}{fit_mm}(x, interval = TRUE,
  interval_type = c("prediction", "confidence"),
  level = x$interval_level, n_points = 200,
  xlab = "Substrate concentration", ylab = "Reaction velocity",
  main = NULL, pch = 19, col.points = "#1f78b4",
  col.line = "#222222",
  col.band = grDevices::adjustcolor("#9ecae1", alpha.f = 0.45), ...)
}
\arguments{
  \item{x}{Numeric vector of substrate concentrations. For the plot method, a
  fitted object of class \code{"fit_mm"}.}
  \item{y}{Numeric vector of observed reaction velocities.}
  \item{variance}{A built-in working variance model name or a fixed, vectorized
  function of concentration returning strictly positive values. Use
  \code{\link{variance_function}} to attach a label to a custom function.
  Use \code{"power"} with a fixed exponent or \code{"auto"} to compare
  \code{log(S + 1)} against a grid of power functions.}
  \item{power}{Optional exponent used when \code{variance = "power"}.}
  \item{power_selection}{Criterion used when \code{variance = "auto"}:
  \code{"quasi_aic"} or \code{"quasi_bic"}.}
  \item{power_grid}{Candidate exponents considered when
  \code{variance = "auto"}.}
  \item{interval_level}{Confidence level used for interval estimation.}
  \item{km_bounds}{Numeric vector of length two giving the search bounds for
  \eqn{K_m}.}
  \item{allow_zero_substrate}{Logical; if \code{FALSE}, substrate concentrations
  must be strictly positive.}
  \item{object}{A fitted \code{"fit_mm"} object.}
  \item{parm}{Parameter names for \code{confint()}.}
  \item{level}{Confidence level for intervals.}
  \item{method}{Interval construction method: \code{"wald"} or
  \code{"bootstrap"}.}
  \item{B}{Number of bootstrap resamples when \code{method = "bootstrap"}.}
  \item{bootstrap_type}{Bootstrap scheme: \code{"wild"} or \code{"pairs"}.}
  \item{bootstrap_ci}{Bootstrap interval calibration used when
  \code{method = "bootstrap"}: \code{"studentized"}, \code{"basic"}, or
  \code{"percentile"}.}
  \item{wild_weights}{Wild-bootstrap multipliers used when
  \code{bootstrap_type = "wild"}: \code{"mammen"} or
  \code{"rademacher"}.}
  \item{seed}{Deprecated and ignored. Call \code{set.seed()} before
  bootstrap-based summaries or intervals if reproducibility is needed.}
  \item{digits}{Number of printed digits.}
  \item{newdata}{Optional numeric vector of new substrate concentrations.}
  \item{se.fit}{Logical; if \code{TRUE}, return standard errors for mean
  predictions.}
  \item{interval}{Prediction type: none, confidence, or prediction.}
  \item{interval_type}{For the plot method, draw either a prediction or
  confidence band when \code{interval = TRUE}.}
  \item{n_points}{Number of grid points used for plotting the fitted curve.}
  \item{xlab,ylab,main,pch,col.points,col.line,col.band}{Graphical parameters
  for the plot method.}
  \item{...}{Additional arguments passed to S3 methods.}
}
\value{
The main fitting function returns an object of class \code{"fit_mm"}.
Its \code{convergence} component records the standardized score norm,
tolerance, search method, candidate roots, and scaled information condition.
\code{predict()} returns either a numeric vector or a data frame containing
prediction intervals.
}
\description{
Estimate Michaelis-Menten parameters under a constant or
concentration-dependent working variance specification.
}
\details{
Built-in working variance models are \code{"constant"}, \code{"log"},
\code{"sqrt"}, \code{"cuberoot"}, \code{"power"}, and \code{"auto"}.
All single-curve fits use the Michaelis-Menten profile-score
estimator from the companion methodology work, with the constant-variance case
treated as the special case \eqn{h(S) \equiv 1}. Confidence intervals can be
computed either from the default asymptotic Wald approximation or from a
bootstrap refitting scheme. The bootstrap option supports percentile, basic,
and studentized calibration, together with wild or pairs resampling. The
default printed summaries suppress interval output when the fitted curve has
fewer than 20 observations or fewer than 6 distinct substrate concentrations;
in that setting, \code{confint()} remains available as an explicit sensitivity
analysis tool and issues a warning. The automatic option compares
\code{log(S + 1)} with a user-supplied grid of \eqn{S^p} working variance
functions and keeps the best fit according to quasi-AIC or quasi-BIC.

The profile search examines sign-changing brackets on a 400-point log grid
and verifies the mean score equations after solving. The diagnostic projects
the \eqn{K_m} score direction off the eliminated \eqn{V_{\max}} direction
and standardizes its residual inner product by both weighted norms; the
\eqn{V_{\max}} equation is also checked. The maximum absolute standardized
score must not exceed \eqn{10^{-7}}. If no bracketed root is verified,
minimization of the squared standardized profile score supplies additional
candidates, not automatically accepted estimates. Among detected positive
interior roots passing these checks, the smallest weighted residual sum of
squares determines the fitted root. A finite grid does not guarantee finding
every root or a global minimum.

Candidates within \eqn{10^{-6}} of either log search bound are rejected,
as are non-positive or non-finite parameters and numerically zero residual
variance. The mean information is scaled by its gradient-column magnitudes;
a reciprocal condition number below \eqn{10^{-12}} stops covariance
calculation. No ridge term is added to obtain the single-curve Wald covariance.
These checks do not guarantee identifiability or accurate finite-sample coverage.

Profile validation errors have class \code{inferMM_profile_error}, also
inheriting from \code{inferMM_convergence_error}. They identify the failure
stage, reason, search bounds, score norm, tolerance and candidate count, with
information conditioning where relevant. No fit or confidence intervals are
returned on these errors. Use \code{tryCatch()} to inspect the condition.
Automatic variance selection records unsuccessful candidates in
\code{variance_candidate_table} using \code{success = 0} and an
\code{error_message}; they cannot be selected. If every candidate fails,
selection stops. Screening and grouped workflows likewise record failures
without reporting estimates or intervals for those curves. The profile pilot
of \code{\link{fit_mm_moments}} uses the same root-validation checks;
the later joint higher-moment solve has its own diagnostics. These new
single-curve checks do not replace the separate clustered solver.

A custom function specifies \eqn{h(S)} in
\eqn{\operatorname{Var}(Y\mid S)=\gamma h(S)}; see
\code{\link{variance_function}} for its input/output contract. Only \eqn{K_m}
requires numerical search, followed by closed-form updates for \eqn{V_{\max}}
and \eqn{\gamma}. The stored \code{variance_fun} is reused for new-data
prediction and wild or pairs bootstrap refitting. Custom shapes are normalized
internally for numerical scaling, without clipping; the reported \code{scale}
is \eqn{\gamma} for the original supplied shape. Multiplying the shape by a
positive constant changes \eqn{\gamma} inversely, not the fitted variance.
Extremely ill-scaled shapes can still produce numerical errors.

Custom shape parameters are not estimated. Wald covariance and working-Gaussian
prediction bands condition on the fixed supplied shape. Quasi-AIC and quasi-BIC
count three fitted parameters (\eqn{V_{\max}}, \eqn{K_m}, \eqn{\gamma});
they do not account for coefficients fitted externally inside a custom function
or for variance-model selection uncertainty. Free parameters, shapes depending on
the fitted mean parameters, and data-adaptive variance-function estimation are
outside this interface. The feature does not implement conditional third/fourth
moment estimation or imply an efficiency gain over identically weighted
nonlinear least squares. Prediction bands retain the existing working-normal
approximation rather than estimating skewed predictive quantiles.
}
\examples{
one_curve <- subset(sdl_demo, enzyme == "1111")

fit <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "sqrt"
)

fit_power <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "power",
  power = 0.4
)

fit_auto <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "auto",
  power_selection = "quasi_aic"
)

coef(fit)
confint(fit)
set.seed(1)
confint(
  fit,
  method = "bootstrap",
  B = 99,
  bootstrap_ci = "studentized",
  wild_weights = "mammen"
)
head(predict(fit, newdata = seq(0, 80, length.out = 6), interval = "prediction"))

h <- variance_function(function(S) 1 + S / 20,
                       label = "Var = gamma (1 + S/20)")
fit_custom <- fit_mm(one_curve$s_uM, one_curve$v_uM_per_min, variance = h)
summary(fit_custom)
confint(fit_custom)
head(predict(fit_custom, newdata = c(0, 10, 40), interval = "prediction"))

\donttest{
plot(fit)
report_mm(fit, interval_type = "confidence")
}
}
