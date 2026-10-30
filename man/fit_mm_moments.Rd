\name{fit_mm_moments}
\alias{fit_mm_moments}
\alias{confint.fit_mm_moments}
\alias{summary.fit_mm_moments}
\alias{print.fit_mm_moments}
\alias{print.summary.fit_mm_moments}
\alias{predict.fit_mm_moments}
\alias{plot.fit_mm_moments}
\title{Fit Michaelis-Menten Conditional Higher-Moment Estimating Equations}
\description{
Fit independent observations using the first two conditional residual moments,
weighted by a covariance involving working third and fourth moments.
This is an optional joint estimating-equation extension of the default
profile-score estimator, not a parametric gamma or Student-t likelihood fit.
}
\usage{
fit_mm_moments(x, y,
  variance = c("log", "sqrt", "cuberoot", "constant", "power"),
  power = NULL,
  moments = c("estimate", "normal", "gamma", "student_t", "custom"),
  shape = NULL, df = NULL, skewness = NULL, kurtosis = NULL,
  covariance = c("model", "sandwich"), interval_level = 0.95,
  km_bounds = NULL, allow_zero_substrate = TRUE, control = list())

\method{confint}{fit_mm_moments}(object, parm = c("Vmax", "Km"),
  level = object$interval_level, method = c("wald", "bootstrap"),
  B = 399L, bootstrap_type = "pairs",
  bootstrap_ci = c("studentized", "basic", "percentile"), ...)

\method{summary}{fit_mm_moments}(object, level = object$interval_level,
  method = c("wald", "bootstrap"), B = 399L, bootstrap_type = "pairs",
  bootstrap_ci = c("studentized", "basic", "percentile"), ...)

\method{print}{fit_mm_moments}(x,
  digits = max(3L, getOption("digits") - 2L), ...)
\method{print}{summary.fit_mm_moments}(x,
  digits = max(3L, getOption("digits") - 2L), ...)

\method{predict}{fit_mm_moments}(object, newdata = NULL, se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = object$interval_level, ...)

\method{plot}{fit_mm_moments}(x, interval = TRUE,
  interval_type = "confidence", ...)
}
\arguments{
  \item{x,y}{Concentrations and responses. For print and plot methods,
  \code{x} is the fitted object or summary.}
  \item{variance}{A fixed built-in variance-shape name or a vectorized,
  fixed positive function \eqn{h(S)}. See \code{\link{variance_function}}.
  Automatic variance selection is not supported by this estimator.}
  \item{power}{Fixed positive exponent when \code{variance = "power"}.}
  \item{moments}{Working standardized third and fourth moments:
  \code{"estimate"} estimates constant skewness and kurtosis from centered,
  standardized pilot residuals; \code{"normal"} uses 0 and 3;
  \code{"gamma"} uses \eqn{2/\sqrt{a}} and \eqn{3+6/a};
  \code{"student_t"} uses 0 and \eqn{3+6/(\nu-4)};
  \code{"custom"} uses the supplied values or functions.}
  \item{shape}{Fixed gamma shape \eqn{a>0}; defaults to 4 in gamma mode.
  It is not estimated. Only used in gamma mode.}
  \item{df}{Fixed Student-t degrees of freedom \eqn{\nu>4}; defaults to 10
  in Student-t mode. Only used in that mode.}
  \item{skewness,kurtosis}{In custom mode, finite numbers or deterministic,
  pointwise vectorized functions of concentration. Functions must return one
  value per input. Kurtosis is the standardized fourth moment, not excess
  kurtosis. The difference \eqn{k-1-s^2} must be safely positive.
  Keep captured constants unchanged after fitting.}
  \item{covariance}{Default \code{"model"} uses joint inverse working
  information and requires correct higher-moment specifications for its
  interpretation. Optional \code{"sandwich"} uses derivatives and empirical
  outer products of the full joint estimating equations.}
  \item{interval_level,level}{Confidence level strictly between zero and one.}
  \item{km_bounds}{Positive search bounds for \eqn{K_m}, as in \code{fit_mm()}.}
  \item{allow_zero_substrate}{Whether zero concentrations are allowed.}
  \item{control}{Named list: \code{maxit} (100 iterations per continuation
  attempt), \code{tol} (\eqn{10^{-7}} scaled-score tolerance), and
  \code{min_step} (1/64 minimum continuation step).}
  \item{object}{Fitted object.}
  \item{parm}{Mean-parameter names to include in the confidence intervals.}
  \item{method}{Wald or pairs-bootstrap confidence intervals.}
  \item{B}{Number of bootstrap resamples; at least 20.}
  \item{bootstrap_type}{Only \code{"pairs"} is supported. Observation-level
  wild multipliers do not generally preserve the required higher moments.}
  \item{bootstrap_ci}{Studentized, basic, or percentile bootstrap interval.}
  \item{digits}{Printed digits.}
  \item{newdata}{New concentrations for mean prediction.}
  \item{se.fit}{Return standard errors for mean predictions.}
  \item{interval}{For prediction, \code{"none"} or \code{"confidence"}.
  \code{"prediction"} raises an error because four moments do not specify
  predictive quantiles. For plotting, logical indicating whether to show bands.}
  \item{interval_type}{Only pointwise mean-confidence bands are supported.}
  \item{...}{Additional arguments to the underlying summary, bootstrap,
  or graphical methods.}
}
\details{
Write \eqn{r_i=Y_i-\mu(S_i;\beta)}, \eqn{v_i=\gamma h(S_i)},
and \eqn{u_i=(r_i,r_i^2-v_i)^\top}. For working skewness \eqn{s_i}
and raw kurtosis \eqn{k_i}, the conditional covariance of \eqn{u_i} is
\deqn{\Omega_i=\begin{pmatrix}v_i & s_i v_i^{3/2}\\
 s_i v_i^{3/2} & (k_i-1)v_i^2\end{pmatrix}.}
The derivative matrix is
\deqn{D_i=\begin{pmatrix}g_{i1} & g_{i2} & 0\\0 & 0 & h(S_i)\end{pmatrix},}
where \eqn{g_i} is the two-parameter MM mean gradient. The estimator solves
\eqn{\sum_i D_i^\top\Omega_i^{-1}u_i=0} jointly for
\eqn{\theta=(V_{\max},K_m,\gamma)^\top}. This implements the conditional
first/second-moment score of Kim and Ma (2012), Section 3, with working
third and fourth moments. No error density is fitted.

Computation uses log parameters and continuation from the profile-score
pilot, with damped Newton steps. Working skewness is interpolated from zero
and kurtosis from three. A failed solve, boundary solution, or ill-conditioned
covariance causes an error rather than a silent return of the pilot estimate.
Convergence failures have class \code{inferMM_convergence_error}; their
messages identify the fitting stage, failure reason, total iterations,
scaled score norm, tolerance, and suggested checks. Continuation failures
also report the completed fraction. No fitted object or confidence
intervals are returned. The condition fields can be inspected using
\code{tryCatch()}.
Failures during non-Gaussian higher-moment continuation also suggest manually
retrying with \code{moments = "normal"} while keeping the same variance function.
This changes the working moments and guarantees neither convergence nor
efficiency under the original non-Gaussian moments. No automatic fallback
is performed.
Only the pilot has the original scalar/closed-form reduction; the new fit
uses a joint numerical solve. Multiple roots and weak identification remain
possible. Internal variance normalization is undone in the reported scale.
The pilot uses the root validation described in \code{\link{fit_mm}}.
Its validation errors have class \code{inferMM_profile_error}, also inheriting
from \code{inferMM_convergence_error}, with profile-specific diagnostic fields.

The optional sandwich is \eqn{A^{-1}BA^{-\top}/n}, with \eqn{A} the negative
sample-average derivative of all three equations and \eqn{B} their empirical
second-moment matrix. Numerical derivatives are evaluated in log parameters
and the covariance is transformed back. The default model covariance is the inverse of
\eqn{\sum_i D_i^\top\Omega_i^{-1}D_i}. Both include the covariance with the
estimated scale; mean-parameter intervals use the leading two-by-two block
of the joint inverse, not an inverse that treats the scale as known.

With correct conditional mean and variance, suitable regularity and a
consistent root, working higher moments need not be correct for consistency.
Efficiency in the conditional mean/variance model additionally requires
appropriate higher-moment estimation. The constant standardized-moment
pilot model is not nonparametric estimation of arbitrary concentration-varying
skewness/kurtosis. Its estimation is first-order orthogonal under the model,
but can materially affect finite-sample performance. Neither the sandwich
nor this feature establishes efficiency or nominal coverage for every sample.
Unlike ordinary weighted mean equations, misspecification of \eqn{h(S)} can
bias this estimator because it also uses the squared-residual equation.

Normal working moments recover the profile mean/variance equations.
Constant symmetric Student-t moments also recover those point estimates;
heavy tails alone do not imply an efficiency gain. Gamma moments describe
centered, standardized gamma errors, not strictly positive raw gamma errors.
At least ten observations are required; estimating higher moments below
50 observations issues an exploratory-use warning, not a universal threshold
for reliable inference. The usual small-sample summary policy also applies.

Pairs bootstrap refits this estimator and re-estimates pilot moments in
\code{"estimate"} mode. Predictive bands are deliberately unavailable;
mean-confidence bands are pointwise. Quasi-AIC and quasi-BIC are not reported
because these roots do not maximize the working Gaussian likelihood.
Variance-selection uncertainty and externally estimated variance-shape
parameters are not covered. Screening, grouped selection, and the clustered
routine retain their existing profile estimators; this extension is for
independent single curves only.
}
\value{
An object inheriting from \code{"fit_mm_moments"} and \code{"fit_mm"}.
\code{coef()}, \code{vcov()}, and \code{confint()} concern \eqn{V_{\max},K_m}.
Additional fields include \code{scale}, \code{vcov_joint} (all three
parameters), \code{vcov_joint_model}, \code{vcov_joint_sandwich},
\code{standardized_moments}, \code{moment_spec}, \code{covariance_type},
and \code{convergence} (scaled score norm, iterations, continuation stages,
and derivative reciprocal condition number).
}
\references{
Kim, M. and Ma, Y. (2012). The efficiency of the second-order nonlinear
least squares estimator and its extension. \emph{Annals of the Institute
of Statistical Mathematics}, 64, 751--764.
\doi{10.1007/s10463-011-0332-y}.
}
\examples{
set.seed(19)
dat <- simulate_mm_data(n = 200, error = "skewed", skew_shape = 4)
h <- function(S) 1 + 9 * S / (20 + S)
fit <- fit_mm_moments(dat$x, dat$y, variance = h)
summary(fit)
confint(fit, level = 0.95)
fit$vcov_joint

# Fixed working moments, not a gamma likelihood fit.
fit_gamma <- fit_mm_moments(dat$x, dat$y, variance = h,
                            moments = "gamma", shape = 4)
coef(fit_gamma)
predict(fit_gamma, c(10, 20, 80), interval = "confidence")
\donttest{
set.seed(20)
confint(fit, method = "bootstrap", B = 99)
plot(fit)
}
}
