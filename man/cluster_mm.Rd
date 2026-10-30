\name{cluster_mm}
\alias{cluster_mm}
\alias{coef.cluster_mm}
\alias{confint.cluster_mm}
\alias{plot.cluster_mm}
\alias{predict.cluster_mm}
\alias{print.cluster_mm}
\alias{summary.cluster_mm}
\alias{print.summary.cluster_mm}
\alias{vcov.cluster_mm}
\title{Fit a Cluster-Aware Michaelis-Menten Model}
\usage{
cluster_mm(data, s, v, cluster,
  variance = c("log", "sqrt", "cuberoot", "constant", "power", "auto"),
  power = NULL, power_selection = c("quasi_aic", "quasi_bic"),
  power_grid = default_power_grid(), interval_level = 0.95,
  invalid_response_values = c(-9999),
  invalid_substrate_values = numeric(0), allow_zero_substrate = TRUE,
  sort_x = TRUE, method = c("profile", "projected"), km_bounds = NULL)

\method{coef}{cluster_mm}(object, ...)

\method{vcov}{cluster_mm}(object, ...)

\method{confint}{cluster_mm}(object, parm = c("Vmax", "Km"),
  level = object$interval_level, method = c("wald", "bootstrap"),
  B = 399L, bootstrap_type = c("pairs"),
  bootstrap_ci = c("studentized", "basic", "percentile"),
  seed = NULL, ...)

\method{print}{cluster_mm}(x,
  digits = max(3L, getOption("digits") - 2L), ...)

\method{summary}{cluster_mm}(object, level = object$interval_level,
  method = c("wald", "bootstrap"), B = 399L,
  bootstrap_type = c("pairs"),
  bootstrap_ci = c("studentized", "basic", "percentile"),
  seed = NULL, ...)

\method{print}{summary.cluster_mm}(x,
  digits = max(3L, getOption("digits") - 2L), ...)

\method{predict}{cluster_mm}(object, newdata = NULL, se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = object$interval_level, ...)

\method{plot}{cluster_mm}(x, interval = TRUE,
  interval_type = c("prediction", "confidence"),
  level = x$interval_level, n_points = 200,
  xlab = "Substrate concentration", ylab = "Reaction velocity",
  main = NULL, pch = 19, col.points = "#1f78b4",
  col.line = "#222222",
  col.band = grDevices::adjustcolor("#9ecae1", alpha.f = 0.45), ...)
}
\arguments{
  \item{data}{A data frame containing clustered Michaelis-Menten observations.}
  \item{s}{Name of the substrate concentration column.}
  \item{v}{Name of the response column.}
  \item{cluster}{Character vector naming one or more clustering columns.
  Rows with missing values in any clustering column are removed with a warning
  before fitting, cluster counting, and variance-model screening.}
  \item{variance}{Character string naming a built-in residual working variance
  model. Use \code{"power"} with a fixed exponent or \code{"auto"} to compare
  \code{log(S + 1)} against a grid of power functions.}
  \item{power}{Optional exponent used when \code{variance = "power"}.}
  \item{power_selection}{Criterion used when \code{variance = "auto"}:
  \code{"quasi_aic"} or \code{"quasi_bic"}.}
  \item{power_grid}{Candidate exponents considered when
  \code{variance = "auto"}.}
  \item{interval_level}{Confidence level used for Wald intervals.}
  \item{invalid_response_values}{Values to remove before fitting.}
  \item{invalid_substrate_values}{Values to remove before fitting.}
  \item{allow_zero_substrate}{Logical; if \code{FALSE}, substrate
  concentrations must be strictly positive.}
  \item{sort_x}{Logical; sort observations within each cluster by substrate
  concentration.}
  \item{object}{A fitted \code{"cluster_mm"} object.}
  \item{parm}{Parameter names for \code{confint()}.}
  \item{level}{Confidence level for intervals.}
  \item{method}{In \code{cluster_mm()}, the estimator: \code{"profile"}
  retains the existing profile/moment algorithm; \code{"projected"} uses
  conditional projected scores with an identity instrument and a stacked
  cluster-sandwich covariance. In \code{summary()} and \code{confint()},
  the interval method: \code{"wald"} or a cluster-pairs \code{"bootstrap"}
  sensitivity analysis.}
  \item{km_bounds}{Optional positive increasing search bounds for
  \eqn{K_m}, available only with \code{method = "projected"}.}
  \item{B}{Number of bootstrap resamples when
  \code{method = "bootstrap"}.}
  \item{bootstrap_type}{Bootstrap scheme for clustered fits. The current
  implementation supports \code{"pairs"} only.}
  \item{bootstrap_ci}{Bootstrap confidence-interval calibration:
  \code{"studentized"}, \code{"basic"}, or \code{"percentile"}.}
  \item{seed}{Deprecated and ignored. Call \code{set.seed()} before
  bootstrap-based summaries or intervals if reproducibility is needed.}
  \item{x}{A fitted \code{"cluster_mm"} object for the print and plot methods.}
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
The main fitting function returns an object of class \code{"cluster_mm"}.
\code{predict()} returns either a numeric vector or a data frame containing
pointwise confidence or prediction intervals together with working and
marginal variance estimates.
}
\description{
Estimate shared Michaelis-Menten mean parameters from repeated or clustered
assay data using a random effect on \eqn{V_{\max}} and a residual working
variance function. An opt-in Ma-Genton-style conditional projected-score
construction leaves the random-effect density unspecified.
}
\details{
The current implementation treats repeated measurements through a working
covariance of the form \eqn{Z_i B Z_i^\top + \gamma A_i}, where \eqn{A_i}
is induced by the chosen residual variance function and the random-effect
component is specialized to a first-order perturbation of \eqn{V_{\max}}.
This provides a lightweight cluster-aware extension of the single-curve
profile-score estimator without requiring a fully parametric nonlinear mixed
effects specification.

Both methods remove observations with missing values (including numeric
\code{NaN}) in any clustering column before constructing clusters. A warning
reports the number of removed observations. Cluster counts, sizes, residuals,
and fitted values refer only to the retained data. Unused factor levels do not
count as clusters; a literal character ID \code{"NA"} is not a missing value.
If no observations or too few clusters remain, fitting stops with an error.

With \code{method = "profile"}, the routine initializes the fit from the corresponding
pooled single-curve estimator and then alternates a one-dimensional update for
\eqn{K_m} with plug-in updates for \eqn{V_{\max}} and moment-based updates for
the variance components. The pooled start is obtained from the same
single-curve profile-score machinery used by \code{fit_mm()}, including the
constant-variance special case. Fits that reach the maximum number of outer
iterations without meeting the convergence tolerance remain available but
carry an explicit convergence warning in the object summary.
The profile method permits singleton clusters alongside repeated-measurement
clusters; singleton covariance blocks retain their one-by-one matrix shape.
The projected method still requires at least two readings per retained cluster.

With \code{method = "projected"}, the assumed conditional model is
\eqn{Y_i = z_i(K_m) T_i + e_i}, where \eqn{T_i = V_{\max} + b_i},
\eqn{E(b_i\mid S_i)=0}, and residuals are conditionally independent Gaussian
with covariance \eqn{\gamma H_i}, \eqn{H_i=\mathrm{diag}(h(S_{ij}))}.
The random-effect density is not specified. Clusters must be independent,
the residual variance shape must be correct and fixed, and the random effect
must have finite fourth moments for the reported joint covariance. At least five clusters, at least two
readings per cluster, a positive loading in every cluster, and within-cluster
variation in positive concentrations are required.

For a candidate \eqn{k}, define \eqn{c_i=z_i^\top H_i^{-1}z_i},
\eqn{W_i=z_i^\top H_i^{-1}Y_i/c_i},
\eqn{u_i=Y_i-z_i W_i}, and \eqn{q_i=z_i'(k)^\top H_i^{-1}u_i}.
This is an MM-specific adaptation of the conditional-centering construction
of Ma and Genton (2010, Section 3.2), using the identity instrument
\eqn{a(W)=W}, not a general GLLVM implementation.
Conditional Gaussianity makes \eqn{W_i} and \eqn{u_i} independent given the
design, regardless of the random-effect density, so \eqn{E(W_iq_i\mid S_i)=0}.
The scalar equation \eqn{\sum_i W_iq_i=0} estimates \eqn{K_m}.
Then \eqn{\hat V_{\max}=\mathrm{mean}(W_i)} and
\eqn{\hat\gamma=\sum_i u_i^\top H_i^{-1}u_i/\sum_i(n_i-1)}.
The random-effect variance moment is
\eqn{\mathrm{mean}\{(W_i-\hat V_{\max})^2-\hat\gamma/c_i\}}.
A negative value is truncated to zero with a warning, only for reporting
the variance component and marginal variance; it does not enter the mean fit.

The new path searches 400 log-spaced points for scalar roots, compares
verified interior candidates by within-cluster weighted residual sum of
squares, and rejects a standardized score above \eqn{10^{-7}}. Invalid roots,
degenerate variance and singular derivative/covariance calculations stop the
fit, with class \code{inferMM_cluster_projected_error}. No silent Gaussian
fallback or covariance ridge is applied. The
\code{stage}, \code{reason} and \code{diagnostics} fields on failure conditions
to identify the numerical stage and available score/conditioning values.
An incomplete grid/refinement search can retain a verified interior root
but emits \code{inferMM_cluster_projected_warning}, with reason
\code{"incomplete_search"}; the message is stored in
\code{convergence$search_warning}. A non-positive variance moment instead
has warning reason \code{"variance_component_boundary"}. Neither warning
silently changes the estimation method. Both scalar and joint scaled
score residuals are checked against \eqn{10^{-7}}. Unrepresentable
original-unit scales or non-positive covariance diagonals stop fitting.
The joint covariance uses the
derivatives and empirical covariance of the four cluster-level equations for
\eqn{V_{\max}}, \eqn{K_m}, \eqn{\gamma}, and the unconstrained variance moment.
Its mean block drives standard errors and Wald intervals; it is a stacked
sandwich, not inverse efficient information. \code{vcov_joint} labels the
last coordinate \code{tau2_unconstrained}. General semiparametric efficiency
is not claimed for the identity instrument. Variance-model screening and
its intervals condition on the selected variance shape. Quasi-AIC/BIC are
retained only as working-normal marginal screening scores, not likelihood
information criteria for the unspecified latent-density model.

Projected fits provide mean-confidence intervals. Prediction intervals are
not available because an unspecified random-effect density does not determine
predictive quantiles. The default plot for a projected fit uses a confidence
band instead. This differs from the legacy working-normal prediction band.

Because small numbers of clusters can make first-order interval calibration
fragile, printed summaries suppress interval output by default when fewer than
four clusters are available or when there are fewer than six distinct
substrate concentrations. In that regime, \code{confint()} remains available
as an explicit sensitivity analysis and issues a warning.
}
\references{
Ma, Y. and Genton, M. G. (2010). Explicit estimating equations for
semiparametric generalized linear latent variable models.
\emph{Journal of the Royal Statistical Society: Series B}, 72, 475--495.
\doi{10.1111/j.1467-9868.2010.00741.x}.
}
\examples{
cluster_fit <- cluster_mm(
  data = subset(alves_demo, enzyme == "BG"),
  s = "substrate_conc",
  v = "activity",
  cluster = "core",
  variance = "sqrt"
)

coef(cluster_fit)
summary(cluster_fit)
confint(cluster_fit)
head(
  predict(
    cluster_fit,
    newdata = seq(0, 700, length.out = 6),
    interval = "confidence"
  )
)

set.seed(351)
S <- c(2, 5, 10, 20, 50, 100, 200)
b <- 12 * (rgamma(60, shape = 2) - 2) / sqrt(2)
dat <- data.frame(id = rep(seq_len(60), each = length(S)), S = rep(S, 60))
dat$Y <- (100 + b[dat$id]) * dat$S / (20 + dat$S) +
  rnorm(nrow(dat), sd = sqrt(2 * sqrt(dat$S)))
projected <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt",
                        method = "projected")
summary(projected)
confint(projected)
predict(projected, c(10, 20, 80), interval = "confidence")

\donttest{
plot(cluster_fit, interval_type = "confidence")
report_mm(cluster_fit, interval_type = "confidence")
}
}
