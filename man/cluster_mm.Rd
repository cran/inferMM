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
  sort_x = TRUE)

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
  \item{cluster}{Character vector naming one or more clustering columns.}
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
  \item{method}{Interval construction method. The current implementation
  supports \code{"wald"} or a cluster-pairs \code{"bootstrap"}
  sensitivity analysis.}
  \item{B}{Number of bootstrap resamples when
  \code{method = "bootstrap"}.}
  \item{bootstrap_type}{Bootstrap scheme for clustered fits. The current
  implementation supports \code{"pairs"} only.}
  \item{bootstrap_ci}{Bootstrap confidence-interval calibration:
  \code{"studentized"}, \code{"basic"}, or \code{"percentile"}.}
  \item{seed}{Deprecated and ignored. Call \code{stats::set.seed()} before
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
assay data using a Ma-Genton-style working covariance with a random effect on
\eqn{V_{\max}} and a low-dimensional residual working variance function.
}
\details{
The current implementation treats repeated measurements through a working
covariance of the form \eqn{Z_i B Z_i^\top + \gamma A_i}, where \eqn{A_i}
is induced by the chosen residual variance function and the random-effect
component is specialized to a first-order perturbation of \eqn{V_{\max}}.
This provides a lightweight cluster-aware extension of the single-curve
profile-score estimator without requiring a fully parametric nonlinear mixed
effects specification.

Numerically, the routine initializes the clustered fit from the corresponding
pooled single-curve estimator and then alternates a one-dimensional update for
\eqn{K_m} with plug-in updates for \eqn{V_{\max}} and moment-based updates for
the variance components. The pooled start is obtained from the same
single-curve profile-score machinery used by \code{fit_mm()}, including the
constant-variance special case. Fits that reach the maximum number of outer
iterations without meeting the convergence tolerance remain available but
carry an explicit convergence warning in the object summary.

Because small numbers of clusters can make first-order interval calibration
fragile, printed summaries suppress interval output by default when fewer than
four clusters are available or when there are fewer than six distinct
substrate concentrations. In that regime, \code{confint()} remains available
as an explicit sensitivity analysis and issues a warning.
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

\donttest{
plot(cluster_fit, interval_type = "confidence")
report_mm(cluster_fit, interval_type = "confidence")
}
}
