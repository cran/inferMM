\name{screen_mm}
\alias{screen_mm}
\alias{print.screen_mm}
\title{Screen Working Variance Models}
\usage{
screen_mm(x, y, variance_models = default_variance_models(),
  power_values = NULL, include_auto = FALSE,
  power_selection = c("quasi_aic", "quasi_bic"),
  power_grid = default_power_grid(), interval_level = 0.95,
  allow_zero_substrate = TRUE, km_bounds = NULL, quiet = FALSE,
  variance_functions = NULL)

\method{print}{screen_mm}(x, digits = max(3L, getOption("digits") - 2L), ...)
}
\arguments{
  \item{x}{Numeric vector of substrate concentrations, or a screening object for
  the print method.}
  \item{y}{Numeric vector of observed reaction velocities.}
  \item{variance_models}{Character vector of built-in variance models to
  compare. The special value \code{"auto"} adds a candidate that compares
  \code{log(S + 1)} with a grid of power functions and keeps the better fit.}
  \item{power_values}{Optional numeric vector of fixed exponents for
  additional \code{S^power} candidates.}
  \item{include_auto}{Logical; if \code{TRUE}, add one automatic candidate that
  selects between \code{log(S + 1)} and the supplied \code{power_grid}.}
  \item{power_selection}{Criterion used by the automatic candidate:
  \code{"quasi_aic"} or \code{"quasi_bic"}.}
  \item{power_grid}{Candidate exponents used by the automatic candidate.}
  \item{interval_level}{Confidence level forwarded to \code{\link{fit_mm}}.}
  \item{allow_zero_substrate}{Logical; if \code{FALSE}, substrate concentrations
  must be strictly positive.}
  \item{km_bounds}{Optional search bounds for \eqn{K_m}.}
  \item{quiet}{Logical; if \code{TRUE}, suppress progress messages.}
  \item{variance_functions}{Optional named list of fixed, vectorized positive
  functions \eqn{h(S)}. Names must be unique, non-empty, and distinct from
  built-in candidate names. Set \code{variance_models = character(0)} to compare
  only these custom candidates.}
  \item{digits}{Number of printed digits.}
  \item{...}{Additional arguments passed to S3 methods.}
}
\value{
An object of class \code{"screen_mm"} with a fitted-model table and the
successful fits.
}
\description{
Fit several working variance models to the same Michaelis-Menten curve and
order them primarily by quasi-AIC, while also reporting quasi-BIC, RMSE, and
weighted residual sums of squares for reference.
}
\details{
Custom functions follow the contract in \code{\link{variance_function}} and
are fitted by \code{\link{fit_mm}}. Their names identify rows in the ranking
table and entries in \code{fits}. A candidate whose function fails validation
is recorded with \code{success = 0} and an \code{error_message}; other candidates
are still evaluated. Malformed or duplicate candidate names cause an error.
Custom functions are added to the outer comparison; they do not change the
log/power search performed by \code{include_auto}.

All custom shapes are treated as fixed and add no estimated shape parameters.
Criteria and intervals do not account for externally estimated function
coefficients or for selection uncertainty. Adding more candidates does not
guarantee better estimation or interval coverage.
}
\examples{
one_curve <- subset(sdl_demo, enzyme == "1111")
screen <- screen_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  power_values = c(0.4, 0.6),
  include_auto = TRUE,
  quiet = TRUE
)

screen$table[, c("model", "selected_model", "quasi_aic", "quasi_bic", "rmse")]

custom_screen <- screen_mm(
  one_curve$s_uM, one_curve$v_uM_per_min,
  variance_models = c("constant", "sqrt"),
  variance_functions = list(
    linear_floor = function(S) 1 + S / 20,
    saturation = function(S) 1 + 9 * S / (20 + S)
  ), quiet = TRUE
)
custom_screen$table[, c("model", "quasi_aic", "success")]
}
