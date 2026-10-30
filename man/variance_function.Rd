\name{variance_function}
\alias{variance_function}
\title{Create Working Variance Functions}
\usage{
variance_function(model = c("constant", "log", "sqrt", "cuberoot", "power"),
  power = NULL, label = NULL)
}
\arguments{
  \item{model}{A built-in model name or a vectorized function of concentration
  returning one finite, strictly positive numeric variance-shape value per input.}
  \item{power}{Optional exponent used when \code{model = "power"}.}
  \item{label}{Optional non-empty descriptive label for a user-defined function.}
}
\value{
A function of substrate concentration \code{x}.
}
\description{
Construct a working variance function for variance-aware Michaelis-Menten
fitting.
}
\details{
A custom function represents a fixed shape \eqn{h(S)} in
\eqn{\operatorname{Var}(Y\mid S)=\gamma h(S)}. It must be deterministic,
pointwise, and vectorized: a value at a concentration must not depend on other
elements of the input vector. It must return a numeric vector of the same length
as its input, not a scalar recycled across observations. Values must be strictly
positive and finite at all fitting and prediction concentrations, including zero
if it is used. Invalid custom outputs cause an error; they are not clipped.
For example, use \code{function(S) rep(1, length(S))} for a constant shape.

All constants captured by a custom function must be fixed, not parameters to be
estimated by \code{fit_mm()}. Keep its enclosing environment unchanged between
fitting and prediction/refitting; use a local closure to capture external constants
when needed. \code{label} describes the function but does not change its values.
Custom functions are supported by \code{fit_mm()}, \code{screen_mm()}, and
\code{group_mm()}, not currently by \code{cluster_mm()}.
}
\examples{
h_log <- variance_function("log")
h_sqrt <- variance_function("sqrt")
h_log(c(0, 5, 10))
h_sqrt(c(0, 5, 10))
h_custom <- variance_function(function(S) 1 + S / 20,
                              label = "Var = gamma (1 + S/20)")
h_custom(c(0, 5, 10))
}
