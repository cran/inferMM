\name{report_mm}
\alias{report_mm}
\title{Print and Plot a Michaelis-Menten Fit}
\usage{
report_mm(object, level = object$interval_level,
  method = c("wald", "bootstrap"), B = 399L,
  bootstrap_type = c("wild", "pairs"),
  bootstrap_ci = c("studentized", "basic", "percentile"),
  wild_weights = c("mammen", "rademacher"), seed = NULL,
  interval_type = c("confidence", "prediction", "none"),
  digits = max(3L, getOption("digits") - 2L), ...)
}
\arguments{
  \item{object}{A fitted object returned by \code{\link{fit_mm}} or
  \code{\link{cluster_mm}}.}
  \item{level}{Confidence level used in the printed summary and plot.}
  \item{method}{Interval method used in the printed summary:
  \code{"wald"} or \code{"bootstrap"}.}
  \item{B}{Number of bootstrap resamples when
  \code{method = "bootstrap"}.}
  \item{bootstrap_type}{Bootstrap scheme used when
  \code{method = "bootstrap"}: \code{"wild"} or \code{"pairs"}.}
  \item{bootstrap_ci}{Bootstrap interval calibration used when
  \code{method = "bootstrap"}: \code{"studentized"}, \code{"basic"}, or
  \code{"percentile"}.}
  \item{wild_weights}{Wild-bootstrap multipliers used when
  \code{bootstrap_type = "wild"}: \code{"mammen"} or
  \code{"rademacher"}.}
  \item{seed}{Deprecated and ignored. Call \code{stats::set.seed()} before
  bootstrap-based summaries if reproducibility is needed.}
  \item{interval_type}{Character string indicating whether the plotted band
  should be \code{"confidence"}, \code{"prediction"}, or \code{"none"}.}
  \item{digits}{Number of printed digits.}
  \item{...}{Additional graphical arguments passed to \code{\link{plot}}.}
}
\value{
The original fitted object, returned invisibly. When
\code{method = "bootstrap"} and the object inherits from \code{"fit_mm"},
the printed summary uses bootstrap standard errors and confidence intervals;
the plotted band remains the usual pointwise confidence or prediction band
from \code{plot.fit_mm()}. For very small curves (fewer than 20 observations
or fewer than 6 distinct substrate concentrations), the printed summary
reports point estimates only by default and leaves interval construction to an
explicit \code{confint()} call. For sparse clustered fits, the printed summary
similarly suppresses interval output by default and treats \code{confint()}
as an explicit sensitivity analysis; bootstrap summaries are not currently
available for \code{"cluster_mm"} objects.
}
\description{
Convenience wrapper that prints a fitted-model summary and draws the fitted
curve with a confidence or prediction band.
}
\examples{
one_curve <- subset(sdl_demo, enzyme == "1111")
fit <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "sqrt"
)

\donttest{
report_mm(fit, interval_type = "confidence")
set.seed(1)
report_mm(
  fit,
  method = "bootstrap",
  B = 99,
  bootstrap_ci = "studentized",
  wild_weights = "mammen",
  interval_type = "confidence"
)

cluster_fit <- cluster_mm(
  data = subset(alves_demo, enzyme == "BG"),
  s = "substrate_conc",
  v = "activity",
  cluster = "core",
  variance = "sqrt"
)
report_mm(cluster_fit, interval_type = "confidence")
}
}
