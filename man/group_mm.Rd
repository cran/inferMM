\name{group_mm}
\alias{group_mm}
\alias{plot.group_mm}
\alias{print.group_mm}
\title{Fit Multiple Michaelis-Menten Curves by Group}
\usage{
group_mm(data, s, v, groups = NULL,
  variance_models = default_variance_models(), power_values = NULL,
  include_auto = FALSE, power_selection = c("quasi_aic", "quasi_bic"),
  power_grid = default_power_grid(), interval_level = 0.95,
  invalid_response_values = c(-9999),
  invalid_substrate_values = numeric(0), allow_zero_substrate = TRUE,
  sort_x = TRUE, quiet = FALSE)

\method{print}{group_mm}(x,
  digits = max(3L, getOption("digits") - 2L), ...)

\method{plot}{group_mm}(x, which_groups = NULL, model = NULL,
  interval = TRUE, interval_type = c("prediction", "confidence"),
  level = NULL, ncol = NULL, main = NULL,
  xlab = "Substrate concentration", ylab = "Reaction velocity", ...)
}
\arguments{
  \item{data}{A data frame containing Michaelis-Menten observations.}
  \item{s}{Name of the substrate concentration column.}
  \item{v}{Name of the response column.}
  \item{groups}{Optional character vector of grouping columns.}
  \item{variance_models}{Character vector of built-in variance models. The
  special value \code{"auto"} adds a candidate that compares \code{log(S + 1)}
  with a grid of power functions and keeps the better fit.}
  \item{power_values}{Optional numeric vector of fixed exponents for
  additional \code{S^power} candidates.}
  \item{include_auto}{Logical; if \code{TRUE}, add one automatic candidate that
  selects between \code{log(S + 1)} and the supplied \code{power_grid}.}
  \item{power_selection}{Criterion used by the automatic candidate:
  \code{"quasi_aic"} or \code{"quasi_bic"}.}
  \item{power_grid}{Candidate exponents used by the automatic candidate.}
  \item{interval_level}{Confidence level forwarded to \code{\link{fit_mm}}.}
  \item{invalid_response_values}{Values to remove before fitting.}
  \item{invalid_substrate_values}{Values to remove before fitting.}
  \item{allow_zero_substrate}{Logical; if \code{FALSE}, substrate concentrations
  must be strictly positive.}
  \item{sort_x}{Logical; sort observations within each group by concentration.}
  \item{quiet}{Logical; if \code{TRUE}, suppress progress messages.}
  \item{x}{A grouped fit object for the print and plot methods.}
  \item{which_groups}{Optional subset of groups to plot, matched against
  \code{group_id} or \code{group_label}.}
  \item{model}{Optional model name to use for every group. When \code{NULL},
  the plotting method uses the best model recorded for each group.}
  \item{interval}{Logical; if \code{TRUE}, add interval bands to each panel.}
  \item{interval_type}{Either \code{"prediction"} or \code{"confidence"}.}
  \item{level}{Optional interval level. Defaults to the level stored in each fit.}
  \item{ncol}{Optional number of panel columns in the grouped display.}
  \item{main}{Optional panel title or titles.}
  \item{xlab}{Label for the x-axis.}
  \item{ylab}{Label for the y-axis.}
  \item{digits}{Number of printed digits.}
  \item{...}{Additional arguments passed to S3 methods.}
}
\value{
An object of class \code{"group_mm"} containing cleaned data, fit
summaries, augmented residual output, grouped rankings, and fitted objects.
}
\description{
Clean a data frame, fit one or more working variance models within each group,
and compare the competing fits using quasi-AIC as the primary criterion, with
quasi-BIC and RMSE reported alongside it.
}
\examples{
grouped <- group_mm(
  data = sdl_demo,
  s = "s_uM",
  v = "v_uM_per_min",
  groups = "enzyme",
  variance_models = c("constant", "sqrt"),
  power_values = 0.4,
  include_auto = TRUE,
  quiet = TRUE
)

grouped
head(
  grouped$comparison$best_by_group[
    ,
    c("group_label", "model", "selected_model", "quasi_aic", "quasi_bic", "rmse")
  ]
)

plot(grouped, interval_type = "confidence")
}
