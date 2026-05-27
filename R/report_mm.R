#' Print and Plot a Michaelis-Menten Fit
#'
#' Convenience wrapper that prints a fitted-model summary and draws the fitted
#' curve with a confidence or prediction band.
#'
#' @param object A fitted object returned by `fit_mm()` or `cluster_mm()`.
#' @param level Confidence level used in the printed summary and plot.
#' @param method Interval method used in the printed summary: `"wald"` or
#'   `"bootstrap"`.
#' @param B Number of bootstrap resamples when `method = "bootstrap"`.
#' @param bootstrap_type Bootstrap scheme used when `method = "bootstrap"`:
#'   `"wild"` or `"pairs"`.
#' @param bootstrap_ci Bootstrap confidence-interval calibration used when
#'   `method = "bootstrap"`: `"studentized"`, `"basic"`, or `"percentile"`.
#' @param wild_weights Wild-bootstrap multipliers used when
#'   `bootstrap_type = "wild"`: `"mammen"` or `"rademacher"`.
#' @param seed Deprecated and ignored. Call [stats::set.seed()] before
#'   `report_mm()` when reproducible bootstrap output is needed.
#' @param interval_type Character string indicating whether the plotted band
#'   should be `"confidence"`, `"prediction"`, or `"none"`.
#' @param digits Number of printed digits.
#' @param ... Additional graphical arguments passed to `plot()`.
#'
#' @return The original fitted object, invisibly. When `method = "bootstrap"`
#'   and `object` inherits from `"fit_mm"`, the printed summary uses bootstrap
#'   standard errors and confidence intervals; the plotted band remains the
#'   usual pointwise confidence or prediction band from `plot.fit_mm()`.
#' @export
report_mm <- function(object,
                      level = object$interval_level,
                      method = c("wald", "bootstrap"),
                      B = 399L,
                      bootstrap_type = c("wild", "pairs"),
                      bootstrap_ci = c("studentized", "basic", "percentile"),
                      wild_weights = c("mammen", "rademacher"),
                      seed = NULL,
                      interval_type = c("confidence", "prediction", "none"),
                      digits = max(3L, getOption("digits") - 2L),
                      ...) {
  if (!inherits(object, "fit_mm") && !inherits(object, "cluster_mm")) {
    stop("object must be a fitted 'fit_mm' or 'cluster_mm' object.")
  }

  method <- match.arg(method)
  if (inherits(object, "cluster_mm") && missing(bootstrap_type)) {
    bootstrap_type <- "pairs"
  }

  bootstrap_type <- if (inherits(object, "cluster_mm")) {
    match.arg(bootstrap_type, c("pairs"))
  } else {
    match.arg(bootstrap_type)
  }
  bootstrap_ci <- match.arg(bootstrap_ci)
  wild_weights <- match.arg(wild_weights)
  interval_type <- match.arg(interval_type)
  summary_obj <- summary(
    object,
    level = level,
    method = method,
    B = B,
    bootstrap_type = bootstrap_type,
    bootstrap_ci = bootstrap_ci,
    wild_weights = wild_weights,
    seed = seed
  )
  print(summary_obj, digits = digits)

  if (!identical(interval_type, "none")) {
    plot(
      object,
      interval = TRUE,
      interval_type = interval_type,
      level = level,
      ...
    )
  } else {
    plot(
      object,
      interval = FALSE,
      level = level,
      ...
    )
  }

  invisible(object)
}
