default_variance_models <- function() {
  c("constant", "log", "sqrt", "cuberoot")
}

make_screen_row <- function(model,
                            label,
                            success,
                            runtime,
                            fit = NULL,
                            error_message = NA_character_) {
  if (!isTRUE(success)) {
    return(data.frame(
      model = model,
      model_label = label,
      selected_model = NA_character_,
      selected_model_label = NA_character_,
      success = 0,
      method = NA_character_,
      Vmax_hat = NA_real_,
      Km_hat = NA_real_,
      scale_hat = NA_real_,
      rmse = NA_real_,
      weighted_rss = NA_real_,
      quasi_aic = NA_real_,
      quasi_bic = NA_real_,
      runtime = runtime,
      error_message = error_message,
      stringsAsFactors = FALSE
    ))
  }

  selected_model <- selected_variance_model_name(
    fit$variance_model,
    fit$variance_power
  )
  selected_model_label <- selected_variance_label(
    fit$variance_model,
    fit$variance_power
  )

  data.frame(
    model = model,
    model_label = label,
    selected_model = selected_model,
    selected_model_label = selected_model_label,
    success = 1,
    method = fit$method,
    Vmax_hat = unname(coef(fit)["Vmax"]),
    Km_hat = unname(coef(fit)["Km"]),
    scale_hat = fit$scale,
    rmse = fit$rmse,
    weighted_rss = fit$weighted_rss,
    quasi_aic = fit$quasi_aic,
    quasi_bic = fit$quasi_bic,
    runtime = runtime,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Screen Working Variance Models
#'
#' Fit several working variance models to the same Michaelis-Menten curve and
#' order them primarily by quasi-AIC, while also reporting quasi-BIC, RMSE, and
#' weighted residual sums of squares for reference.
#'
#' @param x Numeric vector of substrate concentrations.
#' @param y Numeric vector of observed reaction velocities.
#' @param variance_models Character vector of built-in variance models to screen.
#'   The special value `"auto"` adds a candidate that compares `log(S + 1)` with
#'   a grid of power functions and keeps the better fit.
#' @param power_values Optional numeric vector of fixed exponents for additional
#'   `S^power` candidates.
#' @param include_auto Logical; if `TRUE`, add one automatic candidate that
#'   selects between `log(S + 1)` and the supplied `power_grid`.
#' @param power_selection Criterion used by the automatic candidate:
#'   `"quasi_aic"` or `"quasi_bic"`.
#' @param power_grid Candidate exponents used by the automatic candidate.
#' @param interval_level Confidence level forwarded to `fit_mm()`.
#' @param allow_zero_substrate Logical; if `FALSE`, substrate concentrations must
#'   be strictly positive.
#' @param km_bounds Optional search bounds for `K_m`.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return An object of class `"screen_mm"`.
#' @export
screen_mm <- function(x,
                           y,
                           variance_models = default_variance_models(),
                           power_values = NULL,
                           include_auto = FALSE,
                           power_selection = c("quasi_aic", "quasi_bic"),
                           power_grid = default_power_grid(),
                           interval_level = 0.95,
                           allow_zero_substrate = TRUE,
                           km_bounds = NULL,
                           quiet = FALSE) {
  validated <- validate_mm_inputs(x, y, allow_zero_substrate = allow_zero_substrate)
  power_selection <- match.arg(power_selection)
  candidates <- prepare_variance_candidates(
    variance_models = variance_models,
    power_values = power_values,
    include_auto = include_auto,
    power_selection = power_selection,
    power_grid = power_grid
  )
  fits <- list()
  rows <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    model <- candidate$id
    label <- candidate$label

    if (!isTRUE(quiet)) {
      message("Fitting model: ", model)
    }

    start_time <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      do.call(
        fit_mm,
        c(
          list(
            x = validated$x,
            y = validated$y,
            interval_level = interval_level,
            km_bounds = km_bounds,
            allow_zero_substrate = allow_zero_substrate
          ),
          candidate$fit_args
        )
      ),
      error = function(e) e
    )
    runtime <- proc.time()[["elapsed"]] - start_time

    if (inherits(fit, "error")) {
      rows[[i]] <- make_screen_row(
        model = model,
        label = label,
        success = FALSE,
        runtime = runtime,
        error_message = conditionMessage(fit)
      )
      next
    }

    fits[[model]] <- fit
    rows[[i]] <- make_screen_row(
      model = model,
      label = fit$variance_label,
      success = TRUE,
      runtime = runtime,
      fit = fit
    )
  }

  table <- do.call(rbind, rows)
  ok <- table$success == 1
  table$rank_quasi_aic <- NA_real_
  table$rank_quasi_bic <- NA_real_

  if (any(ok)) {
    table$rank_quasi_aic[ok] <- rank_metric(table$quasi_aic[ok])
    table$rank_quasi_bic[ok] <- rank_metric(table$quasi_bic[ok])
    table <- table[order(table$quasi_aic, table$quasi_bic, table$rmse), , drop = FALSE]
    row.names(table) <- NULL
  }

  structure(
    list(
      x = validated$x,
      y = validated$y,
      fits = fits,
      table = table,
      call = match.call()
    ),
    class = "screen_mm"
  )
}

#' @export
print.screen_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nModel screening summary:\n\n")

  show <- x$table
  preferred <- c(
    "model", "selected_model", "quasi_aic", "quasi_bic", "rmse",
    "success", "error_message"
  )
  keep <- c(preferred, setdiff(names(show), preferred))
  show <- show[, keep, drop = FALSE]
  numeric_cols <- vapply(show, is.numeric, logical(1))
  show[numeric_cols] <- lapply(show[numeric_cols], round, digits = digits)
  print(show, row.names = FALSE)

  invisible(x)
}
