weighted_score_function <- function(km, x, y, h_vals) {
  w <- 1 / h_vals
  left_1 <- sum(w * x * y / (km + x))
  right_1 <- sum(w * x^2 / (km + x)^3)
  left_2 <- sum(w * x * y / (km + x)^2)
  right_2 <- sum(w * x^2 / (km + x)^2)
  left_1 * right_1 - left_2 * right_2
}

weighted_vmax_hat <- function(km, x, y, h_vals) {
  w <- 1 / h_vals
  numerator <- sum(w * x * y / (km + x))
  denominator <- sum(w * x^2 / (km + x)^2)
  numerator / denominator
}

fit_mm_weighted <- function(x,
                            y,
                            variance_spec,
                            interval_level = 0.95,
                            km_bounds = default_km_bounds(x),
                            call = NULL,
                            allow_zero_substrate = TRUE,
                            variance_selection_mode = NULL,
                            variance_selection = NULL,
                            power_grid = NULL,
                            candidate_table = NULL) {
  h_raw <- evaluate_variance_function(variance_spec$fun, x)
  custom <- identical(variance_spec$model, "custom")
  # Rescaling a fixed custom shape changes gamma, not the fitted variance.
  h_normalization <- if (custom) exp(mean(log(h_raw))) else 1
  h_vals <- h_raw / h_normalization
  if (custom && (any(!is.finite(h_vals)) || any(h_vals <= 0) ||
                 any(!is.finite(1 / h_vals)))) {
    stop("The custom variance shape has too large a dynamic range for numerical fitting.",
         call. = FALSE)
  }
  profile <- find_km_root(x, y, h_vals, km_bounds)
  km_hat <- profile$root
  vmax_hat <- profile$vmax
  fitted <- vmax_hat * (x / (km_hat + x))
  residuals <- y - fitted
  normalized_scale <- mean(residuals^2 / h_vals)
  gamma_hat <- normalized_scale / h_normalization
  if (any(!is.finite(c(normalized_scale, gamma_hat))) ||
      normalized_scale <= 0 || gamma_hat <= 0) {
    profile_fit_error("variance-scale validation", "non_positive_parameters", km_bounds,
                      profile$convergence$score_norm, profile$convergence$tolerance,
                      profile$convergence$candidate_count)
  }
  grad <- mm_grad(x, vmax_hat, km_hat)
  covariance <- profile_covariance(grad, h_vals, normalized_scale, km_bounds,
                                    profile$convergence)
  vcov_beta <- covariance$covariance
  profile$convergence$information_rcond <- covariance$information_rcond
  beta <- c(Vmax = unname(vmax_hat), Km = unname(km_hat))
  variance_hat <- normalized_scale * h_vals
  if (custom && (!all(is.finite(c(beta, vcov_beta, gamma_hat))) || gamma_hat <= 0)) {
    stop("The custom variance fit produced non-finite estimates or a non-positive scale; check the data and variance shape.",
         call. = FALSE)
  }
  criteria <- quasi_metrics(y, fitted, variance_hat, k = 3L, clamp = !custom)

  structure(
    c(
      list(
        coefficients = beta,
        vcov = vcov_beta,
        scale = gamma_hat,
        normalized_scale = normalized_scale,
        variance_normalization = h_normalization,
        fitted.values = fitted,
        residuals = residuals,
        variance_hat = variance_hat,
        x = x,
        y = y,
        variance_model = variance_spec$model,
        variance_label = variance_spec$label,
        variance_fun = variance_spec$fun,
        variance_power = attr(variance_spec$fun, "variance_power", exact = TRUE),
        variance_selection_mode = variance_selection_mode,
        variance_selection = variance_selection,
        power_grid = power_grid,
        variance_candidate_table = candidate_table,
        method = "profile-score",
        interval_level = interval_level,
        km_bounds = km_bounds,
        convergence = profile$convergence,
        allow_zero_substrate = allow_zero_substrate,
        call = call
      ),
      criteria
    ),
    class = "fit_mm"
  )
}

fit_mm_from_spec <- function(x,
                             y,
                             variance,
                             power = NULL,
                             interval_level = 0.95,
                             km_bounds = NULL,
                             allow_zero_substrate = TRUE,
                             call = NULL,
                             variance_selection_mode = NULL,
                             variance_selection = NULL,
                             power_grid = NULL,
                             candidate_table = NULL) {
  variance_spec <- resolve_variance_spec(
    variance = variance,
    power = power
  )
  if (is.null(km_bounds)) {
    km_bounds <- default_km_bounds(x)
  }

  fit_mm_weighted(
    x = x,
    y = y,
    variance_spec = variance_spec,
    interval_level = interval_level,
    km_bounds = km_bounds,
    call = call,
    allow_zero_substrate = allow_zero_substrate,
    variance_selection_mode = variance_selection_mode,
    variance_selection = variance_selection,
    power_grid = power_grid,
    candidate_table = candidate_table
  )
}

search_auto_variance <- function(x,
                                 y,
                                 power_selection = c("quasi_aic", "quasi_bic"),
                                 power_grid = default_power_grid(),
                                 interval_level = 0.95,
                                 km_bounds = NULL,
                                 allow_zero_substrate = TRUE,
                                 call = NULL) {
  power_selection <- match.arg(power_selection)
  power_grid <- validate_power_grid(power_grid)

  candidates <- c(
    list(list(model = "log", power = NULL)),
    lapply(power_grid, function(power) list(model = "power", power = power))
  )

  fits <- lapply(candidates, function(spec) {
    tryCatch(fit_mm_from_spec(
      x = x,
      y = y,
      variance = spec$model,
      power = spec$power,
      interval_level = interval_level,
      km_bounds = km_bounds,
      allow_zero_substrate = allow_zero_substrate,
      call = call
    ), error = function(e) e)
  })

  candidate_table <- do.call(
    rbind,
    lapply(seq_along(fits), function(i) {
      fit <- fits[[i]]
      spec <- candidates[[i]]
      failed <- inherits(fit, "error")
      data.frame(
        candidate = if (identical(spec$model, "power")) {
          power_model_name(spec$power)
        } else {
          spec$model
        },
        label = resolve_variance_spec(spec$model, spec$power)$label,
        power = if (identical(spec$model, "power")) spec$power else NA_real_,
        quasi_aic = if (failed) NA_real_ else fit$quasi_aic,
        quasi_bic = if (failed) NA_real_ else fit$quasi_bic,
        rmse = if (failed) NA_real_ else fit$rmse,
        weighted_rss = if (failed) NA_real_ else fit$weighted_rss,
        success = as.integer(!failed),
        error_message = if (failed) conditionMessage(fit) else NA_character_,
        stringsAsFactors = FALSE
      )
    })
  )
  row.names(candidate_table) <- NULL

  successful <- which(candidate_table$success == 1L & is.finite(candidate_table[[power_selection]]))
  if (!length(successful)) {
    profile_fit_error("variance-model selection", "no_successful_candidate",
                      if (is.null(km_bounds)) default_km_bounds(x) else km_bounds,
                      detail = paste(candidate_table$error_message, collapse = " | "))
  }
  best_idx <- successful[which.min(candidate_table[[power_selection]][successful])]
  best_fit <- fits[[best_idx]]
  best_fit$variance_selection_mode <- "auto"
  best_fit$variance_selection <- power_selection
  best_fit$power_grid <- power_grid
  best_fit$variance_candidate_table <- candidate_table
  best_fit$variance_label <- auto_selected_label(best_fit$variance_label, power_selection)
  best_fit
}

#' Fit a Variance-Aware Michaelis-Menten Model
#'
#' Estimate Michaelis-Menten parameters under a constant or
#' concentration-dependent working variance specification.
#'
#' @param x Numeric vector of substrate concentrations.
#' @param y Numeric vector of observed reaction velocities.
#' @param variance A built-in model name or a fixed, vectorized, strictly positive
#'   function of substrate concentration. Use `variance_function()` to label a
#'   custom function. Its shape parameters are not estimated.
#'   Use `"power"` with a fixed exponent or `"auto"` to compare `log(S + 1)`
#'   against a grid of power functions.
#' @param power Optional exponent used when `variance = "power"`.
#' @param power_selection Criterion used when `variance = "auto"`.
#' @param power_grid Candidate exponents considered when `variance = "auto"`.
#' @param interval_level Confidence level used for Wald intervals.
#' @param km_bounds Numeric vector of length two giving the search bounds for
#'   `K_m`.
#' @param allow_zero_substrate Logical; if `FALSE`, substrate concentrations must
#'   be strictly positive.
#'
#' @return An object of class `"fit_mm"`.
#' @export
fit_mm <- function(x,
                   y,
                   variance = c("log", "sqrt", "cuberoot", "constant", "power", "auto"),
                   power = NULL,
                   power_selection = c("quasi_aic", "quasi_bic"),
                   power_grid = default_power_grid(),
                   interval_level = 0.95,
                   km_bounds = NULL,
                   allow_zero_substrate = TRUE) {
  validated <- validate_mm_inputs(x, y, allow_zero_substrate = allow_zero_substrate)
  if (!is.function(variance)) {
    variance <- normalize_variance_model(match.arg(variance))
  }
  power_selection <- match.arg(power_selection)

  if (is.null(km_bounds)) {
    km_bounds <- default_km_bounds(validated$x)
  }
  km_bounds <- as.numeric(km_bounds)

  if (length(km_bounds) != 2L || any(!is.finite(km_bounds)) || km_bounds[1] <= 0 ||
      km_bounds[2] <= km_bounds[1]) {
    stop("km_bounds must be a finite increasing vector of length two.")
  }

  if (identical(variance, "power")) {
    if (is.null(power) || !is.finite(power) || length(power) != 1L || power <= 0) {
      stop("A positive finite scalar 'power' must be supplied when variance = 'power'.")
    }
  }

  if (identical(variance, "auto")) {
    return(
      search_auto_variance(
        x = validated$x,
        y = validated$y,
        power_selection = power_selection,
        power_grid = power_grid,
        interval_level = interval_level,
        km_bounds = km_bounds,
        allow_zero_substrate = allow_zero_substrate,
        call = match.call()
      )
    )
  }

  fit_mm_from_spec(
    x = validated$x,
    y = validated$y,
    variance = variance,
    power = power,
    interval_level = interval_level,
    km_bounds = km_bounds,
    allow_zero_substrate = allow_zero_substrate,
    call = match.call(),
    variance_selection_mode = if (identical(variance, "power")) "manual" else NULL,
    variance_selection = if (identical(variance, "power")) "manual" else NULL,
    power_grid = if (identical(variance, "power")) power else NULL
  )
}
