weighted_score_function <- function(km, x, y, h_vals) {
  w <- 1 / safe_pos(h_vals)
  left_1 <- sum(w * x * y / (km + x))
  right_1 <- sum(w * x^2 / (km + x)^3)
  left_2 <- sum(w * x * y / (km + x)^2)
  right_2 <- sum(w * x^2 / (km + x)^2)
  left_1 * right_1 - left_2 * right_2
}

weighted_vmax_hat <- function(km, x, y, h_vals) {
  w <- 1 / safe_pos(h_vals)
  numerator <- sum(w * x * y / (km + x))
  denominator <- sum(w * x^2 / (km + x)^2)
  numerator / denominator
}

find_km_root <- function(score_fun, bounds, grid_n = 400L) {
  grid <- exp(seq(log(bounds[1]), log(bounds[2]), length.out = grid_n))
  values <- vapply(
    grid,
    function(k) {
      val <- tryCatch(score_fun(k), error = function(e) NA_real_)
      if (!is.finite(val)) {
        return(NA_real_)
      }
      val
    },
    numeric(1)
  )

  ok_pair <- is.finite(values[-length(values)]) & is.finite(values[-1])
  sign_change <- values[-length(values)] * values[-1] <= 0
  idx <- which(ok_pair & sign_change)

  if (length(idx) > 0L) {
    hit <- idx[1]
    return(stats::uniroot(score_fun, interval = c(grid[hit], grid[hit + 1]))$root)
  }

  objective <- function(log_k) {
    val <- tryCatch(score_fun(exp(log_k)), error = function(e) Inf)
    if (!is.finite(val)) {
      return(Inf)
    }
    val^2
  }

  exp(stats::optimize(objective, interval = log(bounds))$minimum)
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
  h_raw <- variance_spec$fun(x)
  h_vals <- safe_pos(h_raw)
  km_hat <- find_km_root(
    score_fun = function(k) weighted_score_function(k, x, y, h_vals),
    bounds = km_bounds
  )
  vmax_hat <- weighted_vmax_hat(km_hat, x, y, h_vals)
  fitted <- mm_mean(x, vmax_hat, km_hat)
  residuals <- y - fitted
  gamma_hat <- mean(residuals^2 / h_vals)
  grad <- mm_grad(x, vmax_hat, km_hat)
  info <- crossprod(grad, grad / h_vals) / length(x)
  vcov_beta <- gamma_hat * safe_inverse(info) / length(x)
  beta <- c(Vmax = unname(vmax_hat), Km = unname(km_hat))
  variance_hat <- gamma_hat * h_vals
  criteria <- quasi_metrics(y, fitted, variance_hat, k = 3L)

  structure(
    c(
      list(
        coefficients = beta,
        vcov = vcov_beta,
        scale = gamma_hat,
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
        allow_zero_substrate = allow_zero_substrate,
        call = call
      ),
      criteria
    ),
    class = "fit_mm"
  )
}

fit_mm_nls <- function(x,
                       y,
                       interval_level = 0.95,
                       call = NULL,
                       allow_zero_substrate = TRUE,
                       variance_selection_mode = NULL,
                       variance_selection = NULL,
                       power_grid = NULL,
                       candidate_table = NULL) {
  dat <- data.frame(x = x, y = y)
  start_v <- max(y, na.rm = TRUE)
  if (!is.finite(start_v) || start_v <= 0) {
    start_v <- 1
  }

  fit <- stats::nls(
    y ~ (V * x) / (K + x),
    data = dat,
    start = list(V = 1.05 * start_v, K = stats::median(x)),
    control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
  )

  cf <- stats::coef(fit)
  beta <- c(Vmax = unname(cf[["V"]]), Km = unname(cf[["K"]]))
  fitted <- mm_mean(x, beta["Vmax"], beta["Km"])
  residuals <- y - fitted
  sigma2_hat <- mean(residuals^2)
  variance_hat <- rep(sigma2_hat, length(x))
  vcov_beta <- as.matrix(stats::vcov(fit))
  dimnames(vcov_beta) <- list(names(beta), names(beta))
  criteria <- quasi_metrics(y, fitted, variance_hat, k = 3L)

  structure(
    c(
      list(
        coefficients = beta,
        vcov = vcov_beta,
        scale = sigma2_hat,
        fitted.values = fitted,
        residuals = residuals,
        variance_hat = variance_hat,
        x = x,
        y = y,
        variance_model = "constant",
        variance_label = model_display_label("constant"),
        variance_fun = variance_function("constant"),
        variance_power = NULL,
        variance_selection_mode = variance_selection_mode,
        variance_selection = variance_selection,
        power_grid = power_grid,
        variance_candidate_table = candidate_table,
        method = "nls",
        interval_level = interval_level,
        km_bounds = default_km_bounds(x),
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

  if (identical(variance_spec$model, "constant")) {
    return(
      fit_mm_nls(
        x = x,
        y = y,
        interval_level = interval_level,
        call = call,
        allow_zero_substrate = allow_zero_substrate,
        variance_selection_mode = variance_selection_mode,
        variance_selection = variance_selection,
        power_grid = power_grid,
        candidate_table = candidate_table
      )
    )
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
    fit_mm_from_spec(
      x = x,
      y = y,
      variance = spec$model,
      power = spec$power,
      interval_level = interval_level,
      km_bounds = km_bounds,
      allow_zero_substrate = allow_zero_substrate,
      call = call
    )
  })

  candidate_table <- do.call(
    rbind,
    lapply(fits, function(fit) {
      data.frame(
        candidate = if (identical(fit$variance_model, "power")) {
          power_model_name(fit$variance_power)
        } else {
          fit$variance_model
        },
        label = fit$variance_label,
        power = if (identical(fit$variance_model, "power")) fit$variance_power else NA_real_,
        quasi_aic = fit$quasi_aic,
        quasi_bic = fit$quasi_bic,
        rmse = fit$rmse,
        weighted_rss = fit$weighted_rss,
        stringsAsFactors = FALSE
      )
    })
  )
  row.names(candidate_table) <- NULL

  best_idx <- which.min(candidate_table[[power_selection]])[1]
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
#' @param variance Character string naming a built-in working variance model.
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
  variance <- normalize_variance_model(match.arg(variance))
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
