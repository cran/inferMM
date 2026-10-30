moment_scalar <- function(x, name, lower = 0, strict = TRUE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || if (strict) x <= lower else x < lower) {
    stop(name, " must be a finite number ", if (strict) "> " else ">= ",
         lower, ".", call. = FALSE)
  }
  as.numeric(x)
}

moment_level <- function(level) {
  level <- moment_scalar(level, "confidence level")
  if (level >= 1) stop("confidence level must be less than 1.", call. = FALSE)
  level
}

moment_values <- function(value, x, name) {
  if (is.function(value)) {
    out <- tryCatch(value(x), error = function(e) {
      stop(name, " function failed: ", conditionMessage(e), call. = FALSE)
    })
  } else {
    if (!is.numeric(value) || length(value) != 1L) {
      stop(name, " must be a number or a vectorized function of concentration.",
           call. = FALSE)
    }
    out <- rep(value, length(x))
  }
  if (!is.numeric(out) || is.complex(out) || !is.null(dim(out)) ||
      length(out) != length(x) || any(!is.finite(out))) {
    stop(name, " must supply one finite numeric value per concentration.",
         call. = FALSE)
  }
  as.numeric(out)
}

prepare_error_moments <- function(pilot, moments, shape, df, skewness, kurtosis) {
  if (!is.null(shape) && moments != "gamma") {
    stop("shape is only used with moments = 'gamma'.", call. = FALSE)
  }
  if (!is.null(df) && moments != "student_t") {
    stop("df is only used with moments = 'student_t'.", call. = FALSE)
  }
  if (moments != "custom" && (!is.null(skewness) || !is.null(kurtosis))) {
    stop("skewness and kurtosis are only used with moments = 'custom'.",
         call. = FALSE)
  }
  if (moments == "estimate") {
    z <- pilot$residuals / sqrt(pilot$variance_hat)
    z <- z - mean(z)
    z <- z / sqrt(mean(z^2))
    s <- mean(z^3)
    k <- mean(z^4)
  } else if (moments == "normal") {
    s <- 0
    k <- 3
  } else if (moments == "gamma") {
    if (is.null(shape)) shape <- 4
    shape <- moment_scalar(shape, "shape")
    s <- 2 / sqrt(shape)
    k <- 3 + 6 / shape
  } else if (moments == "student_t") {
    if (is.null(df)) df <- 10
    df <- moment_scalar(df, "df", lower = 4)
    s <- 0
    k <- 3 + 6 / (df - 4)
  } else {
    if (is.null(skewness) || is.null(kurtosis)) {
      stop("Custom moments require both skewness and kurtosis (not excess kurtosis).",
           call. = FALSE)
    }
    s <- skewness
    k <- kurtosis
  }
  s <- moment_values(s, pilot$x, "skewness")
  k <- moment_values(k, pilot$x, "kurtosis")
  determinant <- k - 1 - s^2
  if (any(!is.finite(determinant)) || any(determinant <= 1e-6 * pmax(1, k))) {
    stop("Moment covariance must be positive definite: kurtosis - 1 - skewness^2 must be safely positive.",
         call. = FALSE)
  }
  list(skewness = s, kurtosis = k,
       spec = list(moments = moments, shape = shape, df = df,
                   skewness = skewness, kurtosis = kurtosis))
}

moment_control <- function(control) {
  defaults <- list(maxit = 100L, tol = 1e-7, min_step = 1 / 64)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults))))) {
    stop("control must be a named list containing maxit, tol, or min_step.",
         call. = FALSE)
  }
  defaults[names(control)] <- control
  maxit <- moment_scalar(defaults$maxit, "control$maxit")
  if (maxit != floor(maxit) || maxit > 1000) {
    stop("control$maxit must be an integer between 1 and 1000.", call. = FALSE)
  }
  defaults$maxit <- as.integer(maxit)
  defaults$tol <- moment_scalar(defaults$tol, "control$tol", lower = 1e-12)
  defaults$min_step <- moment_scalar(defaults$min_step, "control$min_step")
  if (defaults$tol > 1e-3 || defaults$min_step > 0.25) {
    stop("control$tol must not exceed 1e-3 and min_step must not exceed 0.25.",
         call. = FALSE)
  }
  defaults
}

# Optimal conditional-moment score D' Omega^{-1} (r, r^2 - v), in log parameters.
moment_score_state <- function(eta, x, y, h, reference, skewness, kurtosis) {
  theta <- reference * exp(eta)
  mu <- theta[1] * x / (theta[2] + x)
  variance <- theta[3] * h
  z <- (y - mu) / sqrt(variance)
  d1 <- mu / sqrt(variance)
  d2 <- -d1 * theta[2] / (theta[2] + x)
  det <- kurtosis - 1 - skewness^2
  q1 <- ((kurtosis - 1) * z - skewness * (z^2 - 1)) / det
  q2 <- (z^2 - 1 - skewness * z) / det
  scores <- cbind(Vmax = d1 * q1, Km = d2 * q1, gamma = q2)
  w <- (kurtosis - 1) / det
  information <- rbind(
    c(mean(w * d1^2), mean(w * d1 * d2), mean(-skewness * d1 / det)),
    c(mean(w * d1 * d2), mean(w * d2^2), mean(-skewness * d2 / det)),
    c(mean(-skewness * d1 / det), mean(-skewness * d2 / det), mean(1 / det))
  )
  if (any(!is.finite(c(theta, variance, scores, information)))) {
    stop("Non-finite higher-moment score; check data and moment scaling.", call. = FALSE)
  }
  list(theta = theta, mean = mu, variance = variance, scores = scores,
       equation = colMeans(scores), information = information)
}

moment_jacobian <- function(fun, eta) {
  step <- .Machine$double.eps^(1 / 3) * pmax(1, abs(eta))
  vapply(seq_along(eta), function(j) {
    plus <- minus <- eta
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (fun(plus) - fun(minus)) / (2 * step[j])
  }, numeric(length(eta)))
}

moment_newton <- function(fun, eta, lower, upper, control) {
  reason <- "iteration_limit"
  detail <- NULL
  score_norm <- Inf
  for (iteration in 0:control$maxit) {
    f <- tryCatch(fun(eta), error = function(e) e)
    if (inherits(f, "error") || any(!is.finite(f))) {
      reason <- "non_finite_score"
      score_norm <- Inf
      if (inherits(f, "error")) detail <- conditionMessage(f)
      break
    }
    score_norm <- max(abs(f))
    if (score_norm <= control$tol) {
      return(list(ok = TRUE, eta = eta, iterations = iteration,
                  score_norm = score_norm, reason = "converged"))
    }
    if (iteration == control$maxit) break
    jac <- tryCatch(moment_jacobian(fun, eta), error = function(e) e)
    if (inherits(jac, "error") || any(!is.finite(jac))) {
      reason <- "non_finite_jacobian"
      if (inherits(jac, "error")) detail <- conditionMessage(jac)
      break
    }
    merit <- sum(f^2)
    accepted <- FALSE
    direction_found <- FALSE
    for (damping in c(0, 1e-4, 0.01, 1, 100)) {
      direction <- tryCatch({
        if (damping == 0) {
          as.numeric(solve(jac, -f))
        } else {
          as.numeric(solve(crossprod(jac) + diag(damping, 3), -crossprod(jac, f)))
        }
      }, error = function(e) NULL)
      if (is.null(direction) || any(!is.finite(direction))) next
      direction_found <- TRUE
      direction <- direction / max(1, max(abs(direction)))
      for (j in 0:24) {
        fraction <- 2^(-j)
        candidate <- pmin(pmax(eta + fraction * direction, lower + 1e-9), upper - 1e-9)
        trial <- tryCatch(fun(candidate), error = function(e) rep(Inf, 3))
        if (all(is.finite(trial)) && sum(trial^2) < merit * (1 - 1e-4 * fraction)) {
          eta <- candidate
          accepted <- TRUE
          break
        }
      }
      if (accepted) break
    }
    if (!accepted) {
      reason <- if (direction_found) "no_decreasing_step" else "linear_solve_failed"
      break
    }
  }
  list(ok = FALSE, eta = eta, iterations = iteration,
       score_norm = score_norm, reason = reason, detail = detail)
}

moment_fit_error <- function(stage, attempt, moments, control,
                             iterations = attempt$iterations, progress = NULL) {
  reasons <- c(
    iteration_limit = "the iteration limit was reached",
    no_decreasing_step = "no damped Newton step reduced the score norm",
    linear_solve_failed = "the Newton linear systems could not be solved",
    non_finite_score = "the score evaluation was non-finite or failed",
    non_finite_jacobian = "the score derivative was non-finite or failed",
    search_boundary = "the solution is on a search boundary"
  )
  reason <- unname(reasons[attempt$reason])
  if (is.na(reason)) reason <- "the score equations did not converge"
  diagnostic <- sprintf(
    "Total iterations: %d; max absolute scaled score: %.3g; tolerance: %.3g.",
    iterations, attempt$score_norm, control$tol
  )
  if (!is.null(progress)) {
    diagnostic <- paste(diagnostic, sprintf("Completed continuation fraction: %.4g of 1.", progress))
  }
  advice <- if (identical(attempt$reason, "iteration_limit")) {
    "Consider increasing control$maxit and check substrate coverage, outliers, and the variance/moment specifications."
  } else if (identical(attempt$reason, "search_boundary")) {
    "Review km_bounds, substrate coverage, and parameter identifiability."
  } else {
    "Check substrate coverage, outliers, and the variance/moment specifications."
  }
  if (identical(stage, "higher-moment continuation") && !identical(moments, "normal")) {
    advice <- paste(
      advice,
      "Consider retrying with moments='normal' while keeping the same variance function.",
      "This changes the working moments; convergence and efficiency under the original non-Gaussian moments are not guaranteed.",
      "No automatic fallback is performed."
    )
  }
  message <- paste(
    sprintf("fit_mm_moments() did not converge to an acceptable root during %s (moments='%s'): %s.",
            stage, moments, reason),
    diagnostic, "No fit or confidence intervals are returned.", advice
  )
  if (!is.null(attempt$detail)) message <- paste(message, "Evaluation detail:", attempt$detail)
  stop(structure(
    list(message = message, call = NULL, stage = stage, reason = attempt$reason,
         moments = moments, iterations = iterations, score_norm = attempt$score_norm,
         tolerance = control$tol, progress = progress),
    class = c("inferMM_convergence_error", "error", "condition")
  ))
}

moment_inverse <- function(mat, name) {
  if (any(!is.finite(mat)) || !is.finite(rcond(mat)) || rcond(mat) < 1e-12) {
    stop(name, " is singular or poorly conditioned; intervals are not reported.",
         call. = FALSE)
  }
  solve(mat)
}

moment_profile_pilot <- function(x, y, variance, power, interval_level,
                                  km_bounds, allow_zero_substrate) {
  data <- validate_mm_inputs(x, y, allow_zero_substrate = allow_zero_substrate)
  x <- data$x
  y <- data$y
  spec <- resolve_variance_spec(variance, power)
  h_raw <- evaluate_variance_function(spec$fun, x)
  normalization <- exp(mean(log(h_raw)))
  h <- h_raw / normalization
  if (any(!is.finite(1 / h)) || any(!is.finite(h)) || any(h <= 0)) {
    stop("The variance shape has too large a dynamic range for higher-moment fitting.",
         call. = FALSE)
  }
  if (is.null(km_bounds)) km_bounds <- default_km_bounds(x)
  if (!is.numeric(km_bounds) || length(km_bounds) != 2L ||
      any(!is.finite(km_bounds)) || km_bounds[1] <= 0 || km_bounds[2] <= km_bounds[1]) {
    stop("km_bounds must be a finite increasing vector of length two.", call. = FALSE)
  }
  x_scale <- max(x)
  y_scale <- max(abs(y))
  if (!is.finite(y_scale) || y_scale <= 0) {
    stop("Higher-moment fitting requires a nondegenerate response.", call. = FALSE)
  }
  xs <- x / x_scale
  ys <- y / y_scale
  # Only the pilot point estimates are needed; do not invert its ill-scaled raw-unit information.
  profile <- find_km_root(x, y, h, km_bounds)
  km <- profile$root / x_scale
  vmax <- profile$vmax / y_scale
  mu <- mm_mean(xs, vmax, km)
  scale <- mean((ys - mu)^2 / h) * y_scale^2
  beta <- c(Vmax = vmax * y_scale, Km = km * x_scale)
  list(coefficients = beta, normalized_scale = scale, scale = scale / normalization,
       variance_normalization = normalization, variance_hat = scale * h,
       fitted.values = mu * y_scale, residuals = (ys - mu) * y_scale,
       x = x, y = y, variance_model = spec$model, variance_label = spec$label,
       variance_fun = spec$fun, variance_power = attr(spec$fun, "variance_power", exact = TRUE),
       variance_selection_mode = NULL, interval_level = interval_level,
       km_bounds = km_bounds, allow_zero_substrate = allow_zero_substrate)
}

#' Fit Michaelis-Menten Conditional Higher-Moment Estimating Equations
#' @export
fit_mm_moments <- function(x, y,
                           variance = c("log", "sqrt", "cuberoot", "constant", "power"),
                           power = NULL,
                           moments = c("estimate", "normal", "gamma", "student_t", "custom"),
                           shape = NULL, df = NULL, skewness = NULL, kurtosis = NULL,
                           covariance = c("model", "sandwich"), interval_level = 0.95,
                           km_bounds = NULL, allow_zero_substrate = TRUE,
                           control = list()) {
  moments <- match.arg(moments)
  covariance <- match.arg(covariance)
  if (!is.function(variance)) variance <- match.arg(variance)
  interval_level <- moment_level(interval_level)
  control <- moment_control(control)
  validated <- validate_mm_inputs(x, y, allow_zero_substrate = allow_zero_substrate)
  if (length(validated$y) < 10) {
    stop("Higher-moment fitting requires at least 10 observations.", call. = FALSE)
  }
  pilot <- moment_profile_pilot(validated$x, validated$y, variance, power, interval_level,
                                 km_bounds, allow_zero_substrate)
  reference <- c(unname(pilot$coefficients), pilot$normalized_scale)
  if (any(!is.finite(reference)) || any(reference <= 0)) {
    stop("Higher-moment fitting requires a positive, nondegenerate pilot fit.", call. = FALSE)
  }
  h <- evaluate_variance_function(pilot$variance_fun, pilot$x) / pilot$variance_normalization
  normal_state <- function(eta) {
    moment_score_state(eta, pilot$x, pilot$y, h, reference, 0, 3)
  }
  scale <- sqrt(diag(normal_state(rep(0, 3))$information))
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Insufficient information for higher-moment fitting.", call. = FALSE)
  }
  lower <- c(-30, log(pilot$km_bounds[1] / reference[2]), -30)
  upper <- c(30, log(pilot$km_bounds[2] / reference[2]), 30)
  start <- moment_newton(function(eta) normal_state(eta)$equation / scale,
                          rep(0, 3), lower, upper, control)
  if (!start$ok) {
    moment_fit_error("initial mean/variance fitting", start, moments, control)
  }
  initial <- normal_state(start$eta)
  pilot$residuals <- pilot$y - initial$mean
  pilot$variance_hat <- initial$variance
  specs <- prepare_error_moments(pilot, moments, shape, df, skewness, kurtosis)
  state <- function(eta, fraction = 1) {
    moment_score_state(eta, pilot$x, pilot$y, h, reference,
                        fraction * specs$skewness,
                        3 + fraction * (specs$kurtosis - 3))
  }
  eta <- start$eta
  fraction <- 0
  step <- 0.25
  iterations <- start$iterations
  stages <- 0L
  # Continuation follows the pilot root rather than choosing an unrelated score root.
  while (fraction < 1) {
    target <- min(1, fraction + step)
    attempt <- moment_newton(function(par) state(par, target)$equation / scale,
                              eta, lower, upper, control)
    iterations <- iterations + attempt$iterations
    if (!attempt$ok) {
      step <- step / 2
      if (step < control$min_step) {
        moment_fit_error("higher-moment continuation", attempt, moments, control,
                         iterations = iterations, progress = fraction)
      }
    } else {
      eta <- attempt$eta
      fraction <- target
      step <- min(0.25, 2 * step)
      stages <- stages + 1L
    }
  }
  if (any(eta - lower < 1e-6) || any(upper - eta < 1e-6)) {
    moment_fit_error("boundary validation",
                     list(reason = "search_boundary",
                          score_norm = max(abs(state(eta)$equation / scale))),
                     moments, control, iterations = iterations, progress = fraction)
  }
  final <- state(eta)
  bread <- -moment_jacobian(function(par) state(par)$equation, eta)
  inv_bread <- moment_inverse(bread, "Estimating-equation derivative")
  inv_information <- moment_inverse(final$information, "Working information")
  n <- length(pilot$y)
  meat <- crossprod(final$scores) / n
  sandwich_eta <- inv_bread %*% meat %*% t(inv_bread) / n
  model_eta <- inv_information / n
  theta <- final$theta
  theta[3] <- theta[3] / pilot$variance_normalization
  names(theta) <- c("Vmax", "Km", "gamma")
  transform <- diag(theta)
  transform_cov <- function(value) {
    out <- transform %*% value %*% transform
    out <- (out + t(out)) / 2
    dimnames(out) <- list(names(theta), names(theta))
    if (any(!is.finite(out)) || any(diag(out) <= 0)) {
      stop("Invalid higher-moment covariance; intervals are not reported.", call. = FALSE)
    }
    out
  }
  sandwich_cov <- transform_cov(sandwich_eta)
  model_cov <- transform_cov(model_eta)
  full_cov <- if (covariance == "sandwich") sandwich_cov else model_cov
  fit <- pilot
  fit$coefficients <- theta[1:2]
  fit$scale <- unname(theta[3])
  fit$normalized_scale <- unname(final$theta[3])
  fit$fitted.values <- final$mean
  fit$residuals <- fit$y - final$mean
  fit$variance_hat <- final$variance
  fit$vcov <- full_cov[1:2, 1:2, drop = FALSE]
  fit$vcov_joint <- full_cov
  fit$vcov_joint_sandwich <- sandwich_cov
  fit$vcov_joint_model <- model_cov
  fit$covariance_type <- covariance
  fit$moment_spec <- specs$spec
  fit$standardized_moments <- specs[c("skewness", "kurtosis")]
  fit$estimating_scores_log <- final$scores
  fit$information_log <- final$information
  fit$bread_log <- bread
  fit$control <- control
  fit$convergence <- list(converged = TRUE, score_norm = max(abs(final$equation / scale)),
                          iterations = iterations, stages = stages,
                          bread_rcond = rcond(bread))
  fit$method <- "higher-moment estimating equations"
  fit$call <- match.call()
  fit$rss <- sum(fit$residuals^2)
  fit$rmse <- sqrt(mean(fit$residuals^2))
  fit$weighted_rss <- sum(fit$residuals^2 / fit$variance_hat)
  # These roots do not maximize the Gaussian working likelihood.
  fit$quasi_loglik <- fit$quasi_aic <- fit$quasi_bic <- NA_real_
  class(fit) <- c("fit_mm_moments", "fit_mm")
  if (moments == "estimate" && n < 50) {
    warning("Third/fourth-moment estimation is unstable in small samples; n < 50 is exploratory, not a coverage guarantee.",
            call. = FALSE)
  }
  fit
}

refit_mm_moments <- function(object, x, y) {
  args <- c(list(x = x, y = y,
                 variance = if (object$variance_model == "custom") object$variance_fun else object$variance_model,
                 power = object$variance_power, covariance = object$covariance_type,
                 interval_level = object$interval_level, km_bounds = object$km_bounds,
                 allow_zero_substrate = object$allow_zero_substrate,
                 control = object$control), object$moment_spec)
  do.call(fit_mm_moments, args)
}

#' @export
confint.fit_mm_moments <- function(object, parm = c("Vmax", "Km"),
                                   level = object$interval_level,
                                   method = c("wald", "bootstrap"), B = 399L,
                                   bootstrap_type = "pairs",
                                   bootstrap_ci = c("studentized", "basic", "percentile"),
                                   ...) {
  method <- match.arg(method)
  level <- moment_level(level)
  bootstrap_ci <- match.arg(bootstrap_ci)
  if (!identical(bootstrap_type, "pairs")) {
    stop("Higher-moment bootstrap inference supports pairs resampling only.", call. = FALSE)
  }
  confint.fit_mm(object, parm = parm, level = level, method = method, B = B,
                 bootstrap_type = "pairs", bootstrap_ci = bootstrap_ci, ...)
}

#' @export
summary.fit_mm_moments <- function(object, level = object$interval_level,
                                   method = c("wald", "bootstrap"), B = 399L,
                                   bootstrap_type = "pairs",
                                   bootstrap_ci = c("studentized", "basic", "percentile"),
                                   ...) {
  method <- match.arg(method)
  level <- moment_level(level)
  bootstrap_ci <- match.arg(bootstrap_ci)
  if (!identical(bootstrap_type, "pairs")) {
    stop("Higher-moment bootstrap inference supports pairs resampling only.", call. = FALSE)
  }
  out <- summary.fit_mm(object, level = level, method = method, B = B,
                        bootstrap_type = "pairs", bootstrap_ci = bootstrap_ci, ...)
  out$covariance_type <- object$covariance_type
  out$moment_spec <- object$moment_spec
  out$convergence <- object$convergence
  class(out) <- c("summary.fit_mm_moments", "summary.fit_mm")
  out
}

#' @export
print.summary.fit_mm_moments <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\nWorking variance:", x$variance_label,
      "\nWorking moments:", x$moment_spec$moments,
      "\nCovariance:", x$covariance_type, "\nCI method:", x$ci_method, "\n\n")
  print(round(x$coefficients, digits))
  if (identical(x$ci_method, "disabled")) cat("\n", x$small_sample_policy$summary_message, "\n")
  cat("\nVariance scale:", signif(x$scale, digits), " RMSE:", signif(x$rmse, digits),
      "\nScaled score norm:", signif(x$convergence$score_norm, digits), "\n")
  invisible(x)
}

#' @export
print.fit_mm_moments <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}

#' @export
predict.fit_mm_moments <- function(object, newdata = NULL, se.fit = FALSE,
                                   interval = c("none", "confidence", "prediction"),
                                   level = object$interval_level, ...) {
  interval <- match.arg(interval)
  level <- moment_level(level)
  if (interval == "prediction") {
    stop("Higher moments do not specify predictive quantiles; use interval = 'confidence' for the mean curve.",
         call. = FALSE)
  }
  predict.fit_mm(object, newdata = newdata, se.fit = se.fit, interval = interval,
                 level = level, ...)
}

#' @export
plot.fit_mm_moments <- function(x, interval = TRUE, interval_type = "confidence", ...) {
  if (!identical(interval_type, "confidence")) {
    stop("Higher-moment plots support confidence bands for the mean, not prediction bands.",
         call. = FALSE)
  }
  plot.fit_mm(x, interval = interval, interval_type = interval_type, ...)
}
