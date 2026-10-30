profile_fit_error <- function(stage, reason, bounds, score_norm = NA_real_,
                              tolerance = 1e-7, candidates = 0L,
                              information_rcond = NA_real_, detail = NULL) {
  reasons <- c(
    non_finite_profile = "the profile calculation is non-finite or degenerate",
    no_verified_root = "no acceptable positive interior Km root could be verified",
    search_boundary = "the candidate solution is on or too close to a search boundary",
    non_positive_parameters = "the candidate Vmax or variance scale is not positive and finite",
    degenerate_variance = "the residual variance is numerically zero",
    singular_information = "the scaled mean information matrix is singular or poorly conditioned",
    invalid_covariance = "the covariance calculation is non-finite or non-positive",
    no_successful_candidate = "none of the candidate variance models produced a verified fit"
  )
  description <- unname(reasons[reason])
  if (is.na(description)) description <- reason
  message <- paste(
    sprintf("Michaelis-Menten profile fitting failed during %s: %s.", stage, description),
    sprintf("Km search bounds: [%.6g, %.6g]; candidates checked: %d; max absolute standardized score: %.3g; tolerance: %.3g.",
            bounds[1], bounds[2], candidates, score_norm, tolerance),
    "No fit or confidence intervals are returned.",
    "Review km_bounds, substrate coverage, outliers, the variance function, and parameter identifiability."
  )
  if (is.finite(information_rcond)) {
    message <- paste(message, sprintf("Scaled information reciprocal condition number: %.3g.", information_rcond))
  }
  if (!is.null(detail)) message <- paste(message, "Details:", detail)
  stop(structure(
    list(message = message, call = NULL, stage = stage, reason = reason,
         km_bounds = bounds, score_norm = score_norm, tolerance = tolerance,
         candidates = candidates, information_rcond = information_rcond,
         detail = detail),
    class = c("inferMM_profile_error", "inferMM_convergence_error", "error", "condition")
  ))
}

profile_problem <- function(x, y, h_vals, bounds) {
  x_scale <- max(x)
  y_scale <- max(abs(y))
  w <- 1 / h_vals
  if (!is.finite(x_scale) || x_scale <= 0 || !is.finite(y_scale) || y_scale <= 0 ||
      any(!is.finite(w)) || any(w <= 0)) {
    profile_fit_error("profile preparation", "non_finite_profile", bounds)
  }
  w <- w / max(w)
  if (any(!is.finite(w)) || any(w <= 0)) {
    profile_fit_error("profile preparation", "non_finite_profile", bounds,
                      detail = "The variance weights have too large a dynamic range.")
  }
  list(x = x / x_scale, y = y / y_scale, w = w,
       x_scale = x_scale, y_scale = y_scale,
       log_bounds = log(bounds) - log(x_scale))
}

profile_score_state <- function(log_k, problem) {
  k <- exp(log_k)
  z <- problem$x / (k + problem$x)
  u <- z / (k + problem$x)
  b <- sum(problem$w * z^2)
  if (!is.finite(b) || b <= 0) {
    return(list(score = NA_real_, score_norm = Inf, vmax = NA_real_, objective = Inf))
  }
  vmax <- sum(problem$w * z * problem$y) / b
  residuals <- problem$y - vmax * z
  objective <- sum(problem$w * residuals^2)
  # Project out the eliminated Vmax direction before scaling the Km score.
  projected <- u - z * sum(problem$w * z * u) / b
  curvature_scale <- sum(problem$w * projected^2)
  zero_variance <- objective <= 100 * .Machine$double.eps^2 *
    sum(problem$w * problem$y^2)
  if (!is.finite(objective) || !is.finite(curvature_scale) || curvature_scale <= 0) {
    score <- NA_real_
    score_norm <- Inf
  } else if (zero_variance) {
    score <- 0
    score_norm <- 0
  } else {
    score <- sum(problem$w * projected * residuals) /
      (sqrt(curvature_scale) * sqrt(objective))
    vmax_score <- sum(problem$w * z * residuals) / (sqrt(b) * sqrt(objective))
    score_norm <- max(abs(c(score, vmax_score)))
  }
  list(score = score, score_norm = score_norm, vmax = vmax,
       objective = objective, zero_variance = zero_variance)
}

find_km_root <- function(x, y, h_vals, bounds, grid_n = 400L,
                         score_tol = 1e-7, boundary_tol = 1e-6) {
  problem <- profile_problem(x, y, h_vals, bounds)
  state <- function(log_k) profile_score_state(log_k, problem)
  score <- function(log_k) state(log_k)$score
  grid <- seq(problem$log_bounds[1], problem$log_bounds[2], length.out = grid_n)
  values <- vapply(grid, score, numeric(1))
  if (!any(is.finite(values))) {
    profile_fit_error("Km root search", "non_finite_profile", bounds,
                      tolerance = score_tol)
  }
  candidates <- grid[is.finite(values) & values == 0]
  methods <- rep("grid-root", length(candidates))
  endpoints <- c(1L, grid_n)
  near_boundary <- endpoints[is.finite(values[endpoints]) &
                              values[endpoints] != 0 & abs(values[endpoints]) <= score_tol]
  candidates <- c(candidates, grid[near_boundary])
  methods <- c(methods, rep("boundary-candidate", length(near_boundary)))
  finite_pair <- is.finite(values[-grid_n]) & is.finite(values[-1])
  brackets <- which(finite_pair & values[-grid_n] != 0 & values[-1] != 0 &
                      sign(values[-grid_n]) != sign(values[-1]))
  for (i in brackets) {
    root <- tryCatch(stats::uniroot(score, interval = grid[c(i, i + 1L)],
                                    tol = 1e-12, maxiter = 1000L,
                                    check.conv = TRUE)$root,
                     error = function(e) NA_real_, warning = function(w) NA_real_)
    if (is.finite(root)) {
      candidates <- c(candidates, root)
      methods <- c(methods, "bracketed-root")
    }
  }
  objective <- function(log_k) {
    value <- score(log_k)
    if (is.finite(value)) value^2 else Inf
  }
  acceptable_root <- function(k) {
    value <- state(k)
    is.finite(value$score_norm) && value$score_norm <= score_tol &&
      is.finite(value$vmax) && value$vmax > 0 &&
      k - problem$log_bounds[1] > boundary_tol &&
      problem$log_bounds[2] - k > boundary_tol
  }
  if (!length(candidates) || !any(vapply(candidates, acceptable_root, logical(1)))) {
    # Minimization supplies candidates only; a small objective is not a fit.
    local <- which(is.finite(values[2:(grid_n - 1L)]) &
                     abs(values[2:(grid_n - 1L)]) < abs(values[1:(grid_n - 2L)]) &
                     abs(values[2:(grid_n - 1L)]) < abs(values[3:grid_n])) + 1L
    intervals <- c(list(problem$log_bounds), lapply(local, function(i) grid[c(i - 1L, i + 1L)]))
    for (interval in intervals) {
      candidate <- tryCatch(stats::optimize(objective, interval, tol = 1e-12)$minimum,
                            error = function(e) NA_real_, warning = function(w) NA_real_)
      if (is.finite(candidate)) {
        candidates <- c(candidates, candidate)
        methods <- c(methods, "verified-score-minimum")
      }
    }
  }
  states <- lapply(candidates, state)
  norms <- vapply(states, function(s) s$score_norm, numeric(1))
  verified <- is.finite(norms) & norms <= score_tol
  positive <- vapply(states, function(s) is.finite(s$vmax) && s$vmax > 0 &&
                      is.finite(s$objective), logical(1))
  interior <- candidates - problem$log_bounds[1] > boundary_tol &
    problem$log_bounds[2] - candidates > boundary_tol
  valid <- verified & positive & interior
  if (!any(valid)) {
    reason <- if (any(verified & positive & !interior)) {
      "search_boundary"
    } else if (any(verified & !positive)) {
      "non_positive_parameters"
    } else {
      "no_verified_root"
    }
    profile_fit_error("Km root validation", reason, bounds,
                      score_norm = if (length(norms)) min(norms) else Inf,
                      tolerance = score_tol, candidates = length(candidates))
  }
  valid_idx <- which(valid)
  selected <- valid_idx[which.min(vapply(states[valid_idx], function(s) s$objective, numeric(1)))]
  if (isTRUE(states[[selected]]$zero_variance)) {
    profile_fit_error("variance-scale validation", "degenerate_variance", bounds,
                      score_norm = norms[selected], tolerance = score_tol,
                      candidates = length(candidates))
  }
  km <- exp(candidates[selected] + log(problem$x_scale))
  vmax <- states[[selected]]$vmax * problem$y_scale
  if (any(!is.finite(c(km, vmax))) || km <= 0 || vmax <= 0) {
    profile_fit_error("parameter validation", "non_positive_parameters", bounds,
                      score_norm = norms[selected], tolerance = score_tol,
                      candidates = length(candidates))
  }
  candidate_table <- data.frame(
    Km = exp(candidates + log(problem$x_scale)),
    Vmax = vapply(states, function(s) s$vmax, numeric(1)) * problem$y_scale,
    standardized_score = norms,
    scaled_weighted_rss = vapply(states, function(s) s$objective, numeric(1)),
    accepted = valid, selected = seq_along(candidates) == selected,
    method = methods
  )
  list(root = km, vmax = vmax, problem = problem,
       convergence = list(converged = TRUE, score_norm = norms[selected],
                          tolerance = score_tol, method = methods[selected],
                          candidate_count = length(candidates),
                          verified_roots = sum(valid), candidates = candidate_table))
}

profile_covariance <- function(grad, h_vals, scale, bounds, convergence) {
  weighted_grad <- grad / sqrt(h_vals)
  column_scale <- apply(abs(weighted_grad), 2, max)
  if (any(!is.finite(weighted_grad)) || any(!is.finite(column_scale)) ||
      any(column_scale <= 0)) {
    profile_fit_error("covariance calculation", "singular_information", bounds,
                      convergence$score_norm, convergence$tolerance,
                      convergence$candidate_count)
  }
  standardized_grad <- sweep(weighted_grad, 2, column_scale, "/")
  info <- crossprod(standardized_grad) / nrow(grad)
  condition <- rcond(info)
  if (!is.finite(condition) || condition < 1e-12) {
    profile_fit_error("covariance calculation", "singular_information", bounds,
                      convergence$score_norm, convergence$tolerance,
                      convergence$candidate_count, information_rcond = condition)
  }
  inverse <- tryCatch(solve(info), error = function(e) e)
  if (inherits(inverse, "error")) {
    profile_fit_error("covariance calculation", "singular_information", bounds,
                      convergence$score_norm, convergence$tolerance,
                      convergence$candidate_count, information_rcond = condition,
                      detail = conditionMessage(inverse))
  }
  standard_error_scale <- sqrt(scale / nrow(grad)) / column_scale
  covariance <- inverse * outer(standard_error_scale, standard_error_scale)
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    profile_fit_error("covariance calculation", "invalid_covariance", bounds,
                      convergence$score_norm, convergence$tolerance,
                      convergence$candidate_count, information_rcond = condition)
  }
  list(covariance = covariance, information_rcond = condition)
}
