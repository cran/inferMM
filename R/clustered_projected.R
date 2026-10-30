cluster_projected_error <- function(stage, detail, diagnostics = list(), reason = stage) {
  scalar <- vapply(diagnostics, function(x) is.numeric(x) && length(x) == 1L, logical(1))
  values <- if (any(scalar)) paste0(
    " Diagnostics: ", paste(paste0(names(diagnostics)[scalar], "=",
      vapply(diagnostics[scalar], function(x) format(x, digits = 5, trim = TRUE), character(1))),
      collapse = ", "), "."
  ) else ""
  stop(structure(
    list(message = paste0(
      "Projected clustered fitting failed during ", stage, ": ", detail,
      values, " No fit or confidence intervals are returned. Check km_bounds, ",
      "within-cluster substrate coverage, the variance shape and independent clusters; ",
      "rescale concentration/response units if needed. ",
      "No automatic fallback to method = 'profile' is performed."
    ), call = NULL, stage = stage, detail = detail, reason = reason,
      diagnostics = diagnostics),
    class = c("inferMM_cluster_projected_error", "inferMM_convergence_error",
              "error", "condition")
  ))
}

cluster_projected_warning <- function(message, reason, diagnostics = list()) {
  warning(structure(list(message = message, call = NULL, reason = reason,
    diagnostics = diagnostics), class = c("inferMM_cluster_projected_warning", "warning", "condition")))
}

cluster_projected_rms <- function(x) {
  scale <- max(abs(x))
  if (!is.finite(scale) || scale == 0) return(scale)
  scale * sqrt(mean((x / scale)^2))
}

cluster_projected_state <- function(km, x_list, y_list, h_list) {
  rows <- lapply(seq_along(x_list), function(i) {
    x <- x_list[[i]]
    y <- y_list[[i]]
    w <- 1 / h_list[[i]]
    z <- x / (km + x)
    dz <- -x / (km + x)^2
    c_i <- sum(w * z^2)
    if (!is.finite(c_i) || c_i <= 0) {
      cluster_projected_error("conditional projection", "a cluster has no positive substrate loading.")
    }
    amplitude <- sum(w * z * y) / c_i
    residual <- y - z * amplitude
    # Remove the latent-amplitude direction before evaluating its shape score.
    dz_perp <- dz - z * sum(w * z * dz) / c_i
    c(W = amplitude, c = c_i, q = sum(w * dz_perp * residual),
      rss = sum(w * residual^2), df = length(x) - 1)
  })
  out <- do.call(rbind, rows)
  w_rms <- cluster_projected_rms(out[, "W"])
  q_rms <- cluster_projected_rms(out[, "q"])
  score <- if (all(is.finite(c(w_rms, q_rms))) && w_rms > 0 && q_rms > 0) {
    mean((out[, "W"] / w_rms) * (out[, "q"] / q_rms))
  } else {
    NA_real_
  }
  list(clusters = out, score = score, objective = sum(out[, "rss"]))
}

cluster_projected_root <- function(state, bounds, score_tol = 1e-7,
                                   reported_bounds = bounds) {
  log_bounds <- log(bounds)
  if (any(!is.finite(log_bounds)) || diff(log_bounds) <= 2e-6) {
    cluster_projected_error("Km search scaling",
      "the normalized search interval is non-finite or too narrow for an interior root.",
      list(km_lower = reported_bounds[1], km_upper = reported_bounds[2]),
      reason = "invalid_scaled_bounds")
  }
  km_scale <- reported_bounds[2] / bounds[2]
  failures <- character()
  refine <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(tryCatch(expr, error = function(e) {
      messages <<- c(messages, conditionMessage(e)); NA_real_
    }), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    if (length(messages)) failures <<- c(failures, paste(unique(messages), collapse = " | "))
    value
  }
  evaluate <- function(log_k) {
    tryCatch({
      s <- state(exp(log_k))
      if (!is.finite(s$score) || !is.finite(s$objective))
        return(list(score = NA_real_, objective = Inf))
      s
    }, error = function(e) {
      list(score = NA_real_, objective = Inf)
    })
  }
  grid <- seq(log_bounds[1], log_bounds[2], length.out = 400L)
  states <- lapply(grid, evaluate)
  values <- vapply(states, function(s) s$score, numeric(1))
  losses <- vapply(states, function(s) s$objective, numeric(1))
  candidates <- grid[is.finite(values) & values == 0]
  for (j in seq_len(length(grid) - 1L)) {
    if (all(is.finite(values[c(j, j + 1L)])) &&
        sign(values[j]) != sign(values[j + 1L])) {
      root <- refine(stats::uniroot(function(k) evaluate(k)$score,
        c(grid[j], grid[j + 1L]), tol = 1e-12)$root)
      candidates <- c(candidates, root)
    }
  }
  # A local objective minimum is only a candidate; it must also solve the score.
  local_min <- which(is.finite(losses[1:398]) & is.finite(losses[2:399]) &
                       is.finite(losses[3:400]) & losses[2:399] <= losses[1:398] &
                       losses[2:399] <= losses[3:400]) + 1L
  for (j in local_min) {
    minimum <- refine(stats::optimize(function(k) evaluate(k)$objective,
      c(grid[j - 1L], grid[j + 1L]), tol = 1e-12)$minimum)
    candidates <- c(candidates, minimum)
  }
  candidates <- unique(candidates[is.finite(candidates)])
  accepted <- lapply(candidates, function(k) {
    s <- evaluate(k)
    if (min(k - log_bounds[1], log_bounds[2] - k) <= 1e-6 ||
        !is.finite(s$score) || abs(s$score) > score_tol ||
        !is.finite(s$objective)) return(NULL)
    list(km = exp(k), state = s)
  })
  accepted <- Filter(Negate(is.null), accepted)
  finite <- which(is.finite(values))
  best_grid <- if (length(finite)) finite[which.min(abs(values[finite]))] else NA_integer_
  diagnostics <- list(km_bounds = reported_bounds, grid_points = length(grid),
    finite_grid_points = length(finite), candidates = length(candidates),
    accepted_roots = length(accepted), refinement_failures = length(failures),
    min_abs_grid_score = if (length(finite)) abs(values[best_grid]) else NA_real_,
    best_grid_km = if (length(finite)) exp(grid[best_grid] + log(km_scale)) else NA_real_,
    refinement_messages = unique(failures))
  if (!length(accepted)) {
    cluster_projected_error("Km root validation", sprintf(
      "no positive interior root passed the standardized score tolerance %.1e over [%.6g, %.6g].",
      score_tol, reported_bounds[1], reported_bounds[2]), diagnostics,
      reason = "no_acceptable_root")
  }
  best <- accepted[[which.min(vapply(accepted, function(z) z$state$objective, numeric(1)))]]
  boundary_loss <- min(losses[c(1L, length(losses))])
  loss_scale <- max(abs(losses[is.finite(losses)]), .Machine$double.eps)
  if (is.finite(boundary_loss) && boundary_loss < best$state$objective -
      1e-10 * loss_scale) {
    cluster_projected_error("Km objective validation",
      "a search boundary has a lower within-cluster loss than the verified interior candidates.",
      c(diagnostics, list(boundary_loss = boundary_loss, interior_loss = best$state$objective)),
      reason = "boundary_dominates")
  }
  best$candidates <- length(candidates)
  best$diagnostics <- diagnostics
  best
}

fit_cluster_mm_projected <- function(cleaned, variance_spec, interval_level = 0.95,
                                     call = NULL, allow_zero_substrate = TRUE,
                                     variance_selection_mode = NULL,
                                     variance_selection = NULL, power_grid = NULL,
                                     candidate_table = NULL, km_bounds = NULL) {
  if (anyNA(cleaned$.cluster_id)) {
    cluster_projected_error("data validation", "cluster identifiers must not be missing.")
  }
  clusters <- split(cleaned, cleaned$.cluster_id, drop = TRUE)
  m <- length(clusters)
  if (m < 5L || any(lengths(lapply(clusters, function(d) d$.S)) < 2L)) {
    cluster_projected_error("data validation",
      "at least five independent clusters and two readings per cluster are required.",
      list(n_clusters = m, min_cluster_size = min(lengths(lapply(clusters, function(d) d$.S)))))
  }
  loading <- vapply(clusters, function(d) any(d$.S > 0), logical(1))
  if (!all(loading)) cluster_projected_error("conditional projection",
    "a cluster has no positive substrate loading.",
    list(n_clusters = m, clusters_without_loading = sum(!loading)))
  informative <- vapply(clusters, function(d) length(unique(d$.S[d$.S > 0])) >= 2L,
                        logical(1))
  if (!any(informative)) {
    cluster_projected_error("identification", "no cluster spans two distinct positive concentrations.")
  }
  if (!is.numeric(interval_level) || length(interval_level) != 1L ||
      !is.finite(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("interval_level must be a finite number strictly between zero and one.", call. = FALSE)
  }
  if (is.null(km_bounds)) km_bounds <- default_km_bounds(cleaned$.S)
  if (!is.numeric(km_bounds) || length(km_bounds) != 2L ||
      any(!is.finite(km_bounds)) || km_bounds[1] <= 0 || km_bounds[2] <= km_bounds[1]) {
    stop("km_bounds must be a finite positive increasing vector of length two.", call. = FALSE)
  }
  x_scale <- max(cleaned$.S)
  y_scale <- max(abs(cleaned$.V))
  if (!is.finite(y_scale) || y_scale <= 0) {
    cluster_projected_error("data validation", "responses are numerically degenerate.")
  }
  x_list <- lapply(clusters, function(d) d$.S / x_scale)
  y_list <- lapply(clusters, function(d) d$.V / y_scale)
  h_evaluated <- lapply(clusters, function(d) variance_spec$fun(d$.S))
  nonfinite_h <- sum(!is.finite(unlist(h_evaluated, use.names = FALSE)))
  if (nonfinite_h > 0L) {
    cluster_projected_error("variance evaluation",
      "the working variance shape is non-finite; rescale concentrations or choose a numerically representable shape.",
      list(nonfinite_variance_values = nonfinite_h), reason = "nonfinite_variance_shape")
  }
  h_raw <- lapply(h_evaluated, safe_pos)
  h_scale <- exp(mean(log(unlist(h_raw, use.names = FALSE))))
  h_list <- lapply(h_raw, function(h) h / h_scale)
  scaled_bounds <- km_bounds / x_scale
  if (any(!is.finite(scaled_bounds)) || any(scaled_bounds <= 0)) {
    cluster_projected_error("Km search scaling",
      "km_bounds cannot be represented after concentration scaling; use compatible units and bounds.",
      list(x_scale = x_scale, km_lower = km_bounds[1], km_upper = km_bounds[2]),
      reason = "invalid_scaled_bounds")
  }
  state <- function(k) cluster_projected_state(k, x_list, y_list, h_list)
  root <- cluster_projected_root(state, scaled_bounds, reported_bounds = km_bounds)
  km <- root$km
  st <- root$state$clusters
  vmax <- mean(st[, "W"])
  gamma <- sum(st[, "rss"]) / sum(st[, "df"])
  if (!all(is.finite(c(vmax, gamma))) || vmax <= 0 ||
      gamma <= 100 * .Machine$double.eps^2) {
    cluster_projected_error("parameter validation", "Vmax or residual scale is non-positive or degenerate.",
      list(scaled_vmax = vmax, scaled_gamma = gamma, n_clusters = m))
  }
  tau_raw <- mean((st[, "W"] - vmax)^2 - gamma / st[, "c"])
  tau_scale <- max(mean(st[, "W"]^2), gamma * mean(1 / st[, "c"]))
  equation_scale <- c(
    vmax,
    sqrt(mean(st[, "W"]^2)) * sqrt(mean((km * st[, "q"])^2)),
    gamma * mean(st[, "df"]), tau_scale
  )
  if (any(!is.finite(equation_scale)) || any(equation_scale <= 0)) {
    cluster_projected_error("score scaling", "the estimating equations are degenerate.")
  }
  eta <- c(log(vmax), log(km), log(gamma), tau_raw / tau_scale)
  scores <- function(par) {
    v <- exp(par[1])
    k <- exp(par[2])
    g <- exp(par[3])
    t_raw <- par[4] * tau_scale
    s <- state(k)$clusters
    out <- cbind(s[, "W"] - v, k * s[, "W"] * s[, "q"],
                 s[, "rss"] - g * s[, "df"],
                 (s[, "W"] - v)^2 - g / s[, "c"] - t_raw)
    sweep(out, 2L, equation_scale, "/")
  }
  psi <- scores(eta)
  score_norm <- max(abs(colMeans(psi)))
  if (!is.finite(score_norm) || score_norm > 1e-7) {
    cluster_projected_error("joint score validation", "the scaled joint estimating equations are not solved.",
      list(score_norm = score_norm, score_tolerance = 1e-7), reason = "joint_score_residual")
  }
  bread <- tryCatch(-moment_jacobian(function(par) colMeans(scores(par)), eta),
    error = function(e) cluster_projected_error("covariance derivative", conditionMessage(e),
      list(n_clusters = m), reason = "derivative_evaluation"))
  if (is.finite(bread[2, 2]) && bread[2, 2] <= 0) {
    cluster_projected_error("Km curvature validation", "the shape-score root is not a strict local loss minimum.",
      list(scaled_curvature = bread[2, 2]), reason = "nonpositive_curvature")
  }
  if (any(!is.finite(bread)) || !is.finite(rcond(bread)) || rcond(bread) < 1e-12) {
    cluster_projected_error("covariance", "the scaled estimating-equation derivative is singular or poorly conditioned.",
      list(bread_rcond = if (all(is.finite(bread))) rcond(bread) else NA_real_,
        minimum_rcond = 1e-12, n_clusters = m), reason = "singular_bread")
  }
  inv_bread <- tryCatch(solve(bread), error = function(e)
    cluster_projected_error("covariance inversion", conditionMessage(e),
      list(bread_rcond = rcond(bread)), reason = "bread_inversion"))
  centered <- sweep(psi, 2L, colMeans(psi), "-")
  meat <- crossprod(centered) / m
  covariance_eta <- inv_bread %*% meat %*% t(inv_bread) / m
  transformation <- c(vmax * y_scale, km * x_scale,
    exp(log(gamma) + 2 * log(y_scale) - log(h_scale)),
    exp(log(tau_scale) + 2 * log(y_scale)))
  if (any(!is.finite(transformation)) || any(transformation <= 0)) {
    cluster_projected_error("original-unit scaling",
      "parameter scales cannot be represented as finite positive values; rescale the assay units.",
      list(x_scale = x_scale, y_scale = y_scale), reason = "unrepresentable_parameter_scale")
  }
  transform <- diag(transformation)
  covariance <- transform %*% covariance_eta %*% transform
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(c("Vmax", "Km", "gamma", "tau2_unconstrained"),
                               c("Vmax", "Km", "gamma", "tau2_unconstrained"))
  beta_cov <- covariance[1:2, 1:2, drop = FALSE]
  beta_se <- sqrt(pmax(diag(beta_cov), 0))
  beta_cor <- sweep(sweep(beta_cov, 1L, beta_se, "/"), 2L, beta_se, "/")
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0) ||
      any(!is.finite(beta_cor)) || min(eigen(beta_cor, symmetric = TRUE, only.values = TRUE)$values) <= 1e-12) {
    cluster_projected_error("covariance", "the mean-parameter sandwich covariance is non-finite or singular.",
      list(n_clusters = m, bread_rcond = rcond(bread)), reason = "invalid_mean_covariance")
  }
  beta <- c(Vmax = vmax * y_scale, Km = km * x_scale)
  gamma <- transformation[3]
  tau_raw <- if (tau_raw == 0) 0 else sign(tau_raw) * exp(log(abs(tau_raw)) + 2 * log(y_scale))
  tau2 <- max(tau_raw, 0)
  z <- cleaned$.S / (beta["Km"] + cleaned$.S)
  within_var <- gamma * safe_pos(variance_spec$fun(cleaned$.S))
  marginal_var <- tau2 * z^2 + within_var
  fitted <- beta["Vmax"] * z
  reference <- unique(cleaned[, c(".cluster_id", ".cluster_label"), drop = FALSE])
  names(reference) <- c("cluster_id", "cluster_label")
  if (tau_raw <= 0) {
    cluster_projected_warning(sprintf(paste(
      "The projected fit passed its score checks, but the random-effect variance moment",
      "is non-positive (tau2_unconstrained=%.6g) and is truncated to zero.",
      "This truncation is not used in mean estimation or its sandwich covariance."), tau_raw),
      "variance_component_boundary", list(tau2_unconstrained = unname(tau_raw)))
  }
  incomplete_search <- root$diagnostics$finite_grid_points < root$diagnostics$grid_points ||
    root$diagnostics$refinement_failures > 0L
  search_warning <- if (incomplete_search) sprintf(paste(
    "Projected search verified an interior root, but only %d/%d grid scores were finite",
    "and %d refinements failed. Inspect km_bounds and units; this is not a complete search guarantee."),
    root$diagnostics$finite_grid_points, root$diagnostics$grid_points,
    root$diagnostics$refinement_failures) else NULL
  if (incomplete_search) cluster_projected_warning(search_warning, "incomplete_search", root$diagnostics)
  structure(c(list(
    coefficients = beta, vcov = beta_cov, vcov_joint = covariance,
    covariance_method = "cluster-sandwich", estimation_method = "projected",
    scale = unname(gamma), tau2 = unname(tau2), tau2_unconstrained = unname(tau_raw),
    fitted.values = unname(fitted), residuals = unname(cleaned$.V - fitted),
    variance_hat = unname(marginal_var), within_variance = unname(within_var),
    between_variance = unname(tau2 * z^2), x = cleaned$.S, y = cleaned$.V,
    cluster = as.character(cleaned$.cluster_id), cluster_reference = reference,
    n_clusters = m, cluster_sizes = as.integer(table(cleaned$.cluster_id)),
    cleaned_data = cleaned, variance_model = variance_spec$model,
    variance_label = variance_spec$label, variance_fun = variance_spec$fun,
    variance_power = attr(variance_spec$fun, "variance_power", exact = TRUE),
    variance_selection_mode = variance_selection_mode,
    variance_selection = variance_selection, power_grid = power_grid,
    variance_candidate_table = candidate_table,
    method = "clustered-projected-score (Ma-Genton-style, identity instrument)",
    converged = TRUE, iterations = 1L,
    initialization_method = "within-cluster projection; no pooled mean initializer",
    initialization_warning = NULL, convergence_warning = search_warning,
    interval_level = interval_level, allow_zero_substrate = allow_zero_substrate,
    km_bounds = km_bounds, call = call, score_contributions = psi,
    bread = bread, meat = meat,
    convergence = list(score_norm = abs(root$state$score), joint_score_norm = score_norm, tolerance = 1e-7,
      candidates = root$candidates, bread_rcond = rcond(bread),
      search = root$diagnostics, search_warning = search_warning,
      variance_component_boundary = tau_raw <= 0),
    assumptions = c("Conditionally independent Gaussian residuals with variance gamma h(S).",
      "Independent clusters; E(b | S) = 0 and common finite random-effect moments.",
      "Random-effect density unspecified; no general efficiency guarantee.")
  ), quasi_metrics(cleaned$.V, fitted, marginal_var, k = 4L)), class = "cluster_mm")
}
