prepare_cluster_mm_data <- function(data,
                                    s,
                                    v,
                                    cluster,
                                    invalid_response_values = c(-9999),
                                    invalid_substrate_values = numeric(0),
                                    allow_zero_substrate = TRUE,
                                    sort_x = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  if (missing(cluster) || !length(cluster)) {
    stop("cluster must name one or more clustering columns.")
  }

  needed <- c(s, v, cluster)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  dat <- data
  dat$.row_id <- seq_len(nrow(dat))
  dat$.S <- suppressWarnings(as.numeric(dat[[s]]))
  dat$.V <- suppressWarnings(as.numeric(dat[[v]]))

  keep <- is.finite(dat$.S) & is.finite(dat$.V)
  if (length(invalid_substrate_values)) {
    keep <- keep & !(dat$.S %in% invalid_substrate_values)
  }
  if (length(invalid_response_values)) {
    keep <- keep & !(dat$.V %in% invalid_response_values)
  }
  if (isTRUE(allow_zero_substrate)) {
    keep <- keep & dat$.S >= 0
  } else {
    keep <- keep & dat$.S > 0
  }

  dat <- dat[keep, , drop = FALSE]
  if (!nrow(dat)) {
    stop("No valid observations remain after cleaning.")
  }

  cluster_parts <- lapply(cluster, function(col) as.character(dat[[col]]))
  dat$.cluster_id <- interaction(cluster_parts, drop = TRUE, lex.order = TRUE)
  dat$.cluster_label <- apply(
    as.data.frame(cluster_parts, stringsAsFactors = FALSE),
    1,
    function(z) paste(paste(cluster, z, sep = "="), collapse = " | ")
  )

  if (isTRUE(sort_x)) {
    ord <- order(dat$.cluster_id, dat$.S, dat$.row_id)
    dat <- dat[ord, , drop = FALSE]
  }

  row.names(dat) <- NULL

  if (length(unique(dat$.cluster_id)) < 2L) {
    stop("At least two clusters are required for cluster_mm().")
  }

  if (length(unique(dat$.S)) < 3L) {
    stop("At least 3 distinct substrate concentrations are required.")
  }

  dat
}

working_cov_inv_re <- function(z, h_vals, tau2, gamma) {
  d_inv <- 1 / safe_pos(gamma * h_vals)
  denom <- 1 + tau2 * sum(z^2 * d_inv)
  outer_piece <- outer(z * d_inv, z * d_inv)
  diag(d_inv) - (tau2 / denom) * outer_piece
}

working_logdet_re <- function(z, h_vals, tau2, gamma) {
  sum(log(safe_pos(gamma * h_vals))) +
    log1p(tau2 * sum(z^2 / safe_pos(gamma * h_vals)))
}

update_tau2_gamma_clustered <- function(resid_list, z_list, h_list) {
  x_off_all <- numeric(0)
  y_off_all <- numeric(0)

  for (i in seq_along(resid_list)) {
    r_i <- as.numeric(resid_list[[i]])
    z_i <- as.numeric(z_list[[i]])
    n_i <- length(r_i)

    if (n_i <= 1L) {
      next
    }

    idx <- which(lower.tri(matrix(0, n_i, n_i)), arr.ind = TRUE)
    x_off_all <- c(x_off_all, z_i[idx[, 1]] * z_i[idx[, 2]])
    y_off_all <- c(y_off_all, outer(r_i, r_i)[lower.tri(matrix(0, n_i, n_i))])
  }

  tau2_hat <- if (length(x_off_all) && sum(x_off_all^2) > 0) {
    sum(x_off_all * y_off_all) / sum(x_off_all^2)
  } else {
    1e-8
  }
  tau2_hat <- max(unname(tau2_hat), 1e-8)

  diag_target <- unlist(
    Map(
      function(r_i, z_i) r_i^2 - tau2_hat * z_i^2,
      resid_list,
      z_list
    ),
    use.names = FALSE
  )
  h_all <- unlist(h_list, use.names = FALSE)
  gamma_hat <- if (length(h_all) && sum(h_all^2) > 0) {
    sum(h_all * diag_target) / sum(h_all^2)
  } else {
    1e-8
  }
  gamma_hat <- max(unname(gamma_hat), 1e-8)

  c(tau2 = tau2_hat, gamma = gamma_hat)
}

fit_cluster_mm_weighted <- function(cleaned,
                                    variance_spec,
                                    interval_level = 0.95,
                                    max_iter = 8L,
                                    tol = 1e-5,
                                    call = NULL,
                                    allow_zero_substrate = TRUE,
                                    variance_selection_mode = NULL,
                                    variance_selection = NULL,
                                    power_grid = NULL,
                                    candidate_table = NULL) {
  cluster_split <- split(cleaned, cleaned$.cluster_id, drop = TRUE)
  x_list <- lapply(cluster_split, function(dat) dat$.S)
  y_list <- lapply(cluster_split, function(dat) dat$.V)
  h_list <- lapply(x_list, function(x_i) safe_pos(variance_spec$fun(x_i)))
  x_all <- cleaned$.S
  y_all <- cleaned$.V

  pooled_start <- tryCatch(
    fit_mm_from_spec(
      x = x_all,
      y = y_all,
      variance = variance_spec$model,
      power = attr(variance_spec$fun, "variance_power", exact = TRUE),
      interval_level = interval_level,
      km_bounds = default_km_bounds(x_all),
      allow_zero_substrate = allow_zero_substrate
    ),
    error = function(e) e
  )
  initialization_method <- "pooled single-curve profile-score fit"
  initialization_warning <- NULL

  if (inherits(pooled_start, "error")) {
    stop(
      paste(
        "Unable to initialize cluster_mm() from the pooled single-curve",
        "profile-score fit.",
        conditionMessage(pooled_start)
      ),
      call. = FALSE
    )
  }

  Vmax <- unname(coef(pooled_start)["Vmax"])
  Km <- unname(coef(pooled_start)["Km"])
  z_list <- lapply(x_list, function(x_i) x_i / (Km + x_i))
  resid_list <- Map(function(y_i, z_i) y_i - Vmax * z_i, y_list, z_list)
  nuis_start <- update_tau2_gamma_clustered(
    resid_list,
    z_list = z_list,
    h_list = h_list
  )
  tau2_default <- max(stats::var(base::tapply(y_all, cleaned$.cluster_id, mean)), 1e-4)
  gamma_default <- max(pooled_start$scale, 1e-4)
  tau2 <- if (is.finite(nuis_start["tau2"])) max(unname(nuis_start["tau2"]), 1e-4) else tau2_default
  gamma <- if (is.finite(nuis_start["gamma"])) max(unname(nuis_start["gamma"]), 1e-4) else gamma_default
  delta <- Inf
  iter <- 0L
  converged <- FALSE

  objective_for_km <- function(km, tau2, gamma) {
    z_list <- lapply(x_list, function(x_i) x_i / (km + x_i))
    inv_list <- vector("list", length(z_list))
    denom <- 0
    numer <- 0
    logdet <- 0

    for (i in seq_along(z_list)) {
      inv_i <- working_cov_inv_re(z_list[[i]], h_list[[i]], tau2 = tau2, gamma = gamma)
      inv_list[[i]] <- inv_i
      denom <- denom + as.numeric(crossprod(z_list[[i]], inv_i %*% z_list[[i]]))
      numer <- numer + as.numeric(crossprod(z_list[[i]], inv_i %*% y_list[[i]]))
      logdet <- logdet + working_logdet_re(z_list[[i]], h_list[[i]], tau2, gamma)
    }

    if (!is.finite(denom) || denom <= 0 || !is.finite(numer)) {
      return(Inf)
    }

    vmax <- numer / denom
    if (!is.finite(vmax)) {
      return(Inf)
    }
    obj <- 0
    for (i in seq_along(z_list)) {
      resid_i <- y_list[[i]] - vmax * z_list[[i]]
      obj <- obj + as.numeric(crossprod(resid_i, inv_list[[i]] %*% resid_i))
    }

    out <- obj + logdet
    if (!is.finite(out)) {
      return(Inf)
    }
    out
  }

  for (iter in seq_len(max_iter)) {
    km_old <- Km
    vmax_old <- Vmax

    opt <- stats::optimize(
      f = function(log_k) objective_for_km(exp(log_k), tau2 = tau2, gamma = gamma),
      interval = log(default_km_bounds(x_all))
    )
    Km <- exp(opt$minimum)
    if (!is.finite(Km) || Km <= 0) {
      stop("Unable to obtain a finite positive Km update in cluster_mm().")
    }
    z_list <- lapply(x_list, function(x_i) x_i / (Km + x_i))
    inv_list <- lapply(
      seq_along(z_list),
      function(i) working_cov_inv_re(z_list[[i]], h_list[[i]], tau2 = tau2, gamma = gamma)
    )

    denom <- 0
    numer <- 0
    for (i in seq_along(z_list)) {
      denom <- denom + as.numeric(crossprod(z_list[[i]], inv_list[[i]] %*% z_list[[i]]))
      numer <- numer + as.numeric(crossprod(z_list[[i]], inv_list[[i]] %*% y_list[[i]]))
    }
    if (!is.finite(denom) || denom <= 0 || !is.finite(numer)) {
      stop("Clustered Michaelis-Menten update produced a non-finite weighted normal equation.")
    }
    Vmax <- numer / denom
    if (!is.finite(Vmax)) {
      stop("Unable to obtain a finite Vmax update in cluster_mm().")
    }

    resid_list <- Map(function(y_i, z_i) y_i - Vmax * z_i, y_list, z_list)
    nuis <- update_tau2_gamma_clustered(resid_list, z_list = z_list, h_list = h_list)
    tau2 <- nuis["tau2"]
    gamma <- nuis["gamma"]
    if (!is.finite(tau2) || !is.finite(gamma)) {
      stop("Clustered variance-component update produced non-finite values.")
    }

    delta <- max(abs(c(Vmax - vmax_old, Km - km_old)))
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  convergence_warning <- NULL
  if (!converged) {
    convergence_warning <- sprintf(
      paste(
        "cluster_mm() reached the maximum of %d outer iterations",
        "without meeting the tolerance %.2e.",
        "Estimates may still be usable, but treat the fit as numerically fragile."
      ),
      max_iter,
      tol
    )
    warning(convergence_warning, call. = FALSE)
  }

  z_all <- x_all / (Km + x_all)
  within_var <- gamma * safe_pos(variance_spec$fun(x_all))
  marginal_var <- tau2 * z_all^2 + within_var
  fitted <- mm_mean(x_all, Vmax, Km)
  residuals <- y_all - fitted

  info_beta <- matrix(0, nrow = 2L, ncol = 2L)
  for (i in seq_along(cluster_split)) {
    dat_i <- cluster_split[[i]]
    D_i <- mm_grad(dat_i$.S, Vmax, Km)
    z_i <- dat_i$.S / (Km + dat_i$.S)
    inv_i <- working_cov_inv_re(z_i, h_list[[i]], tau2 = tau2, gamma = gamma)
    info_beta <- info_beta + crossprod(D_i, inv_i %*% D_i)
  }
  dimnames(info_beta) <- list(c("Vmax", "Km"), c("Vmax", "Km"))
  vcov_beta <- safe_inverse(info_beta)
  beta <- c(Vmax = unname(Vmax), Km = unname(Km))
  criteria <- quasi_metrics(y_all, fitted, marginal_var, k = 4L)

  cluster_reference <- unique(cleaned[, c(".cluster_id", ".cluster_label"), drop = FALSE])
  names(cluster_reference) <- c("cluster_id", "cluster_label")

  structure(
    c(
      list(
        coefficients = beta,
        vcov = vcov_beta,
        scale = gamma,
        tau2 = tau2,
        fitted.values = fitted,
        residuals = residuals,
        variance_hat = marginal_var,
        within_variance = within_var,
        between_variance = tau2 * z_all^2,
        x = x_all,
        y = y_all,
        cluster = as.character(cleaned$.cluster_id),
        cluster_reference = cluster_reference,
        n_clusters = nrow(cluster_reference),
        cluster_sizes = as.integer(table(cleaned$.cluster_id)),
        cleaned_data = cleaned,
        variance_model = variance_spec$model,
        variance_label = variance_spec$label,
        variance_fun = variance_spec$fun,
        variance_power = attr(variance_spec$fun, "variance_power", exact = TRUE),
        variance_selection_mode = variance_selection_mode,
        variance_selection = variance_selection,
        power_grid = power_grid,
        variance_candidate_table = candidate_table,
        method = "clustered-profile-score",
        converged = converged,
        iterations = iter,
        initialization_method = initialization_method,
        initialization_warning = initialization_warning,
        convergence_warning = convergence_warning,
        interval_level = interval_level,
        allow_zero_substrate = allow_zero_substrate,
        call = call
      ),
      criteria
    ),
    class = "cluster_mm"
  )
}

fit_cluster_mm_from_spec <- function(cleaned,
                                     variance,
                                     power = NULL,
                                     interval_level = 0.95,
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

  fit_cluster_mm_weighted(
    cleaned = cleaned,
    variance_spec = variance_spec,
    interval_level = interval_level,
    call = call,
    allow_zero_substrate = allow_zero_substrate,
    variance_selection_mode = variance_selection_mode,
    variance_selection = variance_selection,
    power_grid = power_grid,
    candidate_table = candidate_table
  )
}

search_auto_cluster_variance <- function(cleaned,
                                         power_selection = c("quasi_aic", "quasi_bic"),
                                         power_grid = default_power_grid(),
                                         interval_level = 0.95,
                                         allow_zero_substrate = TRUE,
                                         call = NULL) {
  power_selection <- match.arg(power_selection)
  power_grid <- validate_power_grid(power_grid)

  candidates <- c(
    list(list(model = "log", power = NULL)),
    lapply(power_grid, function(power) list(model = "power", power = power))
  )

  fits <- lapply(candidates, function(spec) {
    fit_cluster_mm_from_spec(
      cleaned = cleaned,
      variance = spec$model,
      power = spec$power,
      interval_level = interval_level,
      allow_zero_substrate = allow_zero_substrate
    )
  })

  candidate_table <- do.call(
    rbind,
    lapply(fits, function(fit) {
      make_screen_row(
        model = if (identical(fit$variance_model, "power")) {
          power_model_name(fit$variance_power)
        } else {
          fit$variance_model
        },
        label = fit$variance_label,
        success = TRUE,
        runtime = NA_real_,
        fit = fit
      )
    })
  )

  metric_name <- if (identical(power_selection, "quasi_bic")) "quasi_bic" else "quasi_aic"
  best_idx <- which.min(candidate_table[[metric_name]])
  best_fit <- fits[[best_idx]]
  best_fit$variance_selection_mode <- "auto"
  best_fit$variance_selection <- power_selection
  best_fit$power_grid <- power_grid
  best_fit$variance_candidate_table <- candidate_table
  best_fit$variance_label <- auto_selected_label(best_fit$variance_label, power_selection)
  best_fit$call <- call
  best_fit
}

#' Fit a Cluster-Aware Michaelis-Menten Model
#'
#' Estimate Michaelis-Menten mean parameters from repeated or clustered assay
#' data using a Ma-Genton-style working covariance with a random effect on
#' `Vmax` and a low-dimensional residual working variance function.
#'
#' @param data A data frame containing clustered Michaelis-Menten observations.
#' @param s Name of the substrate concentration column.
#' @param v Name of the response column.
#' @param cluster Character vector naming one or more clustering columns.
#' @param variance Character string naming a built-in residual working variance
#'   model. Use `"power"` with a fixed exponent or `"auto"` to compare
#'   `log(S + 1)` against a grid of power functions.
#' @param power Optional exponent used when `variance = "power"`.
#' @param power_selection Criterion used when `variance = "auto"`:
#'   `"quasi_aic"` or `"quasi_bic"`.
#' @param power_grid Candidate exponents considered when `variance = "auto"`.
#' @param interval_level Confidence level used for Wald intervals.
#' @param invalid_response_values Values to remove before fitting.
#' @param invalid_substrate_values Values to remove before fitting.
#' @param allow_zero_substrate Logical; if `FALSE`, substrate concentrations
#'   must be strictly positive.
#' @param sort_x Logical; sort observations within each cluster by substrate
#'   concentration.
#'
#' @return An object of class `"cluster_mm"`.
#' @export
cluster_mm <- function(data,
                       s,
                       v,
                       cluster,
                       variance = c("log", "sqrt", "cuberoot", "constant", "power", "auto"),
                       power = NULL,
                       power_selection = c("quasi_aic", "quasi_bic"),
                       power_grid = default_power_grid(),
                       interval_level = 0.95,
                       invalid_response_values = c(-9999),
                       invalid_substrate_values = numeric(0),
                       allow_zero_substrate = TRUE,
                       sort_x = TRUE) {
  variance <- normalize_variance_model(match.arg(variance))
  power_selection <- match.arg(power_selection)
  cleaned <- prepare_cluster_mm_data(
    data = data,
    s = s,
    v = v,
    cluster = cluster,
    invalid_response_values = invalid_response_values,
    invalid_substrate_values = invalid_substrate_values,
    allow_zero_substrate = allow_zero_substrate,
    sort_x = sort_x
  )

  if (identical(variance, "power")) {
    if (length(power) != 1L || !is.finite(power) || power <= 0) {
      stop("A positive finite scalar 'power' must be supplied when variance = 'power'.")
    }
  } else {
    power <- NULL
  }

  if (identical(variance, "auto")) {
    return(
      search_auto_cluster_variance(
        cleaned = cleaned,
        power_selection = power_selection,
        power_grid = power_grid,
        interval_level = interval_level,
        allow_zero_substrate = allow_zero_substrate,
        call = match.call()
      )
    )
  }

  fit_cluster_mm_from_spec(
    cleaned = cleaned,
    variance = variance,
    power = power,
    interval_level = interval_level,
    allow_zero_substrate = allow_zero_substrate,
    call = match.call(),
    variance_selection_mode = if (identical(variance, "power")) "manual" else NULL,
    variance_selection = if (identical(variance, "power")) "manual" else NULL,
    power_grid = if (identical(variance, "power")) power else NULL
  )
}

refit_cluster_mm <- function(object, cleaned = object$cleaned_data) {
  if (!inherits(object, "cluster_mm")) {
    stop("object must inherit from 'cluster_mm'.")
  }

  if (identical(object$variance_selection_mode, "auto")) {
    return(
      search_auto_cluster_variance(
        cleaned = cleaned,
        power_selection = object$variance_selection,
        power_grid = object$power_grid,
        interval_level = object$interval_level,
        allow_zero_substrate = object$allow_zero_substrate,
        call = object$call
      )
    )
  }

  fit_cluster_mm_from_spec(
    cleaned = cleaned,
    variance = object$variance_model,
    power = object$variance_power,
    interval_level = object$interval_level,
    allow_zero_substrate = object$allow_zero_substrate,
    call = object$call,
    variance_selection_mode = object$variance_selection_mode,
    variance_selection = object$variance_selection,
    power_grid = object$power_grid
  )
}

#' @export
coef.cluster_mm <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.cluster_mm <- function(object, ...) {
  object$vcov
}

small_cluster_ci_policy <- function(object,
                                    min_clusters = 4L,
                                    min_unique_x = 6L) {
  n_clusters <- object$n_clusters
  n_unique_x <- length(unique(object$x))
  disable <- n_clusters < min_clusters || n_unique_x < min_unique_x

  list(
    disable = disable,
    n_clusters = n_clusters,
    n_unique_x = n_unique_x,
    min_clusters = as.integer(min_clusters),
    min_unique_x = as.integer(min_unique_x),
    summary_message = sprintf(
      paste(
        "Interval inference is disabled by default for sparse clustered fits",
        "(%d clusters, %d distinct substrate concentrations).",
        "Default interval reporting requires at least %d clusters and at least %d",
        "distinct concentrations."
      ),
      n_clusters,
      n_unique_x,
      min_clusters,
      min_unique_x
    ),
    confint_warning = sprintf(
      paste(
        "Small-cluster warning: this fit uses %d clusters across %d distinct",
        "substrate concentrations.",
        "Default interval reporting is disabled below %d clusters or below %d",
        "distinct concentrations because coverage may be unstable."
      ),
      n_clusters,
      n_unique_x,
      min_clusters,
      min_unique_x
    )
  )
}

cluster_bootstrap_policy <- function(object,
                                     min_boot_clusters = 8L) {
  n_clusters <- object$n_clusters
  below_recommendation <- n_clusters < min_boot_clusters

  list(
    below_recommendation = below_recommendation,
    n_clusters = n_clusters,
    min_boot_clusters = as.integer(min_boot_clusters),
    summary_message = if (below_recommendation) {
      sprintf(
        paste(
          "Experimental clustered bootstrap sensitivity analysis:",
          "this fit uses %d clusters, below the recommended %d-cluster threshold",
          "for cluster-level pairs resampling."
        ),
        n_clusters,
        min_boot_clusters
      )
    } else {
      sprintf(
        paste(
          "Experimental clustered bootstrap sensitivity analysis:",
          "cluster-level pairs resampling was used with %d clusters.",
          "Treat these intervals as sensitivity checks rather than default inference."
        ),
        n_clusters
      )
    },
    confint_warning = if (below_recommendation) {
      sprintf(
        paste(
          "Experimental clustered bootstrap warning: this fit uses %d clusters,",
          "below the recommended %d-cluster threshold for cluster-level pairs resampling."
        ),
        n_clusters,
        min_boot_clusters
      )
    } else {
      sprintf(
        paste(
          "Experimental clustered bootstrap warning: cluster-level pairs resampling",
          "is provided only as a sensitivity analysis (fit uses %d clusters)."
        ),
        n_clusters
      )
    }
  )
}

#' @export
confint.cluster_mm <- function(object,
                               parm = c("Vmax", "Km"),
                               level = object$interval_level,
                               method = c("wald", "bootstrap"),
                               B = 399L,
                               bootstrap_type = c("pairs"),
                               bootstrap_ci = c("studentized", "basic", "percentile"),
                               seed = NULL,
                               ...) {
  method <- match.arg(method)
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  policy <- small_cluster_ci_policy(object)
  warn_parts <- character(0)

  if (isTRUE(policy$disable)) {
    warn_parts <- c(warn_parts, policy$confint_warning)
  }

  interval <- if (identical(method, "wald")) {
    wald_interval(coef(object), vcov(object), level = level)
  } else {
    boot_policy <- cluster_bootstrap_policy(object)
    warn_parts <- c(warn_parts, boot_policy$confint_warning)
    cluster_bootstrap_interval(
      object = object,
      level = level,
      B = B,
      bootstrap_type = bootstrap_type,
      bootstrap_ci = bootstrap_ci,
      seed = seed
    )
  }

  if (length(warn_parts)) {
    warning(
      paste(c(warn_parts, "Use these intervals only as a sensitivity analysis."), collapse = " "),
      call. = FALSE
    )
  }

  out <- cbind(interval$lower, interval$upper)
  colnames(out) <- c("lower", "upper")
  out[parm, , drop = FALSE]
}

#' @export
print.cluster_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  cat("Clusters:", x$n_clusters, "\n")
  cat("Initialization:", x$initialization_method, "\n")
  cat(
    "Converged:", if (isTRUE(x$converged)) "yes" else "no",
    "  Iterations:", x$iterations, "\n"
  )
  cat("Residual working variance:", x$variance_label, "\n\n")
  policy <- small_cluster_ci_policy(x)

  coef_table <- if (isTRUE(policy$disable)) {
    estimate_only_table(x)
  } else {
    coefficient_table(x, level = x$interval_level, method = "wald")
  }
  print(round(coef_table, digits = digits))

  if (isTRUE(policy$disable)) {
    cat("\n", policy$summary_message, "\n", sep = "")
  }
  if (!is.null(x$initialization_warning)) {
    cat("\nInitialization note:", x$initialization_warning, "\n")
  }
  if (!isTRUE(x$converged) && !is.null(x$convergence_warning)) {
    cat("\nConvergence note:", x$convergence_warning, "\n")
  }

  cat(
    "\nTau^2:", round(x$tau2, digits),
    "  Gamma:", round(x$scale, digits),
    "  RMSE:", round(x$rmse, digits), "\n",
    sep = ""
  )
  cat(
    "Quasi-AIC:", round(x$quasi_aic, digits),
    "  Quasi-BIC:", round(x$quasi_bic, digits), "\n",
    sep = ""
  )

  invisible(x)
}

#' @export
summary.cluster_mm <- function(object,
                               level = object$interval_level,
                               method = c("wald", "bootstrap"),
                               B = 399L,
                               bootstrap_type = c("pairs"),
                               bootstrap_ci = c("studentized", "basic", "percentile"),
                               seed = NULL,
                               ...) {
  method <- match.arg(method)
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  policy <- small_cluster_ci_policy(object)
  boot_policy <- if (!isTRUE(policy$disable) && identical(method, "bootstrap")) {
    cluster_bootstrap_policy(object)
  } else {
    NULL
  }

  out <- list(
    call = object$call,
    method = object$method,
    n_clusters = object$n_clusters,
    converged = object$converged,
    iterations = object$iterations,
    initialization_method = object$initialization_method,
    initialization_warning = object$initialization_warning,
    convergence_warning = object$convergence_warning,
    variance_model = object$variance_model,
    variance_label = object$variance_label,
    coefficients = if (isTRUE(policy$disable)) {
      estimate_only_table(object)
    } else {
      coefficient_table(
        object,
        level = level,
        method = method,
        B = B,
        bootstrap_type = bootstrap_type,
        bootstrap_ci = bootstrap_ci,
        seed = seed
      )
    },
    scale = object$scale,
    tau2 = object$tau2,
    rss = object$rss,
    rmse = object$rmse,
    weighted_rss = object$weighted_rss,
    quasi_aic = object$quasi_aic,
    quasi_bic = object$quasi_bic,
    level = level,
    ci_method = if (isTRUE(policy$disable)) "disabled" else method,
    bootstrap_type = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) bootstrap_type else NULL,
    bootstrap_ci = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) bootstrap_ci else NULL,
    bootstrap_B = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) B else NULL,
    small_cluster_policy = policy,
    cluster_bootstrap_policy = boot_policy
  )
  class(out) <- "summary.cluster_mm"
  out
}

#' @export
print.summary.cluster_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  cat("Clusters:", x$n_clusters, "\n")
  cat("Initialization:", x$initialization_method, "\n")
  cat(
    "Converged:", if (isTRUE(x$converged)) "yes" else "no",
    "  Iterations:", x$iterations, "\n"
  )
  cat("Residual working variance:", x$variance_label, "\n\n")
  if (identical(x$ci_method, "disabled")) {
    cat("CI method: disabled by default for sparse clustered fits\n\n")
  } else if (identical(x$ci_method, "bootstrap")) {
    cat(
      "CI method: bootstrap (cluster pairs, ",
      x$bootstrap_ci,
      ", B = ",
      x$bootstrap_B,
      ")\n\n",
      sep = ""
    )
  } else {
    cat("CI method: wald\n\n")
  }
  print(round(x$coefficients, digits = digits))
  if (identical(x$ci_method, "disabled") && !is.null(x$small_cluster_policy)) {
    cat("\n", x$small_cluster_policy$summary_message, "\n", sep = "")
  } else if (identical(x$ci_method, "bootstrap") && !is.null(x$cluster_bootstrap_policy)) {
    cat("\n", x$cluster_bootstrap_policy$summary_message, "\n", sep = "")
  }
  if (!is.null(x$initialization_warning)) {
    cat("\nInitialization note:", x$initialization_warning, "\n")
  }
  if (!isTRUE(x$converged) && !is.null(x$convergence_warning)) {
    cat("\nConvergence note:", x$convergence_warning, "\n")
  }
  cat(
    "\nTau^2:", round(x$tau2, digits),
    "  Gamma:", round(x$scale, digits),
    "  RMSE:", round(x$rmse, digits), "\n",
    sep = ""
  )
  cat(
    "Quasi-AIC:", round(x$quasi_aic, digits),
    "  Quasi-BIC:", round(x$quasi_bic, digits), "\n",
    sep = ""
  )

  invisible(x)
}

#' @export
predict.cluster_mm <- function(object,
                               newdata = NULL,
                               se.fit = FALSE,
                               interval = c("none", "confidence", "prediction"),
                               level = object$interval_level,
                               ...) {
  interval <- match.arg(interval)
  x_new <- if (is.null(newdata)) {
    object$x
  } else {
    suppressWarnings(as.numeric(newdata))
  }

  if (any(!is.finite(x_new))) {
    stop("newdata must be numeric.")
  }

  beta <- coef(object)
  fit <- mm_mean(x_new, beta["Vmax"], beta["Km"])

  if (!isTRUE(se.fit) && identical(interval, "none")) {
    return(fit)
  }

  grad <- mm_grad(x_new, beta["Vmax"], beta["Km"])
  mean_var <- rowSums((grad %*% vcov(object)) * grad)
  se_mean <- sqrt(pmax(mean_var, 0))
  z_new <- x_new / (beta["Km"] + x_new)
  work_var <- object$scale * safe_pos(object$variance_fun(x_new))
  marginal_var <- object$tau2 * z_new^2 + work_var

  out <- data.frame(
    x = x_new,
    fit = fit,
    se.fit = se_mean,
    working_variance = work_var,
    marginal_variance = marginal_var
  )

  if (!identical(interval, "none")) {
    alpha <- 1 - level
    z_val <- stats::qnorm(1 - alpha / 2)
    band_se <- if (identical(interval, "confidence")) {
      se_mean
    } else {
      sqrt(se_mean^2 + marginal_var)
    }
    out$lwr <- fit - z_val * band_se
    out$upr <- fit + z_val * band_se
  }

  out
}

#' @export
plot.cluster_mm <- function(x,
                            interval = TRUE,
                            interval_type = c("prediction", "confidence"),
                            level = x$interval_level,
                            n_points = 200,
                            xlab = "Substrate concentration",
                            ylab = "Reaction velocity",
                            main = NULL,
                            pch = 19,
                            col.points = "#1f78b4",
                            col.line = "#222222",
                            col.band = grDevices::adjustcolor("#9ecae1", alpha.f = 0.45),
                            ...) {
  dots <- list(...)
  dots[c("x", "y", "xlab", "ylab", "main", "pch", "col")] <- NULL
  interval_type <- match.arg(interval_type)
  x_seq <- seq(min(x$x), max(x$x), length.out = n_points)
  pred <- if (isTRUE(interval)) {
    predict(
      x,
      newdata = x_seq,
      interval = interval_type,
      level = level,
      se.fit = TRUE
    )
  } else {
    predict(
      x,
      newdata = x_seq,
      interval = "none",
      se.fit = FALSE
    )
  }

  fit_values <- if (isTRUE(interval)) pred$fit else pred
  x_range <- range(c(x$x, x_seq), finite = TRUE)
  y_range <- if (isTRUE(interval)) {
    range(c(x$y, pred$lwr, pred$upr), finite = TRUE)
  } else {
    range(c(x$y, fit_values), finite = TRUE)
  }
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- max(abs(y_range), 1)
  }
  y_pad <- 0.04 * y_span

  if (is.null(dots$xlim)) {
    dots$xlim <- x_range
  }
  if (is.null(dots$ylim)) {
    dots$ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  }

  do.call(
    graphics::plot,
    c(
      list(
        x = x$x,
        y = x$y,
        xlab = xlab,
        ylab = ylab,
        main = if (is.null(main)) paste("Clustered fit:", x$variance_label) else main,
        pch = pch,
        col = col.points
      ),
      dots
    )
  )

  if (isTRUE(interval)) {
    graphics::polygon(
      c(pred$x, rev(pred$x)),
      c(pred$lwr, rev(pred$upr)),
      border = NA,
      col = col.band
    )
    graphics::points(x$x, x$y, pch = pch, col = col.points)
  }

  graphics::lines(x_seq, fit_values, lwd = 2, col = col.line)
  invisible(x)
}
