safe_pos <- function(x, lower = 1e-8, upper = 1e10) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- upper
  pmin(pmax(x, lower), upper)
}

ignore_seed_argument <- function(seed) {
  invisible(seed)
}

safe_inverse <- function(mat, ridge = 1e-8, max_tries = 8L) {
  p <- nrow(mat)
  eye <- diag(p)

  for (j in 0:max_tries) {
    out <- tryCatch(
      solve(mat + ridge * 10^j * eye),
      error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) {
      return(out)
    }
  }

  stop("Unable to invert the information matrix.")
}

mm_mean <- function(x, vmax, km) {
  vmax * x / (km + x)
}

mm_grad <- function(x, vmax, km) {
  cbind(
    Vmax = x / (km + x),
    Km = -(vmax * x) / (km + x)^2
  )
}

default_km_bounds <- function(x) {
  upper <- max(x, na.rm = TRUE) * 10 + 1
  c(1e-8, max(upper, 1))
}

validate_mm_inputs <- function(x, y, allow_zero_substrate = TRUE) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }

  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y)

  if (!all(keep)) {
    warning("Dropping rows with non-finite x or y values.", call. = FALSE)
    x <- x[keep]
    y <- y[keep]
  }

  if (!length(x)) {
    stop("No valid observations remain after removing missing values.")
  }

  if (isTRUE(allow_zero_substrate)) {
    if (any(x < 0)) {
      stop("Substrate concentrations must be non-negative.")
    }
  } else if (any(x <= 0)) {
    stop("Substrate concentrations must be strictly positive.")
  }

  if (length(unique(x)) < 3L) {
    stop("At least 3 distinct substrate concentrations are required.")
  }

  list(x = x, y = y)
}

wald_interval <- function(beta, vcov_beta, level = 0.95) {
  alpha <- 1 - level
  z_val <- stats::qnorm(1 - alpha / 2)
  se <- sqrt(pmax(diag(vcov_beta), 0))

  list(
    lower = beta - z_val * se,
    upper = beta + z_val * se,
    se = se
  )
}

refit_fit_mm <- function(object, x = object$x, y = object$y) {
  if (!inherits(object, "fit_mm")) {
    stop("object must inherit from 'fit_mm'.")
  }

  if (identical(object$variance_selection_mode, "auto")) {
    return(
      fit_mm(
        x = x,
        y = y,
        variance = "auto",
        power_selection = object$variance_selection,
        power_grid = object$power_grid,
        interval_level = object$interval_level,
        km_bounds = object$km_bounds,
        allow_zero_substrate = object$allow_zero_substrate
      )
    )
  }

  fit_mm(
    x = x,
    y = y,
    variance = object$variance_model,
    power = object$variance_power,
    interval_level = object$interval_level,
    km_bounds = object$km_bounds,
    allow_zero_substrate = object$allow_zero_substrate
  )
}

resample_cluster_pairs_data <- function(cleaned) {
  cluster_split <- split(cleaned, cleaned$.cluster_id, drop = TRUE)
  sampled <- sample.int(length(cluster_split), size = length(cluster_split), replace = TRUE)
  out <- vector("list", length(sampled))

  for (j in seq_along(sampled)) {
    dat_j <- cluster_split[[sampled[[j]]]]
    source_label <- unique(dat_j$.cluster_label)
    dat_j$.cluster_id <- sprintf("boot_%03d", j)
    dat_j$.cluster_label <- if (length(source_label)) {
      sprintf("bootstrap[%d] <- %s", j, source_label[[1]])
    } else {
      sprintf("bootstrap[%d]", j)
    }
    dat_j$.row_id <- seq_len(nrow(dat_j))
    out[[j]] <- dat_j
  }

  boot <- do.call(rbind, out)
  row.names(boot) <- NULL
  boot
}

draw_wild_multipliers <- function(n, weights = c("mammen", "rademacher")) {
  weights <- match.arg(weights)

  if (identical(weights, "rademacher")) {
    return(sample(c(-1, 1), size = n, replace = TRUE))
  }

  sqrt5 <- sqrt(5)
  values <- c((1 - sqrt5) / 2, (1 + sqrt5) / 2)
  probs <- c((sqrt5 + 1) / (2 * sqrt5), (sqrt5 - 1) / (2 * sqrt5))
  sample(values, size = n, replace = TRUE, prob = probs)
}

bootstrap_fit_mm <- function(object,
                             B = 399L,
                             bootstrap_type = c("wild", "pairs"),
                             wild_weights = c("mammen", "rademacher"),
                             seed = NULL) {
  bootstrap_type <- match.arg(bootstrap_type)
  wild_weights <- match.arg(wild_weights)
  B <- as.integer(B)

  if (length(B) != 1L || !is.finite(B) || B < 20L) {
    stop("B must be a single integer of at least 20.")
  }

  ignore_seed_argument(seed)

  beta_names <- names(coef(object))
  draws <- matrix(NA_real_, nrow = B, ncol = length(beta_names))
  se_draws <- matrix(NA_real_, nrow = B, ncol = length(beta_names))
  colnames(draws) <- beta_names
  colnames(se_draws) <- beta_names

  x <- object$x
  y <- object$y
  fitted <- object$fitted.values
  residuals <- object$residuals
  n <- length(y)

  for (b in seq_len(B)) {
    if (identical(bootstrap_type, "wild")) {
      multipliers <- draw_wild_multipliers(n, weights = wild_weights)
      x_star <- x
      y_star <- fitted + residuals * multipliers
    } else {
      idx <- sample.int(n, size = n, replace = TRUE)
      x_star <- x[idx]
      y_star <- y[idx]
    }

    fit_star <- suppressWarnings(
      tryCatch(
        refit_fit_mm(object, x = x_star, y = y_star),
        error = function(e) NULL
      )
    )

    if (inherits(fit_star, "fit_mm")) {
      draws[b, ] <- coef(fit_star)[beta_names]
      se_draws[b, ] <- sqrt(pmax(diag(vcov(fit_star))[beta_names], 0))
    }
  }

  keep <- stats::complete.cases(draws)
  draws <- draws[keep, , drop = FALSE]

  if (nrow(draws) < max(20L, ceiling(B / 4))) {
    stop("Too few successful bootstrap refits. Consider using method = 'wald'.")
  }

  if (nrow(draws) < B) {
    warning(
      sprintf(
        "Only %d of %d bootstrap refits succeeded.",
        nrow(draws),
        B
      ),
      call. = FALSE
    )
  }

  attr(draws, "requested_B") <- B
  attr(draws, "successful_B") <- nrow(draws)
  attr(draws, "bootstrap_type") <- bootstrap_type
  attr(draws, "wild_weights") <- if (identical(bootstrap_type, "wild")) wild_weights else NULL
  attr(draws, "se_draws") <- se_draws[keep, , drop = FALSE]
  draws
}

bootstrap_fit_cluster_mm <- function(object,
                                     B = 399L,
                                     bootstrap_type = c("pairs"),
                                     seed = NULL) {
  bootstrap_type <- match.arg(bootstrap_type)
  B <- as.integer(B)

  if (length(B) != 1L || !is.finite(B) || B < 20L) {
    stop("B must be a single integer of at least 20.")
  }

  ignore_seed_argument(seed)

  beta_names <- names(coef(object))
  draws <- matrix(NA_real_, nrow = B, ncol = length(beta_names))
  se_draws <- matrix(NA_real_, nrow = B, ncol = length(beta_names))
  colnames(draws) <- beta_names
  colnames(se_draws) <- beta_names

  for (b in seq_len(B)) {
    cleaned_star <- resample_cluster_pairs_data(object$cleaned_data)
    fit_star <- suppressWarnings(
      tryCatch(
        refit_cluster_mm(object, cleaned = cleaned_star),
        error = function(e) NULL
      )
    )

    if (inherits(fit_star, "cluster_mm")) {
      draws[b, ] <- coef(fit_star)[beta_names]
      se_draws[b, ] <- sqrt(pmax(diag(vcov(fit_star))[beta_names], 0))
    }
  }

  keep <- stats::complete.cases(draws)
  draws <- draws[keep, , drop = FALSE]

  if (nrow(draws) < max(20L, ceiling(B / 4))) {
    stop("Too few successful clustered bootstrap refits. Consider using method = 'wald'.")
  }

  if (nrow(draws) < B) {
    warning(
      sprintf(
        "Only %d of %d clustered bootstrap refits succeeded.",
        nrow(draws),
        B
      ),
      call. = FALSE
    )
  }

  attr(draws, "requested_B") <- B
  attr(draws, "successful_B") <- nrow(draws)
  attr(draws, "bootstrap_type") <- bootstrap_type
  attr(draws, "se_draws") <- se_draws[keep, , drop = FALSE]
  draws
}

bootstrap_interval <- function(object,
                               level = 0.95,
                               B = 399L,
                               bootstrap_type = c("wild", "pairs"),
                               bootstrap_ci = c("studentized", "basic", "percentile"),
                               wild_weights = c("mammen", "rademacher"),
                               seed = NULL) {
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  wild_weights <- match.arg(wild_weights)
  draws <- bootstrap_fit_mm(
    object = object,
    B = B,
    bootstrap_type = bootstrap_type,
    wild_weights = wild_weights,
    seed = seed
  )
  se_draws <- attr(draws, "se_draws")
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  theta_hat <- coef(object)[colnames(draws)]
  se_hat <- sqrt(pmax(diag(vcov(object))[colnames(draws)], 0))

  if (identical(bootstrap_ci, "studentized")) {
    ok_student <- stats::complete.cases(draws, se_draws) &
      apply(se_draws > 0, 1, all)

    if (sum(ok_student) < max(20L, ceiling(nrow(draws) / 4))) {
      warning(
        "Too few finite studentized bootstrap refits; falling back to the basic interval.",
        call. = FALSE
      )
      bootstrap_ci <- "basic"
    } else {
      t_star <- sweep(draws[ok_student, , drop = FALSE], 2, theta_hat, "-")
      t_star <- t_star / se_draws[ok_student, , drop = FALSE]
      t_quantiles <- apply(
        t_star,
        2,
        stats::quantile,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )

      if (is.null(dim(t_quantiles))) {
        t_quantiles <- matrix(t_quantiles, nrow = 2L)
      }

      lower <- theta_hat - t_quantiles[2, ] * se_hat
      upper <- theta_hat - t_quantiles[1, ] * se_hat
    }
  }

  if (!identical(bootstrap_ci, "studentized")) {
    quantiles <- apply(
      draws,
      2,
      stats::quantile,
      probs = probs,
      na.rm = TRUE,
      names = FALSE
    )

    if (is.null(dim(quantiles))) {
      quantiles <- matrix(quantiles, nrow = 2L)
    }

    if (identical(bootstrap_ci, "basic")) {
      lower <- 2 * theta_hat - quantiles[2, ]
      upper <- 2 * theta_hat - quantiles[1, ]
    } else {
      lower <- quantiles[1, ]
      upper <- quantiles[2, ]
    }
  }

  attr(draws, "bootstrap_ci") <- bootstrap_ci

  list(
    lower = lower,
    upper = upper,
    se = apply(draws, 2, stats::sd),
    draws = draws
  )
}

cluster_bootstrap_interval <- function(object,
                                       level = 0.95,
                                       B = 399L,
                                       bootstrap_type = c("pairs"),
                                       bootstrap_ci = c("studentized", "basic", "percentile"),
                                       seed = NULL) {
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  draws <- bootstrap_fit_cluster_mm(
    object = object,
    B = B,
    bootstrap_type = bootstrap_type,
    seed = seed
  )
  se_draws <- attr(draws, "se_draws")
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  theta_hat <- coef(object)[colnames(draws)]
  se_hat <- sqrt(pmax(diag(vcov(object))[colnames(draws)], 0))

  if (identical(bootstrap_ci, "studentized")) {
    ok_student <- stats::complete.cases(draws, se_draws) &
      apply(se_draws > 0, 1, all)

    if (sum(ok_student) < max(20L, ceiling(nrow(draws) / 4))) {
      warning(
        "Too few finite studentized clustered bootstrap refits; falling back to the basic interval.",
        call. = FALSE
      )
      bootstrap_ci <- "basic"
    } else {
      t_star <- sweep(draws[ok_student, , drop = FALSE], 2, theta_hat, "-")
      t_star <- t_star / se_draws[ok_student, , drop = FALSE]
      t_quantiles <- apply(
        t_star,
        2,
        stats::quantile,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )

      if (is.null(dim(t_quantiles))) {
        t_quantiles <- matrix(t_quantiles, nrow = 2L)
      }

      lower <- theta_hat - t_quantiles[2, ] * se_hat
      upper <- theta_hat - t_quantiles[1, ] * se_hat
    }
  }

  if (!identical(bootstrap_ci, "studentized")) {
    quantiles <- apply(
      draws,
      2,
      stats::quantile,
      probs = probs,
      na.rm = TRUE,
      names = FALSE
    )

    if (is.null(dim(quantiles))) {
      quantiles <- matrix(quantiles, nrow = 2L)
    }

    if (identical(bootstrap_ci, "basic")) {
      lower <- 2 * theta_hat - quantiles[2, ]
      upper <- 2 * theta_hat - quantiles[1, ]
    } else {
      lower <- quantiles[1, ]
      upper <- quantiles[2, ]
    }
  }

  attr(draws, "bootstrap_ci") <- bootstrap_ci

  list(
    lower = lower,
    upper = upper,
    se = apply(draws, 2, stats::sd),
    draws = draws
  )
}

quasi_metrics <- function(y, fitted, variance_hat, k = 3L) {
  residuals <- y - fitted
  variance_hat <- safe_pos(variance_hat)
  quasi_loglik <- -0.5 * sum(log(2 * pi * variance_hat) + residuals^2 / variance_hat)
  n <- length(y)

  list(
    rss = sum(residuals^2),
    rmse = sqrt(mean(residuals^2)),
    weighted_rss = sum(residuals^2 / variance_hat),
    quasi_loglik = quasi_loglik,
    quasi_aic = -2 * quasi_loglik + 2 * k,
    quasi_bic = -2 * quasi_loglik + log(n) * k
  )
}

rank_metric <- function(x) {
  if (!length(x) || all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  rank(x, ties.method = "average", na.last = "keep")
}

small_sample_ci_policy <- function(object,
                                   min_n = 20L,
                                   min_unique_x = 6L) {
  n_obs <- length(object$y)
  n_unique_x <- length(unique(object$x))
  disable <- n_obs < min_n || n_unique_x < min_unique_x

  list(
    disable = disable,
    n_obs = n_obs,
    n_unique_x = n_unique_x,
    min_n = as.integer(min_n),
    min_unique_x = as.integer(min_unique_x),
    summary_message = sprintf(
      paste(
        "Interval inference is disabled by default for small samples",
        "(n = %d, %d distinct substrate concentrations).",
        "Default intervals require at least n >= %d and at least %d",
        "distinct concentrations."
      ),
      n_obs,
      n_unique_x,
      min_n,
      min_unique_x
    ),
    confint_warning = sprintf(
      paste(
        "Small-sample warning: this fit uses n = %d observations across %d",
        "distinct substrate concentrations.",
        "Default interval reporting is disabled below n = %d or below %d",
        "distinct concentrations because coverage may be unstable."
      ),
      n_obs,
      n_unique_x,
      min_n,
      min_unique_x
    )
  )
}

estimate_only_table <- function(object) {
  out <- cbind(Estimate = coef(object))
  out[rownames(out), , drop = FALSE]
}

coefficient_table <- function(object,
                              level = object$interval_level,
                              method = c("wald", "bootstrap"),
                              B = 399L,
                              bootstrap_type = c("wild", "pairs"),
                              bootstrap_ci = c("studentized", "basic", "percentile"),
                              wild_weights = c("mammen", "rademacher"),
                              seed = NULL) {
  method <- match.arg(method)
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  wild_weights <- match.arg(wild_weights)

  interval <- if (identical(method, "wald")) {
    wald_interval(coef(object), vcov(object), level = level)
  } else if (inherits(object, "cluster_mm")) {
    cluster_bootstrap_interval(
      object = object,
      level = level,
      B = B,
      bootstrap_type = bootstrap_type,
      bootstrap_ci = bootstrap_ci,
      seed = seed
    )
  } else {
    bootstrap_interval(
      object = object,
      level = level,
      B = B,
      bootstrap_type = bootstrap_type,
      bootstrap_ci = bootstrap_ci,
      wild_weights = wild_weights,
      seed = seed
    )
  }

  se <- interval$se

  out <- cbind(
    Estimate = coef(object),
    `Std. Error` = se[names(coef(object))],
    `Lower CI` = interval$lower[names(coef(object))],
    `Upper CI` = interval$upper[names(coef(object))]
  )

  out[rownames(out), , drop = FALSE]
}

model_display_label <- function(model) {
  switch(
    model,
    constant = "NLS",
    log = "Var = gamma log(S + 1)",
    sqrt = "Var = gamma S^(1/2)",
    cuberoot = "Var = gamma S^(1/3)",
    power = "Var = gamma S^a",
    auto = "Auto: best of log/power",
    model
  )
}
