#' @export
coef.fit_mm <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.fit_mm <- function(object, ...) {
  object$vcov
}

#' @export
confint.fit_mm <- function(object,
                           parm = c("Vmax", "Km"),
                           level = object$interval_level,
                           method = c("wald", "bootstrap"),
                           B = 399L,
                           bootstrap_type = c("wild", "pairs"),
                           bootstrap_ci = c("studentized", "basic", "percentile"),
                           wild_weights = c("mammen", "rademacher"),
                           seed = NULL,
                           ...) {
  method <- match.arg(method)
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  wild_weights <- match.arg(wild_weights)
  policy <- small_sample_ci_policy(object)

  if (isTRUE(policy$disable)) {
    warning(
      paste(
        policy$confint_warning,
        "Use these intervals only as a sensitivity analysis."
      ),
      call. = FALSE
    )
  }

  interval <- if (identical(method, "wald")) {
    wald_interval(coef(object), vcov(object), level = level)
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

  out <- cbind(interval$lower, interval$upper)
  colnames(out) <- c("lower", "upper")
  out <- out[parm, , drop = FALSE]
  if (identical(method, "bootstrap")) {
    attr(out, "bootstrap_draws") <- interval$draws
  }
  out
}

#' @export
print.fit_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  cat("Working variance:", x$variance_label, "\n\n")
  policy <- small_sample_ci_policy(x)

  coef_table <- if (isTRUE(policy$disable)) {
    estimate_only_table(x)
  } else {
    coefficient_table(x, level = x$interval_level)
  }
  print(round(coef_table, digits = digits))

  if (isTRUE(policy$disable)) {
    cat(
      "\n",
      policy$summary_message,
      "\nUse confint() manually for a sensitivity analysis.\n",
      sep = ""
    )
  }

  cat(
    "\nQuasi-AIC:", round(x$quasi_aic, digits),
    "  Quasi-BIC:", round(x$quasi_bic, digits),
    "  RMSE:", round(x$rmse, digits), "\n",
    sep = ""
  )

  invisible(x)
}

#' @export
summary.fit_mm <- function(object,
                           level = object$interval_level,
                           method = c("wald", "bootstrap"),
                           B = 399L,
                           bootstrap_type = c("wild", "pairs"),
                           bootstrap_ci = c("studentized", "basic", "percentile"),
                           wild_weights = c("mammen", "rademacher"),
                           seed = NULL,
                           ...) {
  method <- match.arg(method)
  bootstrap_type <- match.arg(bootstrap_type)
  bootstrap_ci <- match.arg(bootstrap_ci)
  wild_weights <- match.arg(wild_weights)
  policy <- small_sample_ci_policy(object)

  out <- list(
    call = object$call,
    method = object$method,
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
        wild_weights = wild_weights,
        seed = seed
      )
    },
    scale = object$scale,
    rss = object$rss,
    rmse = object$rmse,
    weighted_rss = object$weighted_rss,
    quasi_aic = object$quasi_aic,
    quasi_bic = object$quasi_bic,
    level = level,
    ci_method = if (isTRUE(policy$disable)) "disabled" else method,
    bootstrap_type = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) bootstrap_type else NULL,
    bootstrap_ci = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) bootstrap_ci else NULL,
    wild_weights = if (!isTRUE(policy$disable) &&
      identical(method, "bootstrap") &&
      identical(bootstrap_type, "wild")) wild_weights else NULL,
    bootstrap_B = if (!isTRUE(policy$disable) && identical(method, "bootstrap")) B else NULL,
    small_sample_policy = policy
  )
  class(out) <- "summary.fit_mm"
  out
}

#' @export
print.summary.fit_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method, "\n")
  cat("Working variance:", x$variance_label, "\n\n")
  if (identical(x$ci_method, "disabled")) {
    cat("CI method: disabled by default for small samples\n\n")
  } else if (identical(x$ci_method, "bootstrap")) {
    if (identical(x$bootstrap_type, "wild") && !is.null(x$wild_weights)) {
      cat(
        "CI method: bootstrap (",
        x$bootstrap_type,
        ", ",
        x$wild_weights,
        " weights, ",
        x$bootstrap_ci,
        ", B = ",
        x$bootstrap_B,
        ")\n\n",
        sep = ""
      )
    } else {
      cat(
        "CI method: bootstrap (",
        x$bootstrap_type,
        ", ",
        x$bootstrap_ci,
        ", B = ",
        x$bootstrap_B,
        ")\n\n",
        sep = ""
      )
    }
  } else {
    cat("CI method: wald\n\n")
  }
  print(round(x$coefficients, digits = digits))
  if (identical(x$ci_method, "disabled") && !is.null(x$small_sample_policy)) {
    cat(
      "\n",
      x$small_sample_policy$summary_message,
      "\nUse confint() manually for a sensitivity analysis.\n",
      sep = ""
    )
  }
  cat(
    "\nScale:", round(x$scale, digits),
    "  RMSE:", round(x$rmse, digits),
    "  Weighted RSS:", round(x$weighted_rss, digits), "\n",
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
predict.fit_mm <- function(object,
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
  work_var <- object$scale * safe_pos(object$variance_fun(x_new))

  out <- data.frame(
    x = x_new,
    fit = fit,
    se.fit = se_mean,
    working_variance = work_var
  )

  if (!identical(interval, "none")) {
    alpha <- 1 - level
    z_val <- stats::qnorm(1 - alpha / 2)
    band_se <- if (identical(interval, "confidence")) {
      se_mean
    } else {
      sqrt(se_mean^2 + work_var)
    }
    out$lwr <- fit - z_val * band_se
    out$upr <- fit + z_val * band_se
  }

  out
}

#' @export
plot.fit_mm <- function(x,
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
        main = if (is.null(main)) x$variance_label else main,
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
