prepare_grouped_mm_data <- function(data,
                                    s,
                                    v,
                                    groups = NULL,
                                    invalid_response_values = c(-9999),
                                    invalid_substrate_values = numeric(0),
                                    allow_zero_substrate = TRUE,
                                    sort_x = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  needed <- c(s, v, groups)
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

  if (length(groups)) {
    group_parts <- lapply(groups, function(col) as.character(dat[[col]]))
    dat$.group_id <- interaction(group_parts, drop = TRUE, lex.order = TRUE)
    dat$.group_label <- apply(
      as.data.frame(group_parts, stringsAsFactors = FALSE),
      1,
      function(z) paste(paste(groups, z, sep = "="), collapse = " | ")
    )
  } else {
    dat$.group_id <- factor("curve_1")
    dat$.group_label <- "curve_1"
  }

  if (isTRUE(sort_x)) {
    ord <- order(dat$.group_id, dat$.S, dat$.row_id)
    dat <- dat[ord, , drop = FALSE]
  }

  row.names(dat) <- NULL
  dat
}

make_group_comparison_tables <- function(estimates) {
  ok <- estimates[estimates$success == 1, , drop = FALSE]

  if (!nrow(ok)) {
    return(list(rank_by_group = ok, best_by_group = ok))
  }

  rank_by_group <- do.call(
    rbind,
    lapply(split(ok, ok$group_id), function(dat) {
      dat$rank_quasi_aic <- rank_metric(dat$quasi_aic)
      dat$rank_quasi_bic <- rank_metric(dat$quasi_bic)
      dat[order(dat$quasi_aic, dat$quasi_bic, dat$rmse), , drop = FALSE]
    })
  )
  row.names(rank_by_group) <- NULL

  best_by_group <- do.call(
    rbind,
    lapply(split(rank_by_group, rank_by_group$group_id), function(dat) {
      dat[1, , drop = FALSE]
    })
  )
  row.names(best_by_group) <- NULL

  list(rank_by_group = rank_by_group, best_by_group = best_by_group)
}

#' Fit Multiple Michaelis-Menten Curves by Group
#'
#' Clean a data frame, fit one or more working variance models within each
#' group, and compare the competing fits using quasi-AIC as the primary
#' criterion, with quasi-BIC and RMSE reported alongside it.
#'
#' @param data A data frame containing Michaelis-Menten observations.
#' @param s Name of the substrate concentration column.
#' @param v Name of the response column.
#' @param groups Optional character vector of grouping columns.
#' @param variance_models Character vector of built-in working variance models.
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
#' @param invalid_response_values Values to remove before fitting.
#' @param invalid_substrate_values Values to remove before fitting.
#' @param allow_zero_substrate Logical; if `FALSE`, substrate concentrations must
#'   be strictly positive.
#' @param sort_x Logical; sort observations within each group by concentration.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return An object of class `"group_mm"`.
#' @export
group_mm <- function(data,
                            s,
                            v,
                            groups = NULL,
                            variance_models = default_variance_models(),
                            power_values = NULL,
                            include_auto = FALSE,
                            power_selection = c("quasi_aic", "quasi_bic"),
                            power_grid = default_power_grid(),
                            interval_level = 0.95,
                            invalid_response_values = c(-9999),
                            invalid_substrate_values = numeric(0),
                            allow_zero_substrate = TRUE,
                            sort_x = TRUE,
                            quiet = FALSE) {
  power_selection <- match.arg(power_selection)
  cleaned <- prepare_grouped_mm_data(
    data = data,
    s = s,
    v = v,
    groups = groups,
    invalid_response_values = invalid_response_values,
    invalid_substrate_values = invalid_substrate_values,
    allow_zero_substrate = allow_zero_substrate,
    sort_x = sort_x
  )

  curve_split <- split(cleaned, cleaned$.group_id, drop = TRUE)
  estimate_rows <- list()
  augmented_rows <- list()
  fit_objects <- list()
  idx <- 1L

  for (group_name in names(curve_split)) {
    curve_data <- curve_split[[group_name]]

    if (!isTRUE(quiet)) {
      message("Fitting group: ", curve_data$.group_label[1])
    }

    screen <- screen_mm(
      x = curve_data$.S,
      y = curve_data$.V,
      variance_models = variance_models,
      power_values = power_values,
      include_auto = include_auto,
      power_selection = power_selection,
      power_grid = power_grid,
      interval_level = interval_level,
      allow_zero_substrate = allow_zero_substrate,
      quiet = TRUE
    )

    screen_table <- screen$table
    screen_table$group_id <- as.character(curve_data$.group_id[1])
    screen_table$group_label <- curve_data$.group_label[1]
    screen_table$n_obs <- nrow(curve_data)

    screen_table$se_Vmax <- NA_real_
    screen_table$se_Km <- NA_real_
    screen_table$lower_Vmax <- NA_real_
    screen_table$upper_Vmax <- NA_real_
    screen_table$lower_Km <- NA_real_
    screen_table$upper_Km <- NA_real_
    screen_table$quasi_loglik <- NA_real_

    if (length(screen$fits)) {
      for (model_name in names(screen$fits)) {
        fit <- screen$fits[[model_name]]
        hit <- which(screen_table$model == model_name)
        if (!length(hit)) {
          next
        }

        ci <- confint(fit, level = interval_level)
        se <- sqrt(pmax(diag(vcov(fit)), 0))
        screen_table$se_Vmax[hit] <- se["Vmax"]
        screen_table$se_Km[hit] <- se["Km"]
        screen_table$lower_Vmax[hit] <- ci["Vmax", "lower"]
        screen_table$upper_Vmax[hit] <- ci["Vmax", "upper"]
        screen_table$lower_Km[hit] <- ci["Km", "lower"]
        screen_table$upper_Km[hit] <- ci["Km", "upper"]
        screen_table$quasi_loglik[hit] <- fit$quasi_loglik

        augmented_rows[[length(augmented_rows) + 1L]] <- data.frame(
          curve_data,
          model = model_name,
          model_label = fit$variance_label,
          fitted = fit$fitted.values,
          residual = fit$residuals,
          variance_hat = fit$variance_hat,
          stringsAsFactors = FALSE
        )
      }
    }

    estimate_rows[[idx]] <- screen_table[, c(
      "group_id", "group_label", "model", "model_label",
      "selected_model", "selected_model_label", "n_obs", "success",
      "method", "runtime", "Vmax_hat", "Km_hat", "scale_hat", "se_Vmax",
      "se_Km", "lower_Vmax", "upper_Vmax", "lower_Km", "upper_Km",
      "rmse", "weighted_rss", "quasi_loglik", "quasi_aic", "quasi_bic",
      "error_message"
    )]
    fit_objects[[group_name]] <- screen$fits
    idx <- idx + 1L
  }

  estimates <- do.call(rbind, estimate_rows)
  row.names(estimates) <- NULL

  augmented <- if (length(augmented_rows)) {
    out <- do.call(rbind, augmented_rows)
    row.names(out) <- NULL
    out
  } else {
    NULL
  }

  comparison <- make_group_comparison_tables(estimates)
  group_reference <- unique(cleaned[, c(".group_id", ".group_label", groups), drop = FALSE])
  names(group_reference)[names(group_reference) == ".group_id"] <- "group_id"
  names(group_reference)[names(group_reference) == ".group_label"] <- "group_label"

  structure(
    list(
      cleaned_data = cleaned,
      estimates = estimates,
      augmented = augmented,
      comparison = comparison,
      group_reference = group_reference,
      fits = fit_objects,
      call = match.call()
    ),
    class = "group_mm"
  )
}

#' @export
print.group_mm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Call:\n")
  print(x$call)
  cat(
    "\nGrouped Michaelis-Menten analysis for ",
    nrow(x$group_reference), " group(s).\n\n",
    sep = ""
  )

  best <- x$comparison$best_by_group
  if (nrow(best)) {
    show <- best[, c(
      "group_id", "group_label", "model", "selected_model",
      "quasi_aic", "quasi_bic", "rmse"
    )]
    numeric_cols <- vapply(show, is.numeric, logical(1))
    show[numeric_cols] <- lapply(show[numeric_cols], round, digits = digits)
    cat("Best model by group (ordered by quasi-AIC):\n")
    print(show, row.names = FALSE)
  } else {
    cat("No successful fits were produced.\n")
  }

  invisible(x)
}

#' @export
plot.group_mm <- function(x,
                          which_groups = NULL,
                          model = NULL,
                          interval = TRUE,
                          interval_type = c("prediction", "confidence"),
                          level = NULL,
                          ncol = NULL,
                          main = NULL,
                          xlab = "Substrate concentration",
                          ylab = "Reaction velocity",
                          ...) {
  interval_type <- match.arg(interval_type)

  reference <- x$group_reference
  group_ids <- as.character(reference$group_id)
  group_labels <- as.character(reference$group_label)

  if (is.null(which_groups)) {
    selected_ids <- group_ids
  } else {
    requested <- as.character(which_groups)
    matched <- group_ids %in% requested | group_labels %in% requested
    selected_ids <- unique(group_ids[matched])
    if (!length(selected_ids)) {
      stop("No matching groups were found in the grouped fit object.")
    }
  }

  if (!is.null(model) && length(model) != 1L) {
    stop("model must be NULL or a single model name.")
  }

  best_table <- x$comparison$best_by_group
  plot_specs <- vector("list", length(selected_ids))

  for (i in seq_along(selected_ids)) {
    group_id <- selected_ids[i]
    label_hit <- match(group_id, group_ids)
    group_label <- group_labels[label_hit]
    model_name <- if (is.null(model)) {
      best_table$model[match(group_id, best_table$group_id)]
    } else {
      model
    }

    if (!length(model_name) || is.na(model_name)) {
      stop("No selected model is available for group '", group_label, "'.")
    }

    fit <- x$fits[[group_id]][[model_name]]
    if (is.null(fit)) {
      stop(
        "No fitted object is stored for group '", group_label,
        "' and model '", model_name, "'."
      )
    }

    plot_specs[[i]] <- list(
      fit = fit,
      group_label = group_label,
      model = model_name,
      selected_model = best_table$selected_model[match(group_id, best_table$group_id)],
      selected_model_label = best_table$selected_model_label[match(group_id, best_table$group_id)]
    )
  }

  n_panels <- length(plot_specs)
  if (!n_panels) {
    stop("No fitted groups are available to plot.")
  }

  if (is.null(ncol)) {
    ncol <- min(3L, max(1L, ceiling(sqrt(n_panels))))
  }
  nrow <- ceiling(n_panels / ncol)

  main_values <- if (is.null(main)) {
    rep(NA_character_, n_panels)
  } else {
    rep_len(as.character(main), n_panels)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(nrow, ncol))

  for (i in seq_along(plot_specs)) {
    spec <- plot_specs[[i]]
    panel_main <- if (is.na(main_values[i])) {
      paste0(
        spec$group_label, "\n",
        if (!is.na(spec$selected_model_label) && nzchar(spec$selected_model_label)) {
          spec$selected_model_label
        } else if (!is.na(spec$selected_model) && nzchar(spec$selected_model)) {
          spec$selected_model
        } else {
          spec$model
        }
      )
    } else {
      main_values[i]
    }

    plot(
      spec$fit,
      interval = interval,
      interval_type = interval_type,
      level = if (is.null(level)) spec$fit$interval_level else level,
      xlab = xlab,
      ylab = ylab,
      main = panel_main,
      ...
    )
  }

  invisible(x)
}
