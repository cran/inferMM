library(inferMM)

simulate_mm_data <- function(n = 50L,
                             variance_shape = c("mm", "exp", "hill"),
                             error = c("normal", "skewed")) {
  variance_shape <- match.arg(variance_shape)
  error <- match.arg(error)

  x <- seq(1, 100, length.out = n)
  mu <- 100 * x / (20 + x)

  variance_fun <- switch(
    variance_shape,
    mm = function(s) 1 + 9 * s / (20 + s),
    exp = function(s) 1 + 9 * (1 - exp(-0.05 * s)),
    hill = function(s) 1 + 9 * s^2 / (20^2 + s^2)
  )
  z <- if (identical(error, "normal")) {
    stats::rnorm(n)
  } else {
    (stats::rgamma(n, shape = 4, scale = 1) - 4) / sqrt(4)
  }

  list(
    x = x,
    y = mu + sqrt(variance_fun(x)) * z
  )
}

simulate_clustered_mm_data <- function(n_clusters = 6L,
                                       reps = 4L,
                                       error = c("normal", "skewed")) {
  error <- match.arg(error)

  x_levels <- c(10, 25, 50, 100, 200, 400, 700)
  cluster <- rep(seq_len(n_clusters), each = length(x_levels) * reps)
  x <- rep(rep(x_levels, each = reps), times = n_clusters)
  vmax_cluster <- 1400 + stats::rnorm(n_clusters, sd = 120)
  km_true <- 35
  mu <- vmax_cluster[cluster] * x / (km_true + x)
  z <- if (identical(error, "normal")) {
    stats::rnorm(length(x))
  } else {
    (stats::rgamma(length(x), shape = 4, scale = 1) - 4) / sqrt(4)
  }

  data.frame(
    cluster = paste0("C", cluster),
    substrate_conc = x,
    activity = mu + sqrt(40 * sqrt(x)) * z,
    stringsAsFactors = FALSE
  )
}

test_that("fit_mm fits a variance-aware model", {
  set.seed(1)
  sim_dat <- simulate_mm_data(variance_shape = "mm", error = "normal")
  fit <- fit_mm(sim_dat$x, sim_dat$y, variance = "sqrt")

  expect_s3_class(fit, "fit_mm")
  expect_true(all(c("Vmax", "Km") %in% names(coef(fit))))
  expect_true(all(is.finite(coef(fit))))
  expect_equal(dim(confint(fit)), c(2, 2))
})

test_that("fit_mm uses the profile-score path for constant variance", {
  set.seed(3)
  sim_dat <- simulate_mm_data(variance_shape = "mm", error = "normal")
  fit <- fit_mm(sim_dat$x, sim_dat$y, variance = "constant")

  expect_s3_class(fit, "fit_mm")
  expect_identical(fit$variance_model, "constant")
  expect_identical(fit$method, "profile-score")
  expect_true(all(is.finite(coef(fit))))
})

test_that("fit_mm supports fixed-power and automatic variance choices", {
  set.seed(5)
  sim_dat <- simulate_mm_data(variance_shape = "exp", error = "normal")

  fit_power <- fit_mm(sim_dat$x, sim_dat$y, variance = "power", power = 0.4)
  expect_s3_class(fit_power, "fit_mm")
  expect_equal(fit_power$variance_model, "power")
  expect_equal(unname(fit_power$variance_power), 0.4)

  fit_auto <- fit_mm(
    sim_dat$x,
    sim_dat$y,
    variance = "auto",
    power_selection = "quasi_aic"
  )
  expect_s3_class(fit_auto, "fit_mm")
  expect_identical(fit_auto$variance_selection_mode, "auto")
  expect_true(is.data.frame(fit_auto$variance_candidate_table))
  expect_true(nrow(fit_auto$variance_candidate_table) >= 2L)
})

test_that("bootstrap confidence intervals are available", {
  set.seed(7)
  sim_dat <- simulate_mm_data(variance_shape = "exp", error = "skewed")
  fit <- fit_mm(sim_dat$x, sim_dat$y, variance = "sqrt")

  set.seed(1)
  ci_boot <- confint(
    fit,
    method = "bootstrap",
    B = 25,
    bootstrap_ci = "studentized",
    wild_weights = "mammen"
  )
  expect_equal(dim(ci_boot), c(2, 2))
  expect_false(is.null(attr(ci_boot, "bootstrap_draws")))
  expect_identical(attr(attr(ci_boot, "bootstrap_draws"), "bootstrap_ci"), "studentized")

  set.seed(1)
  summ_boot <- summary(
    fit,
    method = "bootstrap",
    B = 25,
    bootstrap_ci = "studentized",
    wild_weights = "mammen"
  )
  expect_identical(summ_boot$ci_method, "bootstrap")
  expect_identical(summ_boot$bootstrap_ci, "studentized")
  expect_identical(summ_boot$wild_weights, "mammen")
  expect_true(all(c("Estimate", "Std. Error", "Lower CI", "Upper CI") %in% colnames(summ_boot$coefficients)))
})

test_that("small-sample fits suppress default interval reporting", {
  set.seed(13)
  sim_dat <- simulate_mm_data(n = 10, variance_shape = "mm", error = "normal")
  fit <- fit_mm(sim_dat$x, sim_dat$y, variance = "sqrt")

  summ_small <- summary(fit)
  expect_identical(summ_small$ci_method, "disabled")
  expect_identical(colnames(summ_small$coefficients), "Estimate")
  expect_true(isTRUE(summ_small$small_sample_policy$disable))

  expect_warning(
    ci_small <- confint(fit),
    "Small-sample warning"
  )
  expect_equal(dim(ci_small), c(2, 2))

  txt_small <- capture.output(print(summ_small))
  expect_true(any(grepl("disabled by default for small samples", txt_small, fixed = TRUE)))
})

test_that("report_mm prints and plots without error", {
  set.seed(11)
  sim_dat <- simulate_mm_data(variance_shape = "mm", error = "normal")
  fit <- fit_mm(sim_dat$x, sim_dat$y, variance = "sqrt")
  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- report_mm(fit, interval_type = "confidence")
  expect_s3_class(out, "fit_mm")

  txt_boot <- capture.output({
    set.seed(1)
    report_mm(
      fit,
      method = "bootstrap",
      B = 25,
      bootstrap_ci = "studentized",
      wild_weights = "mammen",
      interval_type = "none"
    )
  })
  expect_true(any(grepl("CI method: bootstrap", txt_boot, fixed = TRUE)))

  set.seed(17)
  sim_small <- simulate_mm_data(n = 10, variance_shape = "exp", error = "normal")
  fit_small <- fit_mm(sim_small$x, sim_small$y, variance = "sqrt")
  txt_small <- capture.output(report_mm(fit_small, interval_type = "none"))
  expect_true(any(grepl("disabled by default for small samples", txt_small, fixed = TRUE)))
})

test_that("screen_mm returns a ranking table", {
  set.seed(2)
  sim_dat <- simulate_mm_data(variance_shape = "exp", error = "normal")
  out <- screen_mm(
    sim_dat$x,
    sim_dat$y,
    power_values = c(0.4, 0.6),
    include_auto = TRUE,
    quiet = TRUE
  )

  expect_s3_class(out, "screen_mm")
  expect_true(nrow(out$table) >= 2L)
  expect_true(all(c("selected_model", "selected_model_label", "quasi_aic", "quasi_bic", "rmse") %in% names(out$table)))
  ok <- out$table$success == 1
  expect_true(all(diff(out$table$quasi_aic[ok]) >= 0))
  expect_true(any(grepl("^power\\(", out$table$model)))
  expect_true(any(grepl("^auto\\(", out$table$model)))
  auto_hit <- grepl("^auto\\(", out$table$model)
  expect_true(all(!is.na(out$table$selected_model[auto_hit])))
  expect_false(any(grepl("^auto\\(", out$table$selected_model[auto_hit])))
})

test_that("group_mm fits multiple SDL curves", {
  out <- group_mm(
    data = sdl_demo,
    s = "s_uM",
    v = "v_uM_per_min",
    groups = "enzyme",
    variance_models = c("constant", "sqrt"),
    power_values = 0.4,
    include_auto = TRUE,
    quiet = TRUE
  )

  expect_s3_class(out, "group_mm")
  expect_true(nrow(out$comparison$best_by_group) >= 1L)
  expect_true(all(c("group_id", "model", "selected_model") %in% names(out$comparison$best_by_group)))
  auto_rows <- grepl("^auto\\(", out$estimates$model)
  if (any(auto_rows)) {
    expect_false(any(grepl("^auto\\(", out$estimates$selected_model[auto_rows])))
  }
})

test_that("plot.group_mm draws grouped panels without error", {
  out <- group_mm(
    data = sdl_demo,
    s = "s_uM",
    v = "v_uM_per_min",
    groups = "enzyme",
    variance_models = c("constant", "sqrt"),
    quiet = TRUE
  )

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted <- plot(out, which_groups = c("1111", "2222"), interval_type = "confidence")
  expect_s3_class(plotted, "group_mm")
})

test_that("plot.fit_mm expands limits to include interval bands", {
  one_curve <- subset(sdl_demo, enzyme == "4444")
  fit <- fit_mm(one_curve$s_uM, one_curve$v_uM_per_min, variance = "sqrt")
  pred <- predict(
    fit,
    newdata = seq(min(one_curve$s_uM), max(one_curve$s_uM), length.out = 200),
    interval = "prediction"
  )

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(fit, interval_type = "prediction")
  usr <- graphics::par("usr")

  expect_lte(usr[3], min(pred$lwr) + 1e-8)
  expect_gte(usr[4], max(pred$upr) - 1e-8)
})

test_that("cluster_mm fits repeated-measurement data and supports core methods", {
  set.seed(21)
  sim_dat <- simulate_clustered_mm_data(n_clusters = 10L, error = "normal")

  fit <- cluster_mm(
    data = sim_dat,
    s = "substrate_conc",
    v = "activity",
    cluster = "cluster",
    variance = "sqrt"
  )

  expect_s3_class(fit, "cluster_mm")
  expect_true(all(c("Vmax", "Km") %in% names(coef(fit))))
  expect_true(all(is.finite(coef(fit))))
  expect_identical(fit$initialization_method, "pooled single-curve profile-score fit")
  expect_equal(dim(vcov(fit)), c(2, 2))
  expect_equal(dim(confint(fit)), c(2, 2))

  summ <- summary(fit)
  expect_identical(summ$ci_method, "wald")
  expect_true(all(c("Estimate", "Std. Error", "Lower CI", "Upper CI") %in% colnames(summ$coefficients)))

  pred <- predict(
    fit,
    newdata = seq(0, 700, length.out = 6),
    interval = "confidence"
  )
  expect_true(all(c("x", "fit", "se.fit", "working_variance", "marginal_variance", "lwr", "upr") %in% names(pred)))

  fit_auto <- cluster_mm(
    data = sim_dat,
    s = "substrate_conc",
    v = "activity",
    cluster = "cluster",
    variance = "auto",
    power_selection = "quasi_aic"
  )
  expect_s3_class(fit_auto, "cluster_mm")
  expect_identical(fit_auto$variance_selection_mode, "auto")
  expect_true(is.data.frame(fit_auto$variance_candidate_table))

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- report_mm(fit, interval_type = "confidence")
  expect_s3_class(out, "cluster_mm")

  txt_boot <- capture.output({
    set.seed(1)
    report_mm(
      fit,
      method = "bootstrap",
      B = 25,
      bootstrap_type = "pairs",
      bootstrap_ci = "basic",
      interval_type = "none"
    )
  })
  expect_true(any(grepl("CI method: bootstrap \\(cluster pairs", txt_boot)))
})

test_that("plot.cluster_mm expands limits to include interval bands", {
  fit <- cluster_mm(
    data = subset(alves_demo, enzyme == "BG"),
    s = "substrate_conc",
    v = "activity",
    cluster = "core",
    variance = "sqrt"
  )
  pred <- predict(
    fit,
    newdata = seq(min(fit$x), max(fit$x), length.out = 200),
    interval = "prediction"
  )

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(fit, interval_type = "prediction")
  usr <- graphics::par("usr")

  expect_lte(usr[3], min(pred$lwr) + 1e-8)
  expect_gte(usr[4], max(pred$upr) - 1e-8)
})

test_that("cluster_mm applies the sparse-cluster interval policy", {
  fit <- cluster_mm(
    data = subset(alves_demo, enzyme == "BG"),
    s = "substrate_conc",
    v = "activity",
    cluster = "core",
    variance = "sqrt"
  )

  summ <- summary(fit)
  expect_identical(summ$ci_method, "disabled")
  expect_identical(colnames(summ$coefficients), "Estimate")
  expect_true(isTRUE(summ$small_cluster_policy$disable))

  expect_warning(
    ci_small <- confint(fit),
    "Small-cluster warning"
  )
  expect_equal(dim(ci_small), c(2, 2))

  expect_warning(
    {
      set.seed(1)
      ci_boot_small <- confint(
        fit,
        method = "bootstrap",
        B = 25,
        bootstrap_type = "pairs",
        bootstrap_ci = "basic"
      )
    },
    "Experimental clustered bootstrap warning"
  )
  expect_equal(dim(ci_boot_small), c(2, 2))

  txt_small <- capture.output(report_mm(fit, interval_type = "none"))
  expect_true(any(grepl("disabled by default for sparse clustered fits", txt_small, fixed = TRUE)))
})
