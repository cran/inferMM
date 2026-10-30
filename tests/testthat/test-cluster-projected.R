projected_fixture <- function(m = 50L, distribution = "gamma", unbalanced = FALSE) {
  b <- switch(distribution,
    gamma = 12 * (rgamma(m, shape = 2) - 2) / sqrt(2),
    normal = rnorm(m, sd = 12),
    mixture = ifelse(runif(m) < 0.8, -6, 24) + rnorm(m, sd = 3)
  )
  x <- c(2, 5, 10, 20, 50, 100, 200)
  dat <- do.call(rbind, lapply(seq_len(m), function(i) {
    s <- if (unbalanced && i %% 2L == 0L) x[-2L] else x
    data.frame(id = i, S = s, Y = (100 + b[i]) * s / (20 + s) +
                 rnorm(length(s), sd = sqrt(2 * sqrt(s))))
  }))
  row.names(dat) <- NULL
  dat
}

test_that("conditional projected equations are centered without a latent density", {
  x <- c(2, 5, 10, 20, 50, 100, 200)
  k <- 20
  h <- sqrt(x)
  gamma <- 2
  z <- x / (k + x)
  dz <- -x / (k + x)^2
  c_i <- sum(z^2 / h)
  # Symmetric cubature reproduces E(e) = 0 and Cov(e) = gamma H exactly.
  errors <- rbind(diag(sqrt(length(x) * gamma * h)),
                  -diag(sqrt(length(x) * gamma * h)))
  for (amplitude in c(-10, 80, 100, 160)) {
    y <- sweep(errors, 2L, z * amplitude, "+")
    s <- inferMM:::cluster_projected_state(k,
      rep(list(x), nrow(y)), lapply(seq_len(nrow(y)), function(i) y[i, ]),
      rep(list(h), nrow(y)))$clusters
    expect_equal(mean(s[, "W"]), amplitude, tolerance = 1e-12)
    expect_equal(mean(s[, "W"] * s[, "q"]), 0, tolerance = 1e-11)
    expect_equal(mean(s[, "rss"]), gamma * (length(x) - 1), tolerance = 1e-11)
    expect_equal(mean((s[, "W"] - amplitude)^2), gamma / c_i, tolerance = 1e-11)
  }
  y <- 110 * z + seq(-1, 1, length.out = length(x))
  s <- inferMM:::cluster_projected_state(k, list(x), list(y), list(h))$clusters
  u <- y - z * as.numeric(s[, "W"])
  expect_equal(as.numeric(s[, "q"]), sum(dz * u / h), tolerance = 1e-12)
  step <- 1e-4
  loss <- function(km) inferMM:::cluster_projected_state(km, list(x), list(y), list(h))$objective
  expect_equal((loss(k + step) - loss(k - step)) / (2 * step),
               -2 * as.numeric(s[, "W"] * s[, "q"]), tolerance = 1e-7)
})

test_that("projected fits handle normal, skewed and mixture random effects", {
  set.seed(351)
  for (distribution in c("normal", "gamma", "mixture")) {
    dat <- projected_fixture(distribution = distribution, unbalanced = TRUE)
    fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = "projected")
    expect_s3_class(fit, "cluster_mm")
    expect_identical(fit$estimation_method, "projected")
    expect_identical(fit$covariance_method, "cluster-sandwich")
    expect_lt(abs(coef(fit)["Km"] - 20), 3)
    expect_lt(abs(coef(fit)["Vmax"] - 100), 10)
    expect_lt(fit$convergence$score_norm, 1e-7)
    expect_lt(max(abs(colMeans(fit$score_contributions))), 1e-7)
    expect_true(all(eigen(vcov(fit), symmetric = TRUE)$values > 0))
    expect_equal(fit$call[[1]], quote(cluster_mm))
    expect_equal(summary(fit)$covariance_method, "cluster-sandwich")
    ci <- confint(fit)
    se <- sqrt(diag(vcov(fit)))
    expect_equal(unname(ci[, 1]), unname(coef(fit) - qnorm(0.975) * se))
    expect_equal(unname(ci[, 2]), unname(coef(fit) + qnorm(0.975) * se))
    pred <- predict(fit, c(10, 20, 80), interval = "confidence")
    expect_true(all(is.finite(pred$lwr)))
    expect_error(predict(fit, c(10, 20), interval = "prediction"),
                 "specified random-effect distribution")
    inv <- solve(fit$bread)
    cov_eta <- inv %*% fit$meat %*% t(inv) / fit$n_clusters
    transform <- diag(coef(fit))
    expected <- transform %*% cov_eta[1:2, 1:2] %*% transform
    expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-10)
  }
})

test_that("projected root search verifies scores and compares roots", {
  state <- function(k) list(score = (log(k) - log(2)) * (log(k) - log(10)),
                            objective = (log(k) - log(10))^2)
  root <- inferMM:::cluster_projected_root(state, c(0.1, 100))
  expect_equal(root$km, 10, tolerance = 1e-8)
  bad <- function(k) list(score = 1, objective = (log(k) - 1)^2)
  expect_error(inferMM:::cluster_projected_root(bad, c(0.1, 100)),
               class = "inferMM_cluster_projected_error")
  boundary <- function(k) list(score = log(k), objective = log(k)^2)
  expect_error(inferMM:::cluster_projected_root(boundary, c(1, 100)),
               "no positive interior root")
  maximum <- function(k) list(score = -log(k), objective = -log(k)^2)
  expect_error(inferMM:::cluster_projected_root(maximum, c(0.1, 100)),
               "search boundary has a lower")
})

test_that("projected fits preserve units, refit method and bootstrap routing", {
  set.seed(352)
  dat <- projected_fixture(m = 30L)
  fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = "projected",
                    km_bounds = c(0.01, 1000))
  scaled <- transform(dat, S = S * 1e5, Y = Y * 1e-3)
  scaled_fit <- cluster_mm(scaled, "S", "Y", "id", variance = "sqrt", method = "projected",
                           km_bounds = c(0.01, 1000) * 1e5)
  conversion <- diag(c(1e-3, 1e5))
  expect_equal(unname(coef(scaled_fit)), unname(coef(fit) * c(1e-3, 1e5)), tolerance = 1e-7)
  expect_equal(unname(vcov(scaled_fit)),
               unname(conversion %*% vcov(fit) %*% conversion), tolerance = 1e-6)
  refit <- inferMM:::refit_cluster_mm(fit)
  expect_equal(coef(refit), coef(fit))
  expect_equal(vcov(refit), vcov(fit))
  expect_identical(refit$km_bounds, fit$km_bounds)
  resampled <- inferMM:::resample_cluster_pairs_data(fit$cleaned_data)
  boot_fit <- inferMM:::refit_cluster_mm(fit, resampled)
  expect_identical(boot_fit$estimation_method, "projected")
  expect_identical(boot_fit$covariance_method, "cluster-sandwich")
  expect_warning(boot_ci <- confint(fit, method = "bootstrap", B = 20,
                                    bootstrap_ci = "basic"),
                  "Experimental clustered bootstrap")
  expect_true(all(is.finite(boot_ci)))
  auto <- cluster_mm(dat, "S", "Y", "id", variance = "auto", power_grid = c(0.4, 0.5),
                     method = "projected")
  expect_identical(auto$estimation_method, "projected")
  expect_identical(auto$variance_selection_mode, "auto")
  expect_equal(nrow(auto$variance_candidate_table), 3L)
  expect_identical(inferMM:::refit_cluster_mm(auto)$estimation_method, "projected")
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  tryCatch(expect_silent(plot(fit)), finally = grDevices::dev.off())
})

test_that("projected fitting fails explicitly for invalid or unidentified data", {
  set.seed(353)
  dat <- projected_fixture(m = 10L)
  run <- function(d, ...) cluster_mm(d, "S", "Y", "id", variance = "sqrt",
                                     method = "projected", ...)
  expect_error(run(subset(dat, id <= 4)), "at least five independent clusters")
  expect_error(run(dat, km_bounds = c(1, 0)), "positive increasing")
  expect_error(run(dat, interval_level = 1), "strictly between")
  expect_error(run(dat, km_bounds = c(0.01, 0.1)), "no positive interior root")
  zeros <- dat
  zeros$S[zeros$id == 1] <- 0
  expect_error(run(zeros), "no positive substrate loading")
  flat <- dat
  flat$S <- rep(seq_len(10), each = 7)
  expect_error(run(flat), "no cluster spans two distinct positive")
  exact <- transform(dat, Y = 100 * S / (20 + S))
  expect_error(run(exact), class = "inferMM_cluster_projected_error")
  expect_error(cluster_mm(dat, "S", "Y", "id", km_bounds = c(1, 100)),
               "only with method = 'projected'")
})

test_that("legacy clustered default remains the profile estimator", {
  set.seed(354)
  dat <- projected_fixture(m = 12L, distribution = "normal")
  default <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt")
  explicit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = "profile")
  expect_identical(default$estimation_method, "profile")
  expect_equal(coef(default), coef(explicit))
  expect_equal(vcov(default), vcov(explicit))
  expect_identical(inferMM:::refit_cluster_mm(default)$estimation_method, "profile")
})

test_that("projected failures carry actionable numerical diagnostics", {
  state <- function(k) list(score = 1, objective = (log(k) - 1)^2)
  e <- tryCatch(inferMM:::cluster_projected_root(state, c(.1, 100)), error = identity)
  expect_s3_class(e, "inferMM_convergence_error")
  expect_identical(e$reason, "no_acceptable_root")
  expect_equal(e$diagnostics$finite_grid_points, 400)
  expect_equal(e$diagnostics$accepted_roots, 0)
  expect_equal(e$diagnostics$min_abs_grid_score, 1)
  expect_match(conditionMessage(e), "No fit or confidence intervals")
  expect_match(conditionMessage(e), "No automatic fallback")
  expect_match(conditionMessage(e), "finite_grid_points=400")
  expect_error(inferMM:::cluster_projected_root(state, c(20, 20 * (1 + 1e-10))),
    "too narrow")
  expect_equal(inferMM:::cluster_projected_rms(c(-1, 1) * 1e200), 1e200)
})

test_that("invalid tail evaluations do not discard valid projected roots", {
  state <- function(k) {
    if (k > 100) stop("non-finite tail projection")
    list(score = log(k / 10), objective = log(k / 10)^2)
  }
  root <- inferMM:::cluster_projected_root(state, c(.1, 1e200))
  expect_equal(root$km, 10, tolerance = 1e-8)
  expect_lt(root$diagnostics$finite_grid_points, 400)
  expect_gt(root$diagnostics$accepted_roots, 0)
  set.seed(355)
  dat <- projected_fixture(m = 30L)
  fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = "projected")
  expect_warning(wide <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt",
    method = "projected", km_bounds = c(.01, 1e300)), class = "inferMM_cluster_projected_warning")
  expect_equal(coef(wide), coef(fit), tolerance = 1e-7)
  expect_lt(wide$convergence$joint_score_norm, 1e-7)
  expect_match(wide$convergence$search_warning, "grid scores were finite")
  large <- transform(dat, S = S * 1e305)
  e <- tryCatch(cluster_mm(large, "S", "Y", "id", variance = "sqrt",
    method = "projected", km_bounds = c(1e-100, 1e308)), error = identity)
  expect_identical(e$reason, "invalid_scaled_bounds")
  expect_match(conditionMessage(e), "compatible units")
  overflow <- tryCatch(cluster_mm(dat, "S", "Y", "id", variance = "power",
    power = 1000, method = "projected"), error = identity)
  expect_identical(overflow$reason, "nonfinite_variance_shape")
  expect_match(conditionMessage(overflow), "working variance shape is non-finite")
  expect_gt(overflow$diagnostics$nonfinite_variance_values, 0)
})

test_that("a variance-component boundary is a warning, not failed mean convergence", {
  x <- c(2, 5, 10, 20, 50, 100, 200)
  h <- sqrt(x); z <- x / (20 + x)
  errors <- rbind(diag(sqrt(14 * h)), -diag(sqrt(14 * h)))
  amplitude <- as.numeric(errors %*% (z / h) / sum(z^2 / h))
  errors <- errors - .9 * outer(amplitude, z)
  dat <- data.frame(id = rep(seq_len(14), each = 7), S = rep(x, 14),
    Y = as.vector(t(sweep(errors, 2, 100 * z, "+"))))
  expect_warning(fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt",
    method = "projected"), "tau2_unconstrained=")
  expect_true(fit$converged)
  expect_equal(fit$tau2, 0)
  expect_lt(fit$tau2_unconstrained, 0)
  expect_true(all(is.finite(confint(fit))))
})

test_that("legacy nonconvergence warnings do not endorse unverified intervals", {
  set.seed(356)
  dat <- projected_fixture(m = 12L)
  cleaned <- inferMM:::prepare_cluster_mm_data(dat, "S", "Y", "id")
  expect_warning(fit <- inferMM:::fit_cluster_mm_weighted(cleaned,
    inferMM:::resolve_variance_spec("sqrt"), max_iter = 1L, tol = 1e-30),
    "do not treat its estimates or confidence intervals as converged")
  expect_false(fit$converged)
  expect_match(fit$convergence_warning, "Last absolute parameter change")
})
