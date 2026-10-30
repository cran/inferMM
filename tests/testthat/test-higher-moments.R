moment_test_data <- function(n = 200L) {
  x <- seq(1, 100, length.out = n)
  z <- (stats::qgamma((seq_len(n) - 0.5) / n, shape = 4) - 4) / 2
  z <- z[order((seq_len(n) * 71) %% (n + 1))]
  data.frame(x = x, y = 100 * x / (20 + x) + sqrt(1 + x / 20) * z)
}

test_that("the implemented score equals the conditional-moment matrix formula", {
  x <- seq(1, 50, length.out = 30)
  y <- 100 * x / (20 + x) + sin(x)
  h <- 1 + x / 20
  s <- 0.5 + x / 100
  k <- 3 + 1.5 * s^2
  reference <- c(100, 20, 2)
  eta <- c(0.02, -0.01, 0.03)
  value <- inferMM:::moment_score_state(eta, x, y, h, reference, s, k)
  theta <- value$theta
  expected <- t(vapply(seq_along(x), function(i) {
    mu <- theta[1] * x[i] / (theta[2] + x[i])
    v <- theta[3] * h[i]
    g <- c(x[i] / (theta[2] + x[i]), -theta[1] * x[i] / (theta[2] + x[i])^2)
    D <- rbind(c(g, 0), c(0, 0, h[i]))
    Omega <- rbind(c(v, s[i] * v^1.5), c(s[i] * v^1.5, (k[i] - 1) * v^2))
    as.numeric(theta * (t(D) %*% solve(Omega, c(y[i] - mu, (y[i] - mu)^2 - v))))
  }, numeric(3)))
  expect_equal(unname(value$scores), expected, tolerance = 1e-11)

  z <- (y - value$mean) / sqrt(value$variance)
  a <- value$mean / sqrt(value$variance)
  b <- -a * theta[2] / (theta[2] + x)
  det <- k - 1 - s^2
  q1 <- ((k - 1) * z - s * (z^2 - 1)) / det
  c1 <- (k - 1 - 2 * s * z) / det
  c2 <- (2 * z - s) / det
  analytic <- rbind(
    c(mean(a*q1-a^2*c1), mean(b*q1-a*b*c1), mean(-a*(q1+z*c1)/2)),
    c(mean(b*q1-a*b*c1), mean(b*(1-2*theta[2]/(theta[2]+x))*q1-b^2*c1), mean(-b*(q1+z*c1)/2)),
    c(mean(-a*c2), mean(-b*c2), mean(-z*c2/2))
  )
  numerical <- inferMM:::moment_jacobian(function(par) {
    inferMM:::moment_score_state(par, x, y, h, reference, s, k)$equation
  }, eta)
  expect_equal(unname(numerical), analytic, tolerance = 1e-7)
})

test_that("normal and symmetric working moments recover the legacy mean equations", {
  dat <- moment_test_data()
  h <- function(S) 1 + S / 20
  profile <- fit_mm(dat$x, dat$y, variance = h)
  normal <- fit_mm_moments(dat$x, dat$y, variance = h, moments = "normal", covariance = "model")
  symmetric <- fit_mm_moments(dat$x, dat$y, variance = h, moments = "student_t", df = 10, covariance = "model")
  expect_equal(coef(normal), coef(profile), tolerance = 1e-5)
  expect_equal(vcov(normal), vcov(profile), tolerance = 1e-5)
  expect_equal(normal$scale, profile$scale, tolerance = 1e-5)
  expect_equal(coef(symmetric), coef(normal), tolerance = 1e-7)
  expect_equal(vcov(symmetric), vcov(normal), tolerance = 1e-7)
  expect_gt(symmetric$vcov_joint[3, 3], normal$vcov_joint[3, 3])
})

test_that("asymmetric fits solve all three equations and use the full covariance", {
  dat <- moment_test_data()
  fit <- fit_mm_moments(dat$x, dat$y, variance = function(S) 1 + S / 20,
                         moments = "gamma", shape = 4, covariance = "sandwich")
  expect_s3_class(fit, "fit_mm_moments")
  expect_true(fit$convergence$converged)
  expect_lte(fit$convergence$score_norm, fit$control$tol)
  expect_false(isTRUE(all.equal(coef(fit), coef(fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20)))))
  theta <- c(coef(fit), gamma = fit$scale)
  inv <- solve(fit$bread_log)
  expected <- diag(theta) %*% inv %*%
    (crossprod(fit$estimating_scores_log) / length(fit$y)^2) %*% t(inv) %*% diag(theta)
  expect_equal(unname(fit$vcov_joint), unname(expected), tolerance = 1e-9)
  model <- diag(theta) %*% solve(fit$information_log) %*% diag(theta) / length(fit$y)
  expect_equal(unname(fit$vcov_joint_model), unname(model), tolerance = 1e-9)
  expect_equal(vcov(fit), fit$vcov_joint[1:2, 1:2])
  expect_equal(fit$variance_hat, fit$scale * fit$variance_fun(fit$x))
  expect_equal(unname(confint(fit)[, 1]), unname(coef(fit) - qnorm(.975) * sqrt(diag(vcov(fit)))))
  expect_true(is.na(fit$quasi_aic))
  expect_true(is.na(fit$quasi_bic))
  expect_output(print(fit), "higher-moment")
  expect_output(print(summary(fit)), "sandwich")
})

test_that("higher-moment fits default to joint inverse information", {
  dat <- moment_test_data()
  args <- list(x = dat$x, y = dat$y, variance = function(S) 1 + S / 20,
               moments = "gamma", shape = 4)
  fit <- do.call(fit_mm_moments, args)
  model <- do.call(fit_mm_moments, c(args, list(covariance = "model")))
  sandwich <- do.call(fit_mm_moments, c(args, list(covariance = "sandwich")))
  expect_identical(fit$covariance_type, "model")
  expect_identical(fit$vcov_joint, fit$vcov_joint_model)
  expect_equal(vcov(fit), vcov(model))
  expect_equal(confint(fit), confint(model))
  expect_identical(coef(fit), coef(sandwich))
  expect_identical(fit$scale, sandwich$scale)
  expect_false(isTRUE(all.equal(vcov(fit), vcov(sandwich))))
  expect_output(print(summary(fit)), "Covariance: model")
  refit <- inferMM:::refit_fit_mm(fit)
  expect_identical(refit$covariance_type, "model")
  expect_equal(vcov(refit), vcov(fit))
})

test_that("correct skewed moments have the theoretical information advantage in this design", {
  x <- seq(1, 100, length.out = 200)
  h <- 1 + x / 20
  normal <- inferMM:::moment_score_state(rep(0, 3), x, rep(0, 200), h, c(100, 20, 1), 0, 3)
  gamma <- inferMM:::moment_score_state(rep(0, 3), x, rep(0, 200), h, c(100, 20, 1), 1, 4.5)
  gain <- solve(normal$information)[1:2, 1:2] - solve(gamma$information)[1:2, 1:2]
  expect_gte(min(eigen(gain, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
  expect_gt(sum(diag(gain)), 0)
})

test_that("custom moments, shape scaling, and serialization retain their meaning", {
  dat <- moment_test_data()
  h <- function(S) 1 + S / 20
  fit <- fit_mm_moments(dat$x, dat$y, variance = h, moments = "custom",
                         skewness = function(S) rep(1, length(S)), kurtosis = 4.5)
  gamma <- fit_mm_moments(dat$x, dat$y, variance = h, moments = "gamma")
  expect_equal(coef(fit), coef(gamma))
  expect_equal(vcov(fit), vcov(gamma))
  for (multiplier in c(1e-12, 1e12)) {
    shape_fun <- local({ a <- multiplier; function(S) a * (1 + S / 20) })
    scaled <- fit_mm_moments(dat$x, dat$y, variance = shape_fun, moments = "gamma")
    expect_equal(coef(scaled), coef(gamma), tolerance = 1e-7)
    expect_equal(vcov(scaled), vcov(gamma), tolerance = 1e-7)
    expect_equal(scaled$scale * multiplier, gamma$scale, tolerance = 1e-7)
  }
  refit <- inferMM:::refit_fit_mm(fit)
  expect_s3_class(refit, "fit_mm_moments")
  expect_equal(coef(refit), coef(fit))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(predict(restored, c(1, 20, 100), interval = "confidence"),
                 predict(fit, c(1, 20, 100), interval = "confidence"))
})

test_that("estimated moments are re-estimated in pairs refits without internal seeding", {
  dat <- moment_test_data()
  set.seed(981)
  before <- .Random.seed
  fit <- fit_mm_moments(dat$x, dat$y, variance = function(S) 1 + S / 20)
  expect_identical(.Random.seed, before)
  idx <- seq(1, nrow(dat), by = 2)
  refit <- inferMM:::refit_fit_mm(fit, x = dat$x[idx], y = dat$y[idx])
  expect_s3_class(refit, "fit_mm_moments")
  expect_identical(refit$moment_spec$moments, "estimate")
  expect_false(identical(refit$standardized_moments$skewness[1], fit$standardized_moments$skewness[1]))
  ci <- confint(fit, method = "bootstrap", B = 20, bootstrap_ci = "percentile")
  expect_equal(dim(ci), c(2L, 2L))
  expect_true(all(is.finite(ci)))
  expect_equal(attr(attr(ci, "bootstrap_draws"), "bootstrap_type"), "pairs")
  expect_error(confint(fit, method = "bootstrap", bootstrap_type = "wild"), "pairs")
  expect_error(summary(fit, method = "bootstrap", bootstrap_type = "wild"), "pairs")
  expect_error(predict(fit, 20, interval = "prediction"), "predictive quantiles")
  expect_error(plot(fit, interval_type = "prediction"), "confidence bands")
  expect_error(confint(fit, level = 1), "confidence level")
  expect_error(summary(fit, level = -1), "confidence level")
  expect_error(predict(fit, 20, level = NA_real_), "confidence level")
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit({ grDevices::dev.off(); unlink(output) }, add = TRUE)
  expect_output(report_mm(fit), "higher-moment")
  expect_output(report_mm(fit, method = "bootstrap", B = 20,
                           bootstrap_ci = "percentile"), "bootstrap")
})

test_that("changing concentration and response units preserves higher-moment inference", {
  dat <- moment_test_data()
  h <- function(S) 1 + S / 20
  fit <- fit_mm_moments(dat$x, dat$y, variance = h)
  converted <- fit_mm_moments(dat$x * 1e-3, dat$y * 1e5,
                               variance = function(S) h(S * 1e3))
  units <- c(1e5, 1e-3)
  expect_equal(coef(converted) / units, coef(fit), tolerance = 1e-6)
  expect_equal(vcov(converted) / outer(units, units), vcov(fit), tolerance = 1e-6)
  expect_equal(converted$scale / 1e10, fit$scale, tolerance = 1e-6)
})

test_that("invalid moments and unconverged roots cannot produce confidence intervals", {
  dat <- moment_test_data()
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "custom", skewness = 1, kurtosis = 2), "positive definite")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "custom"), "both skewness and kurtosis")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "custom", skewness = function(S) 1, kurtosis = 4.5), "per concentration")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "custom", skewness = 0, kurtosis = Inf), "finite")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "gamma", shape = -1), "shape")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "student_t", df = 4), "df")
  expect_error(fit_mm_moments(dat$x, dat$y, shape = 4), "only used")
  expect_error(fit_mm_moments(dat$x, dat$y, df = 10), "only used")
  expect_error(fit_mm_moments(dat$x, dat$y, skewness = 0), "only used")
  expect_error(fit_mm_moments(dat$x, dat$y, variance = "auto"), "arg")
  expect_error(fit_mm_moments(dat$x, dat$y, interval_level = 1), "less than 1")
  expect_error(fit_mm_moments(dat$x, dat$y, control = list(unknown = 1)), "control")
  expect_error(fit_mm_moments(dat$x, dat$y, control = list(maxit = 1.5)), "integer")
  expect_error(fit_mm_moments(dat$x, dat$y, moments = "gamma", shape = 2,
                              control = list(maxit = 1, min_step = .25, tol = 1e-10)), "converge")
  expect_error(fit_mm_moments(dat$x[1:8], dat$y[1:8]), "at least 10")
  expect_warning(fit_mm_moments(dat$x[1:40], dat$y[1:40], moments = "estimate"), "small samples")
})

test_that("Newton failures distinguish the iteration limit, stagnation, and non-finite evaluations", {
  control <- inferMM:::moment_control(list(maxit = 1, tol = 1e-10))
  lower <- rep(-30, 3)
  upper <- rep(30, 3)
  limit <- inferMM:::moment_newton(function(eta) exp(eta) - 2,
                                  rep(0, 3), lower, upper, control)
  expect_false(limit$ok)
  expect_identical(limit$reason, "iteration_limit")
  expect_gt(limit$score_norm, control$tol)

  stalled <- inferMM:::moment_newton(function(eta) rep(1, 3),
                                    rep(0, 3), lower, upper, control)
  expect_false(stalled$ok)
  expect_identical(stalled$reason, "no_decreasing_step")
  expect_equal(stalled$score_norm, 1)

  invalid <- inferMM:::moment_newton(function(eta) rep(Inf, 3),
                                    rep(0, 3), lower, upper, control)
  expect_identical(invalid$reason, "non_finite_score")
  derivative <- inferMM:::moment_newton(function(eta) {
    if (all(eta == 0)) rep(1, 3) else rep(Inf, 3)
  }, rep(0, 3), lower, upper, control)
  expect_identical(derivative$reason, "non_finite_jacobian")
})

test_that("failed fits provide actionable diagnostics and no fitted object", {
  dat <- moment_test_data()
  failure <- with_mocked_bindings(
    tryCatch(
      fit_mm_moments(dat$x, dat$y, moments = "gamma", shape = 2,
                     control = list(maxit = 1, min_step = .25, tol = 1e-10)),
      inferMM_convergence_error = function(e) e
    ),
    moment_newton = function(...) {
      list(ok = FALSE, iterations = 1L, score_norm = 1e-3, reason = "iteration_limit")
    },
    .package = "inferMM"
  )
  expect_s3_class(failure, "inferMM_convergence_error")
  expect_false(inherits(failure, "fit_mm"))
  expect_identical(failure$moments, "gamma")
  expect_identical(failure$stage, "initial mean/variance fitting")
  expect_gt(failure$score_norm, failure$tolerance)
  expect_gt(failure$iterations, 0)
  expect_match(conditionMessage(failure), "Total iterations:")
  expect_match(conditionMessage(failure), "max absolute scaled score:")
  expect_match(conditionMessage(failure), "No fit or confidence intervals are returned")
  expect_match(conditionMessage(failure), "control\\$maxit")
  expect_false(grepl("Consider retrying with moments='normal'", conditionMessage(failure),
                     fixed = TRUE))

  continuation <- tryCatch(
    fit_mm_moments(dat$x, dat$y, moments = "gamma", shape = 2,
                   control = list(maxit = 1, min_step = .25, tol = 1e-7)),
    inferMM_convergence_error = function(e) e
  )
  expect_s3_class(continuation, "inferMM_convergence_error")
  expect_identical(continuation$stage, "higher-moment continuation")
  expect_identical(continuation$reason, "iteration_limit")
  expect_equal(continuation$progress, 0)
  expect_gt(continuation$score_norm, continuation$tolerance)
  expect_match(conditionMessage(continuation), "Completed continuation fraction:")
  expect_match(conditionMessage(continuation), "Consider retrying with moments='normal'", fixed = TRUE)
  expect_match(conditionMessage(continuation), "keeping the same variance function", fixed = TRUE)
  expect_match(conditionMessage(continuation), "No automatic fallback is performed", fixed = TRUE)
  expect_false(inherits(continuation, "fit_mm"))
})

test_that("normal-moment retry advice is restricted to non-Gaussian continuation", {
  attempt <- list(reason = "no_decreasing_step", iterations = 3L, score_norm = 0.1)
  control <- inferMM:::moment_control(list())
  failure_message <- function(stage, moments) {
    failure <- tryCatch(
      inferMM:::moment_fit_error(stage, attempt, moments, control),
      inferMM_convergence_error = function(e) e
    )
    conditionMessage(failure)
  }
  advice <- "Consider retrying with moments='normal'"
  for (moments in c("estimate", "gamma", "student_t", "custom")) {
    message <- failure_message("higher-moment continuation", moments)
    expect_match(message, advice, fixed = TRUE)
    expect_match(message, "convergence and efficiency under the original non-Gaussian moments are not guaranteed",
                 fixed = TRUE)
  }
  for (stage in c("initial mean/variance fitting", "boundary validation")) {
    expect_false(grepl(advice, failure_message(stage, "gamma"), fixed = TRUE))
  }
  expect_false(grepl(advice, failure_message("higher-moment continuation", "normal"),
                     fixed = TRUE))
})

test_that("Student-t and centered gamma simulations have the advertised standardization", {
  set.seed(731)
  expected <- rt(100, df = 10) * sqrt(8 / 10)
  set.seed(731)
  actual <- simulate_mm_data(error = "student_t", df = 10)
  expect_equal(actual$error, expected)
  expect_equal(actual$y, actual$true_mean + sqrt(actual$true_variance) * expected)
  set.seed(732)
  expected_gamma <- (rgamma(100, shape = 2) - 2) / sqrt(2)
  set.seed(732)
  expect_equal(simulate_mm_data(error = "skewed", skew_shape = 2)$error, expected_gamma)
  expect_error(simulate_mm_data(error = "student_t", df = 4), "df")
})
