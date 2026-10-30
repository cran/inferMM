profile_test_data <- function() {
  x <- seq(1, 100, length.out = 60)
  data.frame(x = x, y = 100 * x / (20 + x) + sin(seq_along(x)))
}

test_that("profile fits verify the score and retain closed-form updates", {
  dat <- profile_test_data()
  h <- function(S) 1 + S / 20
  fit <- fit_mm(dat$x, dat$y, variance = h)
  expect_true(fit$convergence$converged)
  expect_lte(fit$convergence$score_norm, fit$convergence$tolerance)
  expect_identical(fit$convergence$method, "bracketed-root")
  k <- coef(fit)["Km"]
  v <- sum(dat$x * dat$y / ((k + dat$x) * h(dat$x))) /
    sum(dat$x^2 / ((k + dat$x)^2 * h(dat$x)))
  expect_equal(unname(coef(fit)["Vmax"]), unname(v), tolerance = 1e-12)
  expect_equal(fit$scale, mean(fit$residuals^2 / h(dat$x)), tolerance = 1e-12)
  grad <- inferMM:::mm_grad(dat$x, coef(fit)["Vmax"], k)
  expected <- fit$scale * solve(crossprod(grad, grad / h(dat$x)))
  expect_equal(unname(vcov(fit)), unname(expected), tolerance = 1e-12)
  expect_true(is.finite(fit$convergence$information_rcond))
})

test_that("all detected positive roots are compared by weighted residual sum of squares", {
  x <- rep(c(1, 2, 5, 10, 20, 40, 80, 160), 3)
  y <- rep(c(67.019789166268, 32.2317409079629, -16.0786088184372,
             -12.2800820929509, 26.5262068298849, 100.242770989792,
             73.1416945754845, 13.5988482438843), 3)
  fit <- fit_mm(x, y, variance = "constant")
  candidates <- fit$convergence$candidates
  accepted <- candidates[candidates$accepted, ]
  expect_gte(nrow(accepted), 2)
  expect_equal(unname(coef(fit)["Km"]), 14.563336, tolerance = 1e-6)
  expect_gt(unname(coef(fit)["Km"]), min(accepted$Km))
  selected <- candidates[candidates$selected, ]
  expect_equal(selected$scaled_weighted_rss, min(accepted$scaled_weighted_rss))
  expect_true(all(is.finite(confint(fit))))
})

test_that("profile scores and unregularized covariance are invariant to measurement units", {
  dat <- profile_test_data()
  ref <- fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20,
                km_bounds = c(1e-8, 1001))
  for (units in list(c(1e-8, 1e8), c(1e8, 1e-8))) {
    x_unit <- units[1]
    y_unit <- units[2]
    h <- function(S) 1 + S / (20 * x_unit)
    fit <- fit_mm(dat$x * x_unit, dat$y * y_unit, variance = h,
                  km_bounds = c(1e-8, 1001) * x_unit)
    parameter_units <- c(y_unit, x_unit)
    expect_equal(unname(coef(fit) / parameter_units), unname(coef(ref)), tolerance = 1e-9)
    expect_equal(fit$scale / y_unit^2, ref$scale, tolerance = 1e-9)
    expect_equal(unname(vcov(fit) / outer(parameter_units, parameter_units)),
                 unname(vcov(ref)), tolerance = 1e-9)
    expect_equal(fit$convergence$information_rcond,
                 ref$convergence$information_rcond, tolerance = 1e-9)
    expect_lte(fit$convergence$score_norm, fit$convergence$tolerance)
  }
})

test_that("minimizing a score without a verified root returns no fit", {
  x <- seq(1, 100, length.out = 60)
  failure <- tryCatch(fit_mm(x, 2 * x, variance = "constant"),
                      inferMM_profile_error = function(e) e)
  expect_s3_class(failure, "inferMM_profile_error")
  expect_s3_class(failure, "inferMM_convergence_error")
  expect_identical(failure$reason, "no_verified_root")
  expect_gt(failure$score_norm, failure$tolerance)
  expect_match(conditionMessage(failure), "no acceptable positive interior Km root could be verified")
  expect_match(conditionMessage(failure), "No fit or confidence intervals are returned")
  expect_false(inherits(failure, "fit_mm"))
})

test_that("score minimization is accepted only after independent score validation", {
  x <- c(1, 2, 5, 10, 20, 40, 80, 160)
  y <- c(67.019789166268, 32.2317409079629, -16.0786088184372,
         -12.2800820929509, 26.5262068298849, 100.242770989792,
         73.1416945754845, 13.5988482438843)
  profile <- inferMM:::find_km_root(x, y, rep(1, 8), c(1, 1000), grid_n = 3)
  expect_identical(profile$convergence$method, "verified-score-minimum")
  expect_lte(profile$convergence$score_norm, profile$convergence$tolerance)
  expect_equal(profile$root, 14.563336, tolerance = 1e-6)
  residuals <- y - profile$vmax * x / (profile$root + x)
  expect_lt(abs(sum(x * residuals / (profile$root + x)^2)), 1e-6)

  lower <- uniroot(function(k) inferMM:::weighted_score_function(k, x, y, rep(1, 8)),
                    c(1, 5), tol = 1e-12)$root
  interior <- inferMM:::find_km_root(x, y, rep(1, 8), c(lower, 1000), grid_n = 3)
  expect_equal(interior$root, 14.563336, tolerance = 1e-6)
  expect_true(any(!interior$convergence$candidates$accepted))
})

test_that("boundary, non-positive and degenerate fits stop before returning intervals", {
  dat <- profile_test_data()
  fit <- fit_mm(dat$x, dat$y, variance = "constant")
  boundary <- tryCatch(fit_mm(dat$x, dat$y, variance = "constant",
                              km_bounds = c(coef(fit)["Km"], 1001)),
                       inferMM_profile_error = function(e) e)
  expect_s3_class(boundary, "inferMM_profile_error")
  expect_identical(boundary$reason, "search_boundary")
  negative <- tryCatch(fit_mm(dat$x, -dat$y, variance = "constant"),
                       inferMM_profile_error = function(e) e)
  expect_s3_class(negative, "inferMM_profile_error")
  expect_identical(negative$reason, "non_positive_parameters")
  expect_error(fit_mm(dat$x, rep(0, nrow(dat)), variance = "constant"),
               class = "inferMM_profile_error")
  expect_error(fit_mm(dat$x, 100 * dat$x / (20 + dat$x), variance = "constant"),
               class = "inferMM_profile_error")
})

test_that("singular scaled information errors instead of adding ridge", {
  convergence <- list(score_norm = 0, tolerance = 1e-7, candidate_count = 1L)
  grad <- cbind(Vmax = 1:20, Km = -2 * (1:20))
  failure <- tryCatch(inferMM:::profile_covariance(grad, rep(1, 20), 1,
                                                  c(1e-8, 1001), convergence),
                      inferMM_profile_error = function(e) e)
  expect_s3_class(failure, "inferMM_profile_error")
  expect_identical(failure$stage, "covariance calculation")
  expect_identical(failure$reason, "singular_information")
  expect_match(conditionMessage(failure), "reciprocal condition number")

  poorly_conditioned <- grad
  poorly_conditioned[, 2] <- poorly_conditioned[, 2] + 1e-5 * sin(1:20)
  failure <- tryCatch(inferMM:::profile_covariance(poorly_conditioned, rep(1, 20), 1,
                                                  c(1e-8, 1001), convergence),
                      inferMM_profile_error = function(e) e)
  expect_identical(failure$reason, "singular_information")
  expect_gt(failure$information_rcond, 0)
  expect_lt(failure$information_rcond, 1e-12)
})

test_that("screening and grouping retain failures without fabricated estimates", {
  x <- seq(1, 100, length.out = 60)
  screen <- screen_mm(x, 2 * x, variance_models = c("constant", "sqrt"), quiet = TRUE)
  expect_length(screen$fits, 0)
  expect_true(all(screen$table$success == 0))
  expect_true(all(is.na(screen$table$Vmax_hat)))
  expect_true(all(grepl("No fit or confidence intervals", screen$table$error_message)))
  dat <- profile_test_data()
  grouped_data <- rbind(data.frame(x = dat$x, y = dat$y, group = "good"),
                        data.frame(x = x, y = 2 * x, group = "no_root"))
  grouped <- group_mm(grouped_data, s = "x", v = "y", groups = "group",
                      variance_models = "constant", quiet = TRUE)
  failed <- grouped$estimates[grouped$estimates$success == 0, ]
  expect_equal(nrow(failed), 1)
  expect_true(all(is.na(failed$Vmax_hat)))
  expect_true(all(is.na(failed$lower_Km)))
  expect_true(all(grepl("No fit or confidence intervals", failed$error_message)))
})

test_that("automatic variance selection records failed candidates and stops if all fail", {
  x <- seq(1, 100, length.out = 60)
  expect_error(fit_mm(x, 2 * x, variance = "auto", power_grid = c(0.2, 0.5)),
               class = "inferMM_profile_error")
  dat <- profile_test_data()
  fit <- fit_mm(dat$x, dat$y, variance = "auto", power_grid = c(0.2, 0.5))
  expect_true(all(fit$variance_candidate_table$success == 1))
  expect_true(all(is.na(fit$variance_candidate_table$error_message)))

  original <- inferMM:::fit_mm_from_spec
  mixed <- with_mocked_bindings(
    fit_mm(dat$x, dat$y, variance = "auto", power_grid = c(0.2, 0.5)),
    fit_mm_from_spec = function(...) {
      args <- list(...)
      if (identical(args$variance, "power") && identical(args$power, 0.5)) {
        stop("Synthetic candidate failure for selection testing.")
      }
      do.call(original, args, quote = TRUE)
    },
    .package = "inferMM"
  )
  failed <- mixed$variance_candidate_table[mixed$variance_candidate_table$success == 0, ]
  expect_equal(nrow(failed), 1)
  expect_match(failed$error_message, "Synthetic candidate failure")
  expect_true(is.na(failed$quasi_aic))
  expect_true(all(is.finite(confint(mixed))))
  expect_false(identical(mixed$variance_power, 0.5))
})
