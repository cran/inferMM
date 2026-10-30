custom_curve_data <- function() {
  x <- seq(0, 100, length.out = 60)
  data.frame(x = x, y = 100 * x / (20 + x) + sin(seq_along(x)))
}

test_that("custom fixed shapes use the profile score and closed-form updates", {
  dat <- custom_curve_data()
  h <- function(S) 1 + S / 20
  fit <- fit_mm(dat$x, dat$y, variance = h)
  ref <- stats::nls(y ~ Vmax * x / (Km + x), data = dat,
                    weights = 1 / h(dat$x), start = list(Vmax = 100, Km = 20))
  expect_identical(fit$method, "profile-score")
  expect_identical(fit$variance_model, "custom")
  expect_equal(coef(fit), coef(ref), tolerance = 1e-5)
  km <- coef(fit)["Km"]
  vmax <- sum(dat$x * dat$y / ((km + dat$x) * h(dat$x))) /
    sum(dat$x^2 / ((km + dat$x)^2 * h(dat$x)))
  expect_equal(unname(coef(fit)["Vmax"]), unname(vmax))
  expect_equal(fit$scale, mean(fit$residuals^2 / h(dat$x)))
  expect_equal(fit$variance_hat, fit$scale * h(dat$x))
  expect_equal(fit$variance_fun(dat$x), h(dat$x))
  expect_true(all(is.finite(vcov(fit))))
  expect_equal(dim(confint(fit)), c(2L, 2L))
  expect_true(all(is.finite(summary(fit)$coefficients)))
  expect_equal(fit$quasi_aic, -2 * fit$quasi_loglik + 6)
  expect_equal(fit$quasi_bic, -2 * fit$quasi_loglik + 3 * log(nrow(dat)))
})

test_that("custom constant and square-root shapes agree with built-in fits", {
  dat <- custom_curve_data()[-1, ]
  for (model in c("constant", "sqrt")) {
    builtin <- fit_mm(dat$x, dat$y, variance = model)
    custom <- fit_mm(dat$x, dat$y, variance = variance_function(model))
    expect_equal(coef(custom), coef(builtin), tolerance = 1e-6)
    expect_equal(vcov(custom), vcov(builtin), tolerance = 1e-5)
    expect_equal(custom$variance_hat, builtin$variance_hat, tolerance = 1e-6)
    expect_equal(custom$quasi_aic, builtin$quasi_aic, tolerance = 1e-6)
  }
})

test_that("positive shape rescaling preserves fits without clipping", {
  dat <- custom_curve_data()
  ref <- fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20)
  for (multiplier in c(1e-14, 1e14)) {
    h <- local({
      a <- multiplier
      function(S) a * (1 + S / 20)
    })
    fit <- fit_mm(dat$x, dat$y, variance = h)
    expect_equal(coef(fit), coef(ref), tolerance = 1e-7)
    expect_equal(vcov(fit), vcov(ref), tolerance = 1e-7)
    expect_equal(fit$scale * multiplier, ref$scale, tolerance = 1e-7)
    expect_equal(fit$variance_hat, ref$variance_hat, tolerance = 1e-7)
    expect_equal(fit$quasi_aic, ref$quasi_aic, tolerance = 1e-7)
    expect_equal(predict(fit, c(0, 20, 80), interval = "prediction"),
                 predict(ref, c(0, 20, 80), interval = "prediction"),
                 tolerance = 1e-7)
  }
})

test_that("custom functions are validated rather than recycled or clipped", {
  dat <- custom_curve_data()
  invalid <- list(
    function(S) rep(0, length(S)),
    function(S) rep(-1, length(S)),
    function(S) rep(NA_real_, length(S)),
    function(S) rep(NaN, length(S)),
    function(S) rep(Inf, length(S)),
    function(S) rep(TRUE, length(S)),
    function(S) rep(1 + 1i, length(S)),
    function(S) matrix(1, nrow = length(S)),
    function(S) as.list(rep(1, length(S))),
    function(S) "positive",
    function(S) 1,
    function(S) numeric(0),
    function(S) rep(1, length(S) + 1L)
  )
  for (h in invalid) {
    expect_error(fit_mm(dat$x, dat$y, variance = h), "custom variance function")
  }
  expect_error(fit_mm(dat$x, dat$y, variance = function(S) stop("bad shape")),
               "Custom variance function failed: bad shape")
  expect_error(fit_mm(dat$x, dat$y, variance = function(S, a) a + S),
               "Custom variance function failed")
  expect_error(fit_mm(dat$x, dat$y, variance = sqrt), "strictly positive")
  expect_error(fit_mm(dat$x, dat$y, variance = function(S) 1 + S, power = 1),
               "power cannot be supplied")
})

test_that("custom labels and closure contents survive serialization", {
  dat <- custom_curve_data()
  h <- local({
    offset <- 2
    variance_function(function(S) offset + S / 20, label = "Assay variance")
  })
  fit <- fit_mm(dat$x, dat$y, variance = h)
  expect_identical(fit$variance_label, "Assay variance")
  expect_identical(summary(fit)$variance_label, "Assay variance")
  expect_output(print(fit), "Assay variance")
  expect_equal(h(c(0, 20)), c(2, 3))
  expect_equal(variance_function(h, label = "New label")(c(0, 20)), c(2, 3))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(predict(restored, c(0, 10, 200), interval = "prediction"),
               predict(fit, c(0, 10, 200), interval = "prediction"))
  for (label in list("", NA_character_, c("a", "b"), 1)) {
    expect_error(variance_function(function(S) 1 + S, label = label), "label")
  }
  expect_error(variance_function("sqrt", label = "custom"), "only used")
})

test_that("new concentrations use the custom shape and are validated", {
  dat <- custom_curve_data()
  h <- function(S) ifelse(S <= 100, 1 + S / 20, NA_real_)
  fit <- fit_mm(dat$x, dat$y, variance = h)
  pred <- predict(fit, c(0, 10, 40), interval = "prediction")
  expect_equal(pred$working_variance, fit$scale * h(pred$x))
  expect_equal(pred$upr - pred$fit,
               qnorm(0.975) * sqrt(pred$se.fit^2 + pred$working_variance))
  expect_error(predict(fit, c(10, 101), interval = "prediction"), "strictly positive")
  expect_equal(predict(fit, interval = "prediction")$working_variance,
               fit$variance_hat)
})

test_that("custom functions survive wild and pairs bootstrap refits", {
  dat <- custom_curve_data()
  fit <- fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20)
  for (scheme in c("wild", "pairs")) {
    set.seed(44)
    ci <- confint(fit, method = "bootstrap", bootstrap_type = scheme, B = 25)
    expect_true(all(is.finite(ci)))
    draws <- attr(ci, "bootstrap_draws")
    expect_identical(attr(draws, "successful_B"), 25L)
    expect_identical(attr(draws, "bootstrap_type"), scheme)
  }
  refit <- inferMM:::refit_fit_mm(fit, x = rev(dat$x), y = rev(dat$y))
  expect_identical(refit$variance_label, fit$variance_label)
  expect_equal(refit$variance_fun(c(0, 30)), c(1, 2.5))
  expect_equal(coef(refit), coef(fit), tolerance = 1e-7)
})

test_that("custom screening preserves names, labels, and failed candidates", {
  dat <- custom_curve_data()
  out <- screen_mm(dat$x, dat$y, variance_models = "constant",
                   variance_functions = list(
                     linear_floor = function(S) 1 + S / 20,
                     bad = function(S) rep(0, length(S)),
                     saturation = variance_function(function(S) 1 + S / (20 + S),
                                                    label = "Saturating variance")
                   ), quiet = TRUE)
  expect_setequal(out$table$model, c("constant", "linear_floor", "bad", "saturation"))
  expect_setequal(names(out$fits), c("constant", "linear_floor", "saturation"))
  tab <- out$table
  expect_equal(tab$success[tab$model == "bad"], 0)
  expect_match(tab$error_message[tab$model == "bad"], "strictly positive")
  expect_identical(tab$selected_model[tab$model == "linear_floor"], "linear_floor")
  expect_identical(tab$selected_model_label[tab$model == "saturation"], "Saturating variance")
  expect_true(all(diff(tab$quasi_aic[tab$success == 1]) >= 0))

  only <- screen_mm(dat$x, dat$y, variance_models = character(0),
                    variance_functions = list(assay = function(S) 1 + S / 20),
                    quiet = TRUE)
  expect_identical(only$table$model, "assay")
  failed <- screen_mm(dat$x, dat$y, variance_models = character(0),
                      variance_functions = list(bad = function(S) 0), quiet = TRUE)
  expect_length(failed$fits, 0)
  expect_equal(failed$table$success, 0)
  expect_output(print(failed), "bad")
})

test_that("candidate configuration errors are explicit", {
  dat <- custom_curve_data()
  h <- function(S) 1 + S
  bad_lists <- list(list(h), list(a = h, a = h), list(a = 1),
                    stats::setNames(list(h), ""), h, "a", numeric(0))
  for (value in bad_lists) {
    expect_error(screen_mm(dat$x, dat$y, variance_functions = value, quiet = TRUE),
                 "named list")
  }
  expect_error(screen_mm(dat$x, dat$y, variance_functions = list(sqrt = h)),
               "duplicate built-in")
  expect_error(screen_mm(dat$x, dat$y, variance_models = character(0)),
               "at least one")
})

test_that("grouping, reporting, and plotting accept custom fits", {
  dat <- custom_curve_data()
  grouped_data <- rbind(data.frame(dat, group = "A"), data.frame(dat, group = "B"))
  grouped <- group_mm(grouped_data, s = "x", v = "y", groups = "group",
                      variance_models = character(0), quiet = TRUE,
                      variance_functions = list(assay = function(S) 1 + S / 20))
  expect_equal(nrow(grouped$comparison$best_by_group), 2)
  expect_true(all(grouped$estimates$selected_model == "assay"))
  expect_true(all(is.finite(grouped$estimates$se_Vmax)))
  expect_true(all(is.finite(grouped$estimates$lower_Km)))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_s3_class(plot(grouped, model = "assay"), "group_mm")
  fit <- fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20)
  expect_output(report_mm(fit, interval_type = "confidence"), "profile-score")
  expect_error(cluster_mm(grouped_data, "x", "y", "group", variance = function(S) 1 + S),
               "not cluster_mm", fixed = TRUE)
})

test_that("deterministic custom fitting does not change the random seed", {
  dat <- custom_curve_data()
  set.seed(12)
  before <- .Random.seed
  fit_mm(dat$x, dat$y, variance = function(S) 1 + S / 20)
  expect_identical(.Random.seed, before)
})
