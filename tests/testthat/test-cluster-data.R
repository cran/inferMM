cluster_data_fixture <- function() {
  x <- c(2, 5, 10, 20, 50, 100, 200)
  dat <- data.frame(id = rep(seq_len(60L), each = length(x)), S = rep(x, 60L))
  b <- rnorm(60L, sd = 8)
  dat$Y <- (100 + b[dat$id]) * dat$S / (20 + dat$S) +
    rnorm(nrow(dat), sd = sqrt(2 * sqrt(dat$S)))
  dat
}

test_that("working covariance inverses keep singleton matrix dimensions", {
  for (h in c(0.25, 2, 8)) {
    z <- 0.6
    for (tau2 in c(0, 3)) {
      inverse <- inferMM:::working_cov_inv_re(z, h, tau2, gamma = 0.5)
      expect_identical(dim(inverse), c(1L, 1L))
      expect_equal(inverse, matrix(1 / (tau2 * z^2 + 0.5 * h)), tolerance = 1e-12)
    }
  }
  z <- c(0.2, 0.5, 0.8)
  h <- c(1, 2, 4)
  covariance <- 3 * tcrossprod(z) + diag(0.5 * h)
  expect_equal(inferMM:::working_cov_inv_re(z, h, 3, 0.5), solve(covariance),
               tolerance = 1e-12)
})

test_that("profile fitting and refitting support a singleton cluster", {
  set.seed(20261028)
  dat <- cluster_data_fixture()
  dat <- dat[dat$id != 1L | seq_len(nrow(dat)) == 1L, ]
  expect_silent(fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt"))
  expect_true(fit$converged)
  expect_identical(fit$n_clusters, 60L)
  expect_identical(min(fit$cluster_sizes), 1L)
  expect_equal(sum(fit$cluster_sizes), nrow(dat))
  expect_equal(length(fit$y), 414L)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(confint(fit))))
  expect_equal(summary(fit)$n_clusters, 60L)
  expect_silent(refit <- inferMM:::refit_cluster_mm(fit))
  expect_equal(coef(refit), coef(fit))
  expect_equal(vcov(refit), vcov(fit))
  expect_error(cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = "projected"),
               "two readings per cluster")
})

test_that("missing cluster IDs are omitted consistently for both estimators", {
  set.seed(20261028)
  dat <- cluster_data_fixture()
  dat$id[1] <- NA_integer_
  for (method in c("profile", "projected")) {
    expect_warning(fit <- cluster_mm(dat, "S", "Y", "id", variance = "sqrt", method = method),
                   "Removed 1 observation\\(s\\) with missing cluster identifiers")
    reference <- cluster_mm(dat[!is.na(dat$id), ], "S", "Y", "id",
                            variance = "sqrt", method = method)
    expect_true(fit$converged)
    expect_identical(fit$n_clusters, 60L)
    expect_equal(nrow(fit$cluster_reference), 60L)
    expect_equal(sum(fit$cluster_sizes), 419L)
    expect_equal(length(fit$y), 419L)
    expect_equal(length(fitted(fit)), 419L)
    expect_equal(length(residuals(fit)), 419L)
    expect_false(anyNA(fit$cluster))
    expect_false(anyNA(fit$cleaned_data$.cluster_id))
    expect_equal(coef(fit), coef(reference), tolerance = 0)
    expect_equal(vcov(fit), vcov(reference), tolerance = 0)
    expect_equal(confint(fit), confint(reference), tolerance = 0)
    expect_equal(fit$scale, reference$scale, tolerance = 0)
    expect_equal(fit$tau2, reference$tau2, tolerance = 0)
    expect_equal(summary(fit)$n_clusters, 60L)
    expect_output(print(fit), "Clusters: 60")
    expect_silent(refit <- inferMM:::refit_cluster_mm(fit))
    expect_equal(coef(refit), coef(fit))
    expect_equal(refit$n_clusters, 60L)
  }
})

test_that("multi-column IDs drop missing components and unused factor levels", {
  set.seed(20261028)
  dat <- cluster_data_fixture()
  dat$plate <- factor(ifelse(dat$id <= 30, "A", "B"), levels = c("A", "B", "unused"))
  dat$id <- factor(dat$id, levels = seq_len(65L))
  dat$id[1:7] <- NA
  dat$plate[8] <- NA
  for (method in c("profile", "projected")) {
    expect_warning(fit <- cluster_mm(dat, "S", "Y", c("plate", "id"),
                                     variance = "sqrt", method = method, sort_x = FALSE),
                   "Removed 8 observation\\(s\\)")
    expect_identical(fit$n_clusters, 59L)
    expect_equal(sum(fit$cluster_sizes), nrow(dat) - 8L)
    expect_equal(length(fit$cluster_sizes), 59L)
    expect_true(all(fit$cluster_sizes > 0))
    expect_equal(fit$cleaned_data$.row_id, seq.int(9L, nrow(dat)))
    expect_false(anyNA(fit$cleaned_data$.cluster_id))
    expect_false(anyNA(fit$cleaned_data$plate))
    expect_false(anyNA(fit$cleaned_data$id))
  }
})

test_that("missing ID validation uses only retained observations and clusters", {
  set.seed(20261028)
  dat <- cluster_data_fixture()
  nan <- dat
  nan$id[1] <- NaN
  expect_warning(cleaned <- inferMM:::prepare_cluster_mm_data(nan, "S", "Y", "id"),
                 "Removed 1 observation\\(s\\)")
  expect_equal(nrow(cleaned), nrow(dat) - 1L)

  literal <- dat
  literal$id <- as.character(literal$id)
  literal$id[literal$id == "1"] <- "NA"
  expect_silent(cleaned <- inferMM:::prepare_cluster_mm_data(literal, "S", "Y", "id"))
  expect_equal(nrow(cleaned), nrow(dat))
  expect_equal(length(unique(cleaned$.cluster_id)), 60L)

  empty <- transform(dat, id = NA_integer_)
  one <- dat[dat$id <= 2L, ]
  one$id[one$id == 2L] <- NA_integer_
  four <- dat[dat$id <= 5L, ]
  four$id[four$id == 5L] <- NA_integer_
  for (method in c("profile", "projected")) {
    expect_warning(expect_error(cluster_mm(empty, "S", "Y", "id", method = method),
                                "No valid observations remain after cleaning"),
                   "Removed 420 observation\\(s\\)")
    expect_warning(expect_error(cluster_mm(one, "S", "Y", "id", method = method),
                                "At least two clusters are required"),
                   "Removed 7 observation\\(s\\)")
  }
  expect_warning(expect_error(cluster_mm(four, "S", "Y", "id", method = "projected"),
                              "at least five independent clusters"),
                 "Removed 7 observation\\(s\\)")
})
