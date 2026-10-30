variance_shape_mm <- function(x, s0 = 1, t2 = 9, c = 20) {
  s0 + t2 * x / (c + x)
}

variance_shape_exp <- function(x, s0 = 1, t2 = 9, a = 0.05) {
  s0 + t2 * (1 - exp(-a * x))
}

variance_shape_hill <- function(x, s0 = 1, t2 = 9, c = 20, alpha = 2) {
  s0 + t2 * x^alpha / (c^alpha + x^alpha)
}

draw_standardized_error <- function(n, error = c("normal", "skewed", "student_t"),
                                     skew_shape = 4, df = 10) {
  error <- match.arg(error)

  if (identical(error, "normal")) {
    return(stats::rnorm(n))
  }

  if (identical(error, "student_t")) {
    df <- moment_scalar(df, "df", lower = 4)
    return(stats::rt(n, df = df) * sqrt((df - 2) / df))
  }

  skew_shape <- moment_scalar(skew_shape, "skew_shape")
  z <- stats::rgamma(n, shape = skew_shape, scale = 1)
  (z - skew_shape) / sqrt(skew_shape)
}

#' Simulate Michaelis-Menten Data
#'
#' Generate a synthetic Michaelis-Menten dataset with bounded heteroscedastic
#' variance and Gaussian, centered gamma, or standardized Student-t errors.
#'
#' @param n Number of observations when `x` is not supplied.
#' @param x Optional substrate concentration grid.
#' @param vmax True `V_max` value.
#' @param km True `K_m` value.
#' @param variance_shape Character string naming the true variance pattern.
#' @param error Character string naming the error distribution.
#' @param s0,t2,c,a,alpha Variance-shape parameters.
#' @param skew_shape Shape parameter for the skewed gamma-based error.
#' @param seed Deprecated and ignored. Call [set.seed()] before
#'   `simulate_mm_data()` if reproducibility is needed.
#' @param df Degrees of freedom greater than four for Student-t errors.
#'
#' @return A data frame containing the simulated observations.
#' @export
simulate_mm_data <- function(n = 100,
                             x = seq(1, 100, length.out = n),
                             vmax = 100,
                             km = 20,
                             variance_shape = c("mm", "exp", "hill"),
                             error = c("normal", "skewed", "student_t"),
                             s0 = 1,
                             t2 = 9,
                             c = 20,
                             a = 0.05,
                             alpha = 2,
                             skew_shape = 4,
                             seed = NULL,
                             df = 10) {
  variance_shape <- match.arg(variance_shape)
  error <- match.arg(error)

  ignore_seed_argument(seed)

  x <- as.numeric(x)
  variance_fun <- switch(
    variance_shape,
    mm = function(z) variance_shape_mm(z, s0 = s0, t2 = t2, c = c),
    exp = function(z) variance_shape_exp(z, s0 = s0, t2 = t2, a = a),
    hill = function(z) variance_shape_hill(z, s0 = s0, t2 = t2, c = c, alpha = alpha)
  )

  true_mean <- mm_mean(x, vmax, km)
  true_variance <- variance_fun(x)
  z <- draw_standardized_error(length(x), error = error, skew_shape = skew_shape, df = df)
  y <- true_mean + sqrt(true_variance) * z

  data.frame(
    x = x,
    y = y,
    true_mean = true_mean,
    true_variance = true_variance,
    error = z,
    stringsAsFactors = FALSE
  )
}
