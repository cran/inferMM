variance_shape_mm <- function(x, s0 = 1, t2 = 9, c = 20) {
  s0 + t2 * x / (c + x)
}

variance_shape_exp <- function(x, s0 = 1, t2 = 9, a = 0.05) {
  s0 + t2 * (1 - exp(-a * x))
}

variance_shape_hill <- function(x, s0 = 1, t2 = 9, c = 20, alpha = 2) {
  s0 + t2 * x^alpha / (c^alpha + x^alpha)
}

draw_standardized_error <- function(n, error = c("normal", "skewed"), skew_shape = 4) {
  error <- match.arg(error)

  if (identical(error, "normal")) {
    return(stats::rnorm(n))
  }

  z <- stats::rgamma(n, shape = skew_shape, scale = 1)
  (z - skew_shape) / sqrt(skew_shape)
}

#' Simulate Michaelis-Menten Data
#'
#' Generate a synthetic Michaelis-Menten dataset with bounded heteroscedastic
#' variance and either Gaussian or moderately skewed errors.
#'
#' @param n Number of observations when `x` is not supplied.
#' @param x Optional substrate concentration grid.
#' @param vmax True `V_max` value.
#' @param km True `K_m` value.
#' @param variance_shape Character string naming the true variance pattern.
#' @param error Character string naming the error distribution.
#' @param s0,t2,c,a,alpha Variance-shape parameters.
#' @param skew_shape Shape parameter for the skewed gamma-based error.
#' @param seed Deprecated and ignored. Call [stats::set.seed()] before
#'   `simulate_mm_data()` if reproducibility is needed.
#'
#' @return A data frame containing the simulated observations.
#' @export
simulate_mm_data <- function(n = 100,
                             x = seq(1, 100, length.out = n),
                             vmax = 100,
                             km = 20,
                             variance_shape = c("mm", "exp", "hill"),
                             error = c("normal", "skewed"),
                             s0 = 1,
                             t2 = 9,
                             c = 20,
                             a = 0.05,
                             alpha = 2,
                             skew_shape = 4,
                             seed = NULL) {
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
  z <- draw_standardized_error(length(x), error = error, skew_shape = skew_shape)
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
