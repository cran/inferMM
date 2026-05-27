\name{simulate_mm_data}
\alias{simulate_mm_data}
\title{Simulate Michaelis-Menten Data}
\usage{
simulate_mm_data(n = 100, x = seq(1, 100, length.out = n),
  vmax = 100, km = 20, variance_shape = c("mm", "exp", "hill"),
  error = c("normal", "skewed"), s0 = 1, t2 = 9, c = 20,
  a = 0.05, alpha = 2, skew_shape = 4, seed = NULL)
}
\arguments{
  \item{n}{Number of observations when \code{x} is not supplied.}
  \item{x}{Optional substrate concentration grid.}
  \item{vmax}{True \eqn{V_{\max}} value.}
  \item{km}{True \eqn{K_m} value.}
  \item{variance_shape}{Character string naming the true variance pattern.}
  \item{error}{Character string naming the error distribution.}
  \item{s0,t2,c,a,alpha}{Variance-shape parameters.}
  \item{skew_shape}{Shape parameter for the skewed gamma-based error.}
  \item{seed}{Deprecated and ignored. Call \code{stats::set.seed()} before
  \code{simulate_mm_data()} if reproducibility is needed.}
}
\value{
A data frame containing the simulated observations.
}
\description{
Generate a synthetic Michaelis-Menten dataset with bounded heteroscedastic
variance and either Gaussian or moderately skewed errors.
}
\examples{
set.seed(1)
sim_dat <- simulate_mm_data(variance_shape = "hill", error = "skewed")
head(sim_dat)
}
