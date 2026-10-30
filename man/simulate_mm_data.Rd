\name{simulate_mm_data}
\alias{simulate_mm_data}
\title{Simulate Michaelis-Menten Data}
\usage{
simulate_mm_data(n = 100, x = seq(1, 100, length.out = n),
  vmax = 100, km = 20, variance_shape = c("mm", "exp", "hill"),
  error = c("normal", "skewed", "student_t"), s0 = 1, t2 = 9, c = 20,
  a = 0.05, alpha = 2, skew_shape = 4, seed = NULL, df = 10)
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
  \item{seed}{Deprecated and ignored. Call \code{set.seed()} before
  \code{simulate_mm_data()} if reproducibility is needed.}
  \item{df}{Degrees of freedom greater than four for Student-t errors,
  ensuring a finite fourth moment. Ignored for other error choices.}
}
\value{
A data frame containing the simulated observations.
}
\description{
Generate a synthetic Michaelis-Menten dataset with bounded heteroscedastic
variance and Gaussian, centered gamma, or standardized Student-t errors.
}
\details{
The standardized error has mean zero and variance one. For \code{"skewed"},
it is \eqn{(G-a)/\sqrt{a}} with \eqn{G\sim\mathrm{Gamma}(a,1)}
(shape \eqn{a}=\code{skew_shape}, scale one). For \code{"student_t"},
it is \eqn{T\sqrt{(\nu-2)/\nu}} with \eqn{T\sim t_\nu}.
These are generating distributions, not options that select the fitting
algorithm. See \code{\link{fit_mm_moments}} for higher-moment estimation.
}
\examples{
set.seed(1)
sim_dat <- simulate_mm_data(variance_shape = "hill", error = "skewed")
head(sim_dat)
heavy_tail_dat <- simulate_mm_data(error = "student_t", df = 10)
}
