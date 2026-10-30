\name{inferMM-package}
\alias{inferMM-package}
\alias{inferMM}
\title{inferMM: Variance-Aware Michaelis-Menten Estimation and Inference}
\description{
The \pkg{inferMM} package provides variance-aware tools for Michaelis-Menten
estimation, clustered repeated-measurement analysis, model screening, grouped
curve analysis, simulation, and prediction.
Single-curve fits and grouped screening accept user-defined fixed positive
variance functions in addition to the built-in variance shapes.
An optional higher-moment estimator fits independent curves using working
skewness and kurtosis, with joint sandwich or model-based covariance.
}
\details{
The package focuses on enzyme-kinetic settings where measurement variability may
increase with substrate concentration, making standard homoscedastic nonlinear
least squares inadequate for interval estimation and uncertainty quantification.
}
\seealso{
\code{\link{fit_mm}}, \code{\link{fit_mm_moments}},
\code{\link{cluster_mm}}, \code{\link{report_mm}},
\code{\link{screen_mm}}, \code{\link{group_mm}},
\code{\link{simulate_mm_data}}, \code{\link{variance_function}}
}
\keyword{package}
