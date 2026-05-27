\name{variance_function}
\alias{variance_function}
\title{Create Working Variance Functions}
\usage{
variance_function(model = c("constant", "log", "sqrt", "cuberoot", "power"),
  power = NULL)
}
\arguments{
  \item{model}{Character string naming the working variance specification.}
  \item{power}{Optional exponent used when \code{model = "power"}.}
}
\value{
A function of substrate concentration \code{x}.
}
\description{
Construct a working variance function for variance-aware Michaelis-Menten
fitting.
}
\examples{
h_log <- variance_function("log")
h_sqrt <- variance_function("sqrt")
h_log(c(0, 5, 10))
h_sqrt(c(0, 5, 10))
}
