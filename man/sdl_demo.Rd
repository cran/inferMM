\name{sdl_demo}
\alias{sdl_demo}
\title{Self-Driving Laboratory Michaelis-Menten Demo Data}
\format{
A data frame with 168 rows and 5 variables:
\describe{
  \item{\code{enzyme}}{Enzyme identifier.}
  \item{\code{s_uM}}{Substrate concentration in micromolar units.}
  \item{\code{replicate}}{Replicate index.}
  \item{\code{v_uM_per_s}}{Velocity in micromolar units per second.}
  \item{\code{v_uM_per_min}}{Velocity in micromolar units per minute.}
}
}
\description{
A lightweight self-driving-laboratory enzyme-kinetic dataset used in the
package examples and vignette.
}
\source{
Prepared from the SDL study files in the project workspace.
}
\examples{
head(sdl_demo)
table(sdl_demo$enzyme)
}
