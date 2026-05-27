\name{alves_demo}
\alias{alves_demo}
\title{Subset of the Alves Soil Enzyme Kinetics Data}
\format{
A data frame with 247 rows and 6 variables:
\describe{
  \item{\code{core}}{Soil core identifier.}
  \item{\code{depth}}{Depth interval.}
  \item{\code{enzyme}}{Enzyme name.}
  \item{\code{temperature}}{Incubation temperature.}
  \item{\code{substrate_conc}}{Substrate concentration.}
  \item{\code{activity}}{Observed enzyme activity.}
}
}
\description{
A filtered subset of the Alves et al. soil exoenzyme kinetics data, restricted
to depth \code{00-10} and temperature \code{10} for fast reproducible package
examples.
}
\source{
Prepared from the Alves study files in the project workspace.
}
\examples{
head(alves_demo)
table(alves_demo$enzyme)
}
