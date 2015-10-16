\name{MISE}
\alias{MISE}
\title{Approximate Mean Integrated Squared Error}
\description{
MISE for penalized sharpened regression based on trapezoid integration.
}
\usage{
MISE(x, xgrid, sigma2, lambda, h, g, A, B)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{numeric explanatory vector}
  \item{xgrid}{numeric vector}
  \item{sigma2}{numeric vector of variance(s)}
  \item{lambda}{numeric penalty constant}
  \item{h}{numeric bandwidth}
  \item{g}{regression function, numeric-valued}
  \item{A}{numeric matrix, smoother}
  \item{B}{numeric matrix, based on penalty}
}
\value{
  A vector containing the finite sample variance, squared bias,
and mean integrated squared error.
}
\author{W.J. Braun}

\keyword{models}
