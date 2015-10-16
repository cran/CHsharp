\name{sharpen}
\alias{sharpen}
\title{Penalized Data Sharpening Operator for Local Polynomial Regression}
\description{
Data perturbation operator which moves responses a minimal amount subject
to a given penalty.
}
\usage{
sharpen(x, y, lambda, B) 
}
\arguments{
\item{x}{numeric vector of predictor observations}
\item{y}{numeric vector of observed responses}
\item{lambda}{numeric penalty constant}
\item{B}{numeric matrix, based on penalty}
}
\value{
a numeric vector containing the sharpened responses
}
\author{W.J. Braun}
\keyword{models}
