\name{SecondDerivativePenalty}
\alias{SecondDerivativePenalty}
\title{A Roughness Penalty Based on the Squared Second Derivative}
\description{
A roughness penalty function based on squared second derivatives
evaluated numerically.  This is a possible template function for
other types of penalties. 
}
\usage{
SecondDerivativePenalty(xgrid, a)
}
\arguments{
\item{xgrid}{vector of length m, must be increasing}
\item{a}{a function of one numeric variable}
}
\value{
a vector of second derivatives evaluated at the points of xgrid
}
\author{W.J. Braun}

\keyword{models}
