\name{lambda}
\alias{lambda}
\title{Penalty Parameter Selector}
\description{
Data-driven selector of the penalty parameter, given a bandwidth.
}
\usage{lambda(x, y, h, d, xgrid, A, B, niterations=2) 
}
\arguments{ \item{x}{numeric vector of predictor observations} 
\item{y}{numeric vector of observed responses} 
\item{h}{numeric bandwidth}
\item{d}{numeric degree of local polynomial regression}
\item{xgrid}{numeric vector of grid points where regression function is to be evaluated}
\item{A}{numeric matrix, Smoother matrix} 
\item{B}{numeric matrix, based on penalty}
\item{niterations}{number of iterations}
}
\value{
a numeric vector of smoothing parameters, corresponding to successive iterates
}
\author{W.J. Braun}
\keyword{models}
