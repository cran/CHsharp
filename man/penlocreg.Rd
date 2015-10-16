\name{penlocreg}
\alias{penlocreg}
\title{Penalized Local Polynomial Regression}
\description{
Data sharpened local polynomial regression subject to a given penalty.
}
\usage{
penlocreg(x, y, xgrid, degree = 0, h, lambda, L, ...)
}
\arguments{
\item{x}{numeric vector of predictor observations}
\item{y}{numeric vector of observed responses}
\item{xgrid}{numeric vector of grid points where regression function is evaluated}
\item{degree}{numeric vector of local polynomial regression degree}
\item{h}{numeric bandwidth}
\item{lambda}{numeric penalty constant}
\item{L}{function related to penalty}
\item{...}{additional arguments, as required by L}
}
\value{
a list containing the original observed predictor values, the
sharpened responses, the smoother matrix and the penalty matrix
}
\examples{ 
xx <- faithful$waiting 
yy <- faithful$eruptions 
h <- dpill(xx,yy)/2; lam <- 20  # tuning parameter selections 
yy.pen <- penlocreg(xx, yy, seq(min(xx), max(xx), len=401), lambda=lam, degree=1, h = h, L = 
    SecondDerivativePenalty) 
plot(xx, yy, xlab="waiting", ylab="eruptions", col="grey") 
title("Old Faithful")
points(yy.pen, col=2, cex=.6) # sharpened data points
lines(locpoly(xx, yy, bandwidth=h*2, degree=1), lwd=2) # local linear estimate
lines(locpoly(yy.pen$x, yy.pen$y, bandwidth=h, degree=1), col=2, lwd=2) # sharpened estimate 
}
\keyword{models}


