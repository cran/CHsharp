\name{sharp1d}
\alias{sharp1d}
\title{Data Sharpening for Density Estimation}
\description{
Application of Choi and Hall's (1999) data sharpening method for
univariate data, for use prior to density estimation.
}
\usage{
sharp1d(x, h, v = 1)
}
\arguments{
  \item{x}{the x coordinates of the data}
  \item{h}{the bandwidth for sharpening in the direction of the x
axis}
  \item{v}{a positive integer representing the number of iterations to  
perform}
}
\value{
Returns a vector containing the sharpened points x.sharp.
}
\references{
Choi, E. and Hall, P. (1999) Data sharpening as a prelude to density
estimation.  Biometrika 86, 941-947.
}
\author{Douglas G. Woolford, W. John Braun}

\examples{
# Example 1:
y <- c(rnorm(50,-1,1),rnorm(50,2,2), rnorm(100,0,.5))
data.sharp1 <- sharp1d(y,5,1)
data.sharp2 <- sharp1d(y,5,2)
# original data:
plot(density(y, bw=5))
# sharpened data after 1 iterations:
lines(density(data.sharp1, bw=5), col=2)
# sharpened data after 2 iterations:
lines(density(data.sharp2, bw=5), col=4)

x <- rt(100, df=3)
h <- dpik(x)

# Example 2:
curve(dt(x, df=3), from=-4, to=4)
lines(bkde(x, bandwidth=h), col=2, lty=2)

x.sharp <- sharp1d(x, h, 1)
lines(bkde(x.sharp, bandwidth=h), col=3, lty=3)

x.sharp2 <- sharp1d(x, h, 2)
lines(bkde(x.sharp2, bandwidth=h), col=4, lty=4)

x.sharp3 <- sharp1d(x, h, 3)
lines(bkde(x.sharp3, bandwidth=h), col=5, lty=5)
}
\keyword{models}
