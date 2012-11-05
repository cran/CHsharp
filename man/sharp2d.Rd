\name{sharp2d}
\alias{sharp2d}
\title{Identify Cluster Centres for 2-dimensional Data via Data Sharpening}
\description{
Identifies the centres of clusters for 2-dimensional data using a
converged form of Choi and Hall's (1999) data sharpening method.
}
\usage{
sharp2d(x, y, hspace = 1, htime = 1, v = 1)
}
\arguments{
  \item{x}{the x coordinates of the data}
  \item{y}{the y coordinates of the data}
  \item{hspace}{the bandwidth for sharpening in the direction of the x
axis}
  \item{htime}{the bandwidth for sharpening in the y direction}
  \item{v}{a positive integer representing the number of iterations to  
perform}
}
\details{
Identifies the centres of clusters based on a converged form of Choi
and Hall's data sharpening method. This function was originally built for
identifying clusters in space-time where space is the x-y plane and time
is the z-axis.
}
\value{
Returns a (number of data points x 2) data frame containing the sharpened
points x.sharp and y.sharp, respectively.
}
\references{
Woolford, D. G. and Braun, W. J. (2004) Exploring lightning
and fire ignition data as point processes.  2004 Proceeding of the
American Statistical Association, Statistics and the Environment Section
[CD-ROM], Alexandria, VA: American Statistical Association.

Choi, E. and Hall, P. (1999) Data sharpening as a prelude to density
estimation.  Biometrika 86, 941-947.
}
\author{Douglas G. Woolford, W. John Braun}

\examples{
x <- 1:200
y <- c(rnorm(50,-1,1),rnorm(50,2,2), rnorm(100,0,.5))
data.sharp5 <- sharp2d(x,y,5,10,5)
data.sharp10 <- sharp2d(x,y,5,10,10)
# original data:
plot(x,y)
# sharpened data after 5 iterations:
points(data.sharp5$x.sharp, data.sharp5$y.sharp, col=2,pch=19)
# sharpened data after 10 iterations:
points(data.sharp10$x.sharp, data.sharp10$y.sharp, col=4, pch=19)
}
\keyword{cluster}
