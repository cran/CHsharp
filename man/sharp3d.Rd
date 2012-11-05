\name{sharp3d}
\alias{sharp3d}
\title{Identify Cluster Centres for 3-dimensional Data via Data Sharpening}
\description{
Identifies the centres of clusters for 3-dimensional data using a
convergent form of Choi and Hall's (1999) data sharpening method.
}
\usage{
sharp3d(x, y, z, hspace = 1, htime = 1, v = 1)
}
\arguments{
  \item{x}{the x coordinates of the data}
  \item{y}{the y coordinates of the data}
  \item{z}{the z coordinates of the data}
  \item{hspace}{the bandwidth for sharpening in the direction of the x-y
plane}
  \item{htime}{the bandwidth for sharpening in the z direction}
  \item{v}{a positive integer representing the number of iterations to  
perform}
}
\details{
Identifies the centres of clusters based on a convergent form of Choi
and Hall's data sharpening method. This function was originally built for
identifying clusters in space-time where space is the x-y plane and time
is the z-axis.
}
\value{
Returns a (number of data points x 3) data frame containing the sharpened
points x.sharp, y.sharp and z.sharp, respectively.
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

\seealso{sharp3dB}
\examples{
x <- 1:200
y <- c(rnorm(50,-1,1),rnorm(50,2,2), rnorm(100,0,.5))
z <- c(sample(1:50,50), sample(26:75,50), sample(51:150,100))
data.sharp5 <- sharp3d(x,y,z,5,10,5)
data.sharp10 <- sharp3d(x,y,z,5,10,10)
# original data:
dataPlot <- scatterplot3d(x,y,z)
# sharpened data after 5 iterations:
dataPlot$points3d(data.sharp5$x.sharp, data.sharp5$y.sharp, 
data.sharp5$z.sharp, col=2,pch=19)
# sharpened data after 10 iterations:
dataPlot$points3d(data.sharp10$x.sharp, data.sharp10$y.sharp, 
data.sharp10$z.sharp, col=4, pch=19)
}
\keyword{cluster}
