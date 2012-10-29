\name{sharp3dB}
\alias{sharp3dB}
\title{Identify Cluster Centres for 3-dimensional Data via Data Sharpening}
\description{
Identifies the centres of clusters for 3-dimensional data using a 
convergent form of Choi and Hall's (1999) data sharpening method.  For 
use when the data is such that the z coordinates are in increasing 
order.}
\usage{
sharp3dB(x, y, z, hspace = 1, htime = 1, v = 1)
}
\arguments{
  \item{x}{the x coordinates of the data}
  \item{y}{the y coordinates of the data}
  \item{z}{the z coordinates of the data, in increasing order}
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
is the z-axis.  Provided the z-data is in increasing order, this function 
is significantly faster than sharp3d().}
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

\seealso{sharp3d}
\examples{
}
\keyword{cluster}
