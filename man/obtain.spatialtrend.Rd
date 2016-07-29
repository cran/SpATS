\name{obtain.spatialtrend}
\alias{obtain.spatialtrend}
\title{
Predictions of the spatial trend from an \code{SpATS} object 
}
\description{
Takes a fitted \code{SpATS} object produced by \code{SpATS()} and produces predictions of the spatial trend on a regular two-dimensional array.
}
\usage{
obtain.spatialtrend(object, grid = c(100, 100), ...)
}
\arguments{
  \item{object}{an object of class \code{SpATS} as produced by \code{SpATS()}}
  \item{grid}{a numeric vector with the number of grid points along the x- and y- coordinates respectively. Atomic values are recycled. The default is 100.}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
For each spatial coordinate, \code{grid[k]} equally spaced values between the minimum and the maximum are computed (k = 1, 2). The spatial trend is then predicted on the regular two-dimensional array defined by each combination of the x- and y- coordinate values. 
}
\value{
A list with the following components:
\item{col.p}{x-coordinate values at which predictions have been computed.}
\item{row.p}{y-coordinate values at which predictions have been computed}
\item{fit}{a matrix of dimension length(row.p) x length(col.p) with the predicted spatial trend.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}, \code{\link{plot.SpATS}}, \code{\link{predict.SpATS}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, 
 data = wheatdata, control =  list(tolerance = 1e-03))

spat.trend.1 <- obtain.spatialtrend(m0)
spat.trend.2 <- obtain.spatialtrend(m0, grid = c(10, 10))

colors = topo.colors(100)
op <- par(mfrow = c(1,2))
fields::image.plot(spat.trend.1$col.p, spat.trend.1$row.p, t(spat.trend.1$fit), 
 main = "Prediction on a grid of 100 x 100", col = colors, xlab = "Columns", ylab = "Rows")
fields::image.plot(spat.trend.2$col.p, spat.trend.2$row.p, t(spat.trend.2$fit), 
 main = "Prediction on a grid of 10 x 10", col = colors, xlab = "Columns", ylab = "Rows")
par(op)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
