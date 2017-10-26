\name{plot.variogram.SpATS}
\alias{plot.variogram.SpATS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default variogram.SpATS plotting
}
\description{
Takes a fitted \code{variogram.SpATS} object produced by \code{variogram.SpATS()} and plots the associated sample variogram using an RGL 3D perspective plot (package \code{plot3Drgl}).
}
\usage{
\method{plot}{variogram.SpATS}(x, min.length = 30, ...)
}
\arguments{
  \item{x}{an object of class \code{variogram.SpATS} as produced by \code{variogram.SpATS()}.}
  \item{min.length}{numerical value. The sample variogram is depicted including only those pairs with more than \code{min.length} observations (see \code{\link{variogram.SpATS}}).}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
This function as well as function \code{variogram.SpATS()} can only be used for regular two dimensional data.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Stefanova, K.T., Smith, A.B., and Cullis, B.R. (2009). Enhanced Diagnostics for the Spatial Analysis of Field Trials. Journal of Agricultural, Biological, and Environmental Statistics, 14, 392 - 410.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}, \code{\link{variogram.SpATS}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))

# Compute the variogram
var.m0 <- variogram(m0)
# Plot the variogram
plot(var.m0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
