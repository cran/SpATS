\name{variogram}
\alias{variogram}
\alias{variogram.SpATS}
\title{
Sample variogram
}
\description{
Computes the sample variogram from an \code{SpATS} object.
}
\usage{
\method{variogram}{SpATS}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{SpATS} as produced by \code{SpATS()}.}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented}
}
\details{
The present function computes the sample variogram on the basis of the (deviance) residuals of the fitted model. Currently, the function can only be applied for regular two-dimensional data, i.e, when the plots of the field are arranged in a regular two-dimensional array (usually defined by the column and row positions). 
For each pair of (deviance) residuals \eqn{e_i} and \eqn{e_j}, the half-squared difference is computed \deqn{v_{ij} = 0.5(e_i - e_j)^2,} as well as the corresponding column (\eqn{cd_{ij}}) and row displacements (\eqn{rd_{ij}}), with \deqn{cd_{ij} = |c_i - c_j|} and \deqn{rd_{ij} = |r_i - r_j|,} where \eqn{c_k} and \eqn{r_k} denote the column and row position of plot \eqn{k} respectively. The sample variogram is then defined as the triplet \deqn{(cd_{ij}, rd_{ij}, \bar{v}_{ij}),} where \eqn{\bar{v}_{ij}} denotes the average of the \eqn{v_{ij}} that share the same column and row displacements. 
For a more detailed description, see Gilmour et al. (1997).  
}
\value{
An object of class \code{variogram.SpATS} with the following components:
\item{data}{data frame including the following information: ``value'': the value of the sample variogram at each pair of column and row displacements; and ``length'': the number of observations used to compute the sample variogram at the corresponding pair of displacements.}
\item{col.displacement}{numerical vector containing the column displacements}
\item{row.displacement}{numerical vector containing the row displacements}
}
\references{
Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Stefanova, K.T., Smith, A.B. and Cullis, B.R. (2009). Enhanced Diagnostics for the Spatial Analysis of Field Trials. Journal of Agricultural, Biological, and Environmental Statistics, 14, 392 - 410.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpATS}}, \code{\link{plot.variogram.SpATS}}
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
