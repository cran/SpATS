\name{getHeritability}
\alias{getHeritability}
\title{
Calculate heritabilities from \code{SpATS} objects
}
\description{
For the genotype (when random), the function returns the generalized heritability proposed by Oakey (2006).
}
\usage{
getHeritability(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{SpATS} as produced by \code{SpATS()}}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
A numeric vector (usually of length 1) with the heritabilities.
}
\references{
Oakey, H., A. Verbyla, W. Pitchford, B. Cullis, and H. Kuchel (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113, 809 - 819.

Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics, 23, 52 - 71. https://doi.org/10.1016/j.spasta.2017.10.003.
}
\seealso{
\code{\link{SpATS}}, \code{\link{summary.SpATS}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
 genotype = "geno", genotype.as.random = TRUE,
 fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))

getHeritability(m0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
