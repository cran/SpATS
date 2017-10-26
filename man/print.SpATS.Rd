\name{print.SpATS}
\alias{print.SpATS}
\title{
Print method for \code{SpATS} objects
}
\description{
Default print method for objects fitted with \code{SpATS()} function.
}
\usage{
\method{print}{SpATS}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{SpATS} as produced by \code{SpATS()}}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
A short summary is printed including: the variable under study (response), the variable containing the genotypes, the spatial model, and the random and fixed components (when appropriate). Besides this information, the number of observations used to fit the model, as well as of those deleted due to missingness or zero weights are reported. Finally, the effective degrees of freedom (effective dimension) of the fitted model and the (REML) deviance is also printed.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
\code{\link{SpATS}}, \code{\link{summary.SpATS}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))

m0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
