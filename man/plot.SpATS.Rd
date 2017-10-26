\name{plot.SpATS}
\alias{plot.SpATS}
\title{
Default SpATS plotting
}
\description{
Takes a fitted \code{SpATS} object produced by \code{SpATS()} and plots six different graphics (see Details).
}
\usage{
\method{plot}{SpATS}(x, all.in.one = TRUE, main = NULL, annotated = FALSE, depict.missing = FALSE, ...)
}
\arguments{
	\item{x}{an object of class \code{SpATS} as produced by \code{SpATS()}.}
	\item{all.in.one}{logical. If TRUE, the four plots are depicted in one window. Default is TRUE}
	\item{main}{character string specifying the main title to appear on the plot. By default (i.e. when main = NULL), the variable under study is incorporated in the title of the plot.}
	\item{annotated}{logical. If TRUE, the variable under study and the models used is added to the plot. Only applied when argument \code{all.in.one} is TRUE.}
	\item{depict.missing}{logical. If TRUE, the estimated spatial trend is depicted for all plots in the field, even for those with missing values.}
	\item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
The following graphics are depicted: the raw data, the fitted data (on the response scale), the (deviance) residuals, the estimated spatial trend (excluding the intercept), the genotypic BLUEs (or BLUPs) and their histogram. Except for the histogram, the plots are depicted in terms of the spatial coordinates (e.g., the rows and columns of the field).  
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2017). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics (to appear). https://doi.org/10.1016/j.spasta.2017.10.003.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))

# Default plotting
plot(m0)
# Annotated
plot(m0, annotated = TRUE, main = "Wheat data (Gilmour et al., 1997)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
