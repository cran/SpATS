\name{summary.SpATS}
\alias{summary.SpATS}
\title{
Summary method for \code{SpATS} objects
}
\description{
Default summary method for objects fitted with \code{SpATS()} function.
}
\usage{
\method{summary}{SpATS}(object, which = c("dimensions", "variances", "all"), ...)
}
\arguments{
  \item{object}{an object of class \code{SpATS} as produced by \code{SpATS()}}
  \item{which}{character vector indicating the information to be shown on the screen. ``dimensions'': prints the dimensions associated to each component of the fitted model (see details), ``variances'': prints the variance components and standard deviation (SD), the residual variance, and the ratio between the residual variance and the variance components (in log10 scale); and, ``all'': both dimensions and variances are printed. }
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
This function provides relevant information about the fitted model. The dimensions associated to each model term, as well as the variance components (jointly with the standard deviation), can be obtained. Specifically, the following information is shown:
\describe{ 
	\item{dimensions}{dimensions associated to each model term (genotype, spatial, fixed and/or random):
		\itemize{ 
		\item Effective: the effective dimension or effective degrees of freedom. The concept of effective dimension is well known in the smoothing framework. However, its use has been less recognized in the mixed model framework. The effective dimension of a fitted model is defined as the trace of the so-called hat matrix (Hastie and Tibshirani, 1990). For each model term, the effective dimension is then defined as the trace of the block of the hat matrix corresponding to this term. For fixed effects, the effective dimension coincides with the model dimension (see below). For random and smooth terms, however, the effective dimension is usually lower than the model dimension. In these cases, the effective dimension or effective degrees of freedom can be interpreted as the effective number of estimated parameters. It is worth nothing that for the specification of the spatial trend using the \code{\link{SAP}} function, the effective dimension associated to each spatial coordinates is shown (see Rodriguez-Alvarez et al. (2015) for details), as well as the global effective dimension (the sum of these two values). For the ANOVA-type decomposition using \code{\link{SAP}}, besides the former information, the effective dimension associated to each main effect is also shown. This is in concordance with the number of smoothing parameters (or variance components in the mixed model framework), used to model the spatial trend. Accordingly, for the spatial trend specified using the function \code{\link{PSANOVA}}, the effective dimension of each of the five components is shown (and no global summary is reported). More information can be found in Rodriguez-Alvarez et al. (2016).  
		\item Model: the model dimension corresponds to the number of parameters to be estimated. For the specification of the spatial trend using the \code{\link{SpATS}} function, the model and nominal dimensions are only specified for the global term (and for the main effects when using the ANOVA-type decomposition).
		\item Nominal: for the random terms of the model, the nominal dimension corresponds to the model dimension minus one, ``lost'' due to the implicit constraint that ensures the mean of the random effects to be zero.
		\item Ratio: ratio between the nominal dimension and the model dimension. For the genotype (when random), this ratio corresponds to the generalized heritability proposed by Oakey (2006). A deeper discussion can be found in Rodriguez - Alvarez et al. (2016) 
		\item Type: Model term type 'F' Fixed, 'R' Random, and 'S' Smooth (spatial trend). 
	}}
\item{variances}{variance components associated to the random terms and the smooth spatial trend, as well as the ratio between the residual variance and the variance components (in log10 scale). In the smoothing framework, this ratio corresponds to the smoothing parameter.}
}
}
\value{
Returns an object of class \code{summary.SpATS} with the same components as an \code{SpATS} object (see \code{\link{SpATS}}) plus:
\item{p.table.dim}{a matrix containing all the information related to the dimensions of the fitted model.}
\item{p.table.vc}{a matrix containing all the information related to the variance components of the fitted model}
}
\references{
Hastie, T., and Tibshirani, R. (1990). Generalized Additive Models. In: Monographs on Statistics and Applied Probability, Chapman and Hall, London.

Oakey, H., A. Verbyla, W. Pitchford, B. Cullis, and H. Kuchel (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113, 809 - 819.

Rodriguez-Alvarez, M.X., Boer, M.P., van Eeuwijk, F.A, and Eilers, P.H.C. (2016). Spatial Models for Field Trials. Available on arXiv at http://arxiv.org/abs/1607.08255.

Rodriguez-Alvarez, M.X., Lee, D.-J., Kneib, T., Durban, M., and Eilers, P.H.C. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. Statistics and Computing, 25, 941 - 957.
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
summary(wheatdata)

# Create factor variable for row and columns
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))
# Brief summary
m0
# More information: dimensions
summary(m0) # summary(fit.m2, which = "dimensions") 
# More information: variances
summary(m0, which = "variances") 
# More information: all
summary(m0, which = "all") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
