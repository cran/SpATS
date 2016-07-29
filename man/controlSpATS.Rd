\name{controlSpATS}
\alias{controlSpATS}
\title{
Used to set various parameters controlling the fitting process
}
\description{
This function can be used to modify some default parameters that control the estimation of an SpATS model.
}
\usage{
controlSpATS(maxit = 200, tolerance = 0.001, monitoring = 2, update.psi = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxit}{numerical value indicating the maximum number of iterations. Default set to 200 (see Details).}
  \item{tolerance}{numerical value indicating the tolerance for the convergence criterion. Default set to 0.001 (see Details).}
  \item{monitoring}{numerical value determining the level of printing which is done during the estimation. The value of 0 means that no printing is produced, a value of 1 means that only the computing times are printed, and a value of 2 means that, at each iteration, the (REML) deviance and the effective dimensions of the random components are printed. Default set to 2.}
  \item{update.psi}{logical. If TRUE, the dispersion parameter of the exponential family is updated at each iteration of the estimation algorithm. Default is FALSE, except for \code{family = gaussian()}, where the dispersion parameter (or residual variance), is always updated jointly with the rest of variance components.}
}
\details{
The estimation procedure implemented in the SpATS package is an extension of the SAP (Separation of anisotropic penalties) algorithm by Rodriguez - Alvarez et al. (2015). In this case, besides the spatial trend, modelled as a two-dimensional P-spline, the estimation algorithm allows for the incorporation of both fixed and (sets of i.i.d) random effects on the (generalised) linear mixed model. 
For Gaussian response variables, the algorithm is an iterative procedure, with the fixed and random effects as well as the variance components being updated at each iteration. To check the convergence of this iterative procedure, the (REML) deviance is monitored.
For non-Gaussian response variables, estimation is based on Penalized Quasi-likelihood (PQL) methods. Here, the algorithm is a two-loop algorithm: the outer loop corresponds to the Fisher-Scoring algorithm (monitored on the basis of the change in the linear predictor between consecutive iterations), and the inner loop corresponds to that described for the Gaussian case. 
}
\value{
a list with components for each of the possible arguments.
}
\references{
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
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

# Default control parameters
m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, 
 data = wheatdata)

# Modified the number of iterations, the tolerance, and the monitoring
m1 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, 
 data = wheatdata, control = list(maxit = 50, tolerance = 1e-06, monitoring = 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
