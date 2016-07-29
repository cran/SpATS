\name{SpATS}
\alias{SpATS}
\title{
Spatial analysis of field trials with splines
}
\description{
Function specifically designed for the analysis of field trials experiments with the spatial trend being modelled by means of two-dimensional P-splines.
}
\usage{
SpATS(response, genotype, geno.decomp = NULL, genotype.as.random = FALSE, spatial, 
 fixed = NULL, random = NULL, data, family = gaussian(), offset = 0, weights = NULL, 
 control = controlSpATS())
}
\arguments{
  \item{response}{a character string with the name of the variable that contains the response variable of interest.}
  \item{genotype}{a character string with the name of the variable that contains the genotypes or varieties. This variable must be a factor on the data frame.}
  \item{geno.decomp}{an optional character string with the name of a factor variable according which genotypes are grouped, with different genetic variance being assumed for each group. Only applies when genotype is random}
  \item{genotype.as.random}{logical. If TRUE, the genotype is included as random effect in the model. The default is FALSE.}
  \item{spatial}{a right hand \code{\link{formula}} object specifying the spatial P-Spline model. See \code{\link{SAP}} and \code{\link{PSANOVA}} for more details about how to specify the spatial trend.}
  \item{fixed}{an optional right hand \code{\link{formula}} object specifying the fixed effects.}
  \item{random}{an optional right hand \code{\link{formula}} object specifying the random effects. Currently, only sets of independent and identically distributed random effects can be incorporated.}
  \item{data}{data frame containing all needed variables}
  \item{family}{object of class \code{\link{family}} specifying the distribution and link function.}
  \item{offset}{an optional numerical vector containing an a priori known component to be included in the linear predictor during fitting. }
  \item{weights}{an optional numerical vector of weights to be used in the fitting process. By default, the weights are considered to be one.}
  \item{control}{a list of control values to replace the default values returned by the function \code{\link{controlSpATS}}.}
}
\details{
This function fits a (generalised) linear mixed model, with the variance components being estimated by the restricted log-likelihood (REML). The estimation procedure is an extension of the SAP (Separation of anisotropic penalties) algorithm by Rodriguez - Alvarez et al. (2015). In this package, besides the spatial trend, modelled as a two-dimensional P-spline, the implemented estimation algorithm also allows for the incorporation of both fixed and (sets of i.i.d) random effects on the (generalised) linear mixed model.
As far as the implementation of the algorithm is concerned, the sparse structure of the design matrix associated with the genotype has been taken into account. The combination of both the SAP algorithm and this sparse structure makes the algorithm computational efficient, allowing the proposal to be applied for the analysis of large datasets.  
}
\value{
A list with the following components:
\item{call}{the matched call.}
\item{data}{the original supplied data argument with a new column with the weights used during the fitting process.}
\item{model}{a list with the model components: response, spatial, genotype, fixed and/or random.}
\item{fitted}{a numeric vector with the fitted values.}
\item{residuals}{a numeric vector with deviance residuals.}
\item{psi}{a two-length vector with the values of the dispersion parameter at convergence. For Gaussian responses both elements coincide, being the (REML) estimate of dispersion parameter. For non-Gaussian responses, the result depends on the argument \code{update.psi} of the \code{\link{controlSpATS}} function. If this argument was specified to \code{FALSE} (the default), the first component of the vector corresponds to the default value used for the dispersion parameter (usually 1). The second element, correspond to the (REML) estimate of the dispersion parameter at convergence. If the argument \code{update.psi} was specified to \code{TRUE}, both components coincide (as in the Gaussian case).}
\item{var.comp}{a numeric vector with the (REML) variance component estimates. These vector contains the variance components associated with the spatial trend, as well as those related with the random model terms.}
\item{eff.dim}{a numeric vector with the estimated effective dimension (or effective degrees of freedom) for each model component (genotype, spatial, fixed and/or random)}
\item{dim}{a numeric vector with the (model) dimension of each model component (genotype, spatial, fixed and/or random). This value corresponds to the number of parameters to be estimated}
\item{nobs}{number of observations used to fit the model}
\item{niterations}{number of iterations EM-algorithm}
\item{deviance}{the (REML) deviance at convergence (i.e., - 0.5 times the restricted log-likelihood)}
\item{coeff}{a numeric vector with the estimated fixed and random effect coefficients.}
\item{terms}{a list with the model terms: response, spatial, genotype, fixed and/or random. The information provided here is useful for printing and prediction purposes.}
\item{vcov}{inverse of the coefficient matrix of the mixed models equations. The inverse is needed for the computation of standard errors. For computational issues, the inverse is returned as a list: C11_inv corresponds to the coefficient matrix associated with the genotype; C22_inv corresponds to the coefficient matrix associated with the spatial, the fixed and the random components; and C12_inv and C21_inv correspond to the ``combination'' of both.}
}
\references{
Currie, I., and Durban, M. (2002). Flexible smoothing with P-splines: a unified approach. Statistical Modelling, 4, 333 - 349.

Eilers, P.H.C., and Marx, B.D. (2003). Multivariate calibration with temperature interaction using two-dimensional penalised signal regression. Chemometrics and Intelligent Laboratory Systems, 66, 159 - 174.

Eilers, P.H.C., Marx, B.D., Durban, M. (2015). Twenty years of P-splines. SORT, 39(2), 149 - 186.

Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2016). Spatial Models for Field Trials. Available on arXiv at http://arxiv.org/abs/1607.08255.

Rodriguez-Alvarez, M.X., Lee, D.-J., Kneib, T., Durban, M., and Eilers, P.H.C. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. Statistics and Computing, 25, 941 - 957.

}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS-package}}, \code{\link{controlSpATS}}, \code{\link{SAP}}, \code{\link{PSANOVA}}
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

# Plot results
plot(m0)
plot(m0, all.in.one = FALSE)

# Variogram
var.m0 <- variogram(m0)
plot(var.m0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
