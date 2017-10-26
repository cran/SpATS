\name{PSANOVA}
\alias{PSANOVA}
\title{
Define a two-dimensional penalised tensor-product of marginal B-Spline basis functions based on the P-spline ANOVA (PSANOVA) approach.
}
\description{
Auxiliary function used for modelling the spatial or environmental effect as a two-dimensional penalised tensor-product of marginal B-spline basis functions with anisotropic penalties on the basis of the PSANOVA approach by Lee et al. (2013).
}
\usage{
PSANOVA(..., nseg = c(10,10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1))
}
\arguments{
  \item{\dots}{a list of the covariates that this smooth component is a function of. Currently, only two-dimensional tensor-product smoothers are implemented: the first covariate is assumed to define the x-spatial coordinate (e.g., column position) of each plot in the field, and the second argument the y-spatial coordinate (i.e., the row position). Both covariates should be numerical.}
  \item{nseg}{numerical vector of length 2 containing the number of segments for each marginal (strictly \code{nseg} - 1 is the number of internal knots in the domain of the covariate). Atomic values are also valid, being recycled. Default set to 10.}
  \item{pord}{numerical vector of length 2 containing the penalty order for each marginal. Atomic values are also valid, being recycled. Default set to 2 (second order). Currently, only second order penalties are allowed.}
  \item{degree}{numerical vector of length 2 containing the order of the polynomial of the B-spline basis for each marginal. Atomic values are also valid, being recycled. Default set to 3 (cubic B-splines).}
  \item{nest.div}{numerical vector of length 2 containing the divisor of the number of segments (\code{nseg}) to be used for the construction of the nested B-spline basis for the smooth-by-smooth interaction component. In this case, the nested B-spline basis will be constructed assuming a total of \code{nseg}/\code{nest.div} segments. Default set to 1, which implies that nested basis are not used. See \code{\link{SAP}} for more details.}
}
\details{
The approach implemented here represents an alternative method to the \code{\link{SAP}} function. In this case, the smooth bivariate surface (or spatial trend) is decomposed in five different components each of them depending on a single smoothing parameter (see Lee et al., 2013).
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

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}, \code{\link{SAP}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

# Without nested basis
m0 <- SpATS(response = "yield", spatial = ~ PSANOVA(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))

summary(m0)

# With nested basis
m1 <- SpATS(response = "yield", spatial = ~ PSANOVA(col, row, nseg = c(10,20), nest.div = 2), 
genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
control =  list(tolerance = 1e-03))

summary(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
