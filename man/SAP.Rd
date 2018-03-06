\name{SAP}
\alias{SAP}
\title{
Define a two-dimensional penalised tensor-product of marginal B-Spline basis functions.
}
\description{
Auxiliary function used for modelling the spatial or environmental effect as a two-dimensional penalised tensor-product of marginal B-spline basis functions with anisotropic penalties. The approach implemented here corresponds to the Separation of Anisotropic Penalties (SAP) algorithm introduced by Rodriguez-Alvarez et al. (2015). 
}
\usage{
SAP(..., nseg = c(10,10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1), 
ANOVA = FALSE)
}
\arguments{
  \item{\dots}{a list of the covariates that this smooth component is a function of. Currently, only two-dimensional tensor-product smoothers are implemented: the first covariate is assumed to define the x-spatial coordinate (e.g., column position) of each plot in the field, and the second argument the y-spatial coordinate (i.e., the row position). Both covariates should be numerical.}
  \item{nseg}{numerical vector of length 2 containing the number of segments for each marginal (strictly \code{nseg} - 1 is the number of internal knots in the domain of the covariate). Atomic values are also valid, being recycled. Default set to 10.}
  \item{pord}{numerical vector of length 2 containing the penalty order for each marginal. Atomic values are also valid, being recycled. Default set to 2 (second order).}
  \item{degree}{numerical vector of length 2 containing the order of the polynomial of the B-spline basis for each marginal. Atomic values are also valid, being recycled. Default set to 3 (cubic B-splines).}
  \item{nest.div}{numerical vector of length 2 containing the divisor of the number of segments (\code{nseg}) to be used for the construction of the nested B-spline basis for the smooth-by-smooth interaction component. In this case, the nested B-spline basis will be constructed assuming a total of \code{nseg}/\code{nest.div} segments. Default set to 1, which implies that nested basis are not used. See details.}
  \item{ANOVA}{logical. If TRUE, four different smoothing parameters are considered: two for the main effects and two for the smooth interaction. See details.}
}
\details{
Within the P-spline framework, anisotropic low-rank tensor-product smoothers have become the general approach for modelling multidimensional surfaces (Eilers and Marx 2003; Wood 2006). In this package, we propose to model the spatial or environmental effect by means of the tensor-product of B-splines basis functions. In other words, we propose to model the spatial trend as a smooth bivariate surface jointly defined over the the spatial coordinates. Accordingly, the current function has been designed to allow the user to specify the spatial coordinates that the spatial trend is a function of. There is no restriction about how the spatial coordinates shall be specified: these can be the longitude and latitude of the position of the plot on the field or the column and row numbers. The only restriction is that the variables defining the spatial coordinates should be numeric (in contrast to factors).

As far as estimation is concerned, we have used in this package the equivalence between P-splines and linear mixed models (Currie and Durban, 2002). Under this approach, the smoothing parameters are expressed as the ratio between variance components. Moreover, the smooth components are decomposed in two parts: one which is not penalised (and treated as fixed) and one with is penalised (and treated as random). For the two-dimensional case, the mixed model representation leads also to a very interesting decomposition of the penalised part of the bivariate surface in three different components (Lee and Durban, 2011): (a) a component that contains the smooth main effect (smooth trend)	along one of the covariates that the surface is a function of (as, e.g, the x-spatial coordinate or column position of the plot in the field), (b) a component that contains the smooth main effect (smooth trend) along the other covariate (i.e., the y-spatial coordinate or row position); and (c) a smooth interaction component (sum of the linear-by-smooth interaction components and the smooth-by-smooth interaction component). 
The default implementation used in this function (\code{ANOVA} = \code{FALSE}) assumes two different smoothing parameters, i.e., one for each covariate in the smooth component. Accordingly, the same smoothing parameters are used for both, the main effects and the smooth interaction. However, this approach can be extended to deal with the ANOVA-type decomposition presented in Lee and Durban (2011). In their approach, four different smoothing parameters are considered for the smooth surface, that are in concordance with the aforementioned decomposition: (a) two smoothing parameter, one for each of the main effects; and (b) two smoothing parameter for the smooth interaction component. Here, this decomposition can be obtained by specifying the argument \code{ANOVA = TRUE}.  

It should be noted that, the computational burden associated with the estimation of the two-dimensional tensor-product smoother might be prohibitive if the dimension of the marginal bases is large. In these cases, Lee et al. (2013) propose to reduce the computational cost by using nested bases. The idea is to reduce the dimension of the marginal bases (and therefore the associated number of parameters to be estimated), but only for the smooth-by-smooth interaction component. As pointed out by the authors, this simplification can be justified by the fact that the main effects would in fact explain most of the structure (or spatial trend) presented in the data, and so a less rich representation of the smooth-by-smooth interaction component could be needed. In order to ensure that the reduced bivariate surface is in fact nested to the model including only the main effects, Lee et al. (2013) show that the number of segments used for the nested basis should be a divisor of the number of segments used in
the original basis (\code{nseg} argument). In the present function, the divisor of the number of segments is specified through the argument \code{nest.div}. For a more detailed review on this topic, see Lee (2010) and Lee et al. (2013).
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Currie, I., and Durban, M. (2002). Flexible smoothing with P-splines: a unified approach. Statistical Modelling, 4, 333 - 349.

Eilers, P.H.C., and Marx, B.D. (2003). Multivariate calibration with temperature interaction using two-dimensional penalised signal regression. Chemometrics and Intelligent Laboratory Systems, 66, 159 - 174.

Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Lee, D.-J., 2010. Smothing mixed model for spatial and spatio-temporal data. Ph.D. Thesis, Department of Statistics, Universidad Carlos III de Madrid, Spain.

Lee, D.-J., and Durban, M. (2011). P-spline ANOVA-type interaction models for spatio-temporal smoothing. Statistical Modelling 11, 49 - 69.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Rodriguez-Alvarez, M.X., Lee, D.-J., Kneib, T., Durban, M., and Eilers, P.H.C. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. Statistics and Computing, 25, 941 - 957.

Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for generalized additive models. Journal of the Royal Statistical Society: Series B, 70, 495 - 518.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}, \code{\link{PSANOVA}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

# Original anisotropic approach: two smoothing parameters
m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, 
 data = wheatdata, control =  list(tolerance = 1e-03))
summary(m0)

# ANOVA decomposition: four smoothing parameters
\dontrun{
m1 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20), ANOVA = TRUE), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
 control =  list(tolerance = 1e-03))
summary(m1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
