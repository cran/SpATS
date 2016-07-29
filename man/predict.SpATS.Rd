\name{predict.SpATS}
\alias{predict.SpATS}
\title{
Predictions from an \code{SpATS} object 
}
\description{
Takes a fitted \code{SpATS} object produced by \code{SpATS()} and produces predictions.
}
\usage{
\method{predict}{SpATS}(object, newdata = NULL, which = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{SpATS} as produced by \code{SpATS()}.}
  \item{newdata}{an optional data frame to be used for obtaining the predictions.}
  \item{which}{an optional character string with the variables that define the margins of the multiway table to be predicted (see Details).}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\details{
This function allows to produce predictions, either specifying: (1) the data frame on which to obtain the predictions (argument \code{newdata}), or (2) those variables that define the margins of the multiway table to be predicted (argument \code{which}). In the first case, all fixed components (including genotype when fixed) and the spatial coordinates must be present in the data frame. As for the random effects is concerned, they are excluded from the predictions when the value is missing in the data frame. In the second case, predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. For those variables not specified in the argument \code{which}, the following rules have been considered: (a) random factors and the spatial trend are ignored in the predictions, (b) for fixed numeric variables, the mean value is considered; and (c) for fixed factors, the reference level is used.
}
\value{
The data frame used for obtaining the predictions, jointly with the predicted values and the corresponding standard errors. The label ``Excluded'' has been used to indicate those cases where a covariate has been excluded or ignored for the prediction (as for instance the random effect). 
}
\references{
Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{SpATS}}, \code{\link{obtain.spatialtrend}}
}
\examples{
library(SpATS)
data(wheatdata)
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

m0 <- SpATS(response = "yield", spatial = ~ SAP(col, row, nseg = c(10,20)), 
 genotype = "geno", fixed = ~ colcode + rowcode, random = ~ R + C, 
 data = wheatdata, control =  list(tolerance = 1e-03))

# Fitted values: prediction on the dataset used for fitting the model
pred1.m0 <- predict(m0, newdata = wheatdata)
pred1.m0[1:5,]

# Genotype prediction
pred2.m0 <- predict(m0, which = "geno")
pred2.m0[1:5,]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
