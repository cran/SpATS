\name{wheatdata}
\alias{wheatdata}
\docType{data}
\title{
Wheat yield in South Australia
}
\description{
A randomized complete block experiment of wheat in South Australia (Gilmour et al., 1997).
}
\usage{data("wheatdata")}
\format{
  A data frame with 330 observations on the following 7 variables
  \describe{
    \item{\code{yield}}{
    	yield, numeric}
    \item{\code{geno}}{
    	wheat variety, a factor with 107 levels}
    \item{\code{rep}}{
    	replicate factor, 3 levels}
    \item{\code{row}}{
    	row position, numeric}
    \item{\code{col}}{
    	column position, a numeric vector}
    \item{\code{rowcode}}{
    	a factor with 2 levels. This factor is related to the way the trial and the plots were sown (see Details).}
    \item{\code{colcode}}{
    	a factor with 3 levels. This factor is related to the way the plots were trimmed to an assumed equal length before harvest (see Details).}
  }
}
\details{
A near complete block design experiment. In this experiment, a total of 107 varieties were sown in 3 replicates in a near complete block design (3 varieties were sown twice in each replicate). Each replicate consisted of 5 columns and 22 rows. Accordingly, the field consists on 15 columns and 22 rows (i.e. a total of 330 plots or observations). Trimming was done by spraying the wheat with herbicide. The sprayer travelled in a serpentine pattern up and down columns. More specifically, the plots were trimmed according to four possible combinations: up/down, down/down, down/up, and up/up. The \code{colcode} factor describe this 4-phase sequence. With respect to the way the trial was sown, the procedure was done in a serpentine manner with a planter that seeded three rows at a time (Left, Middle, Right). This fact led to a systematic effect on the yield due to sowing. The variable \code{rowcode} describes the position of each plot in the triplet (Left, Middle, Right). Level 1 was used for the middle plot, and 2 for those on the right or on the left of the triplet (see Gilmour et al. (1997) for more details).
}
\source{
Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

Kevin Wright (2015). agridat: Agricultural datasets. R package version 1.12. http://CRAN.R-project.org/package=agridat
}
\examples{
library(SpATS)
data(wheatdata)
summary(wheatdata)
}
\keyword{datasets}
