\name{PH2XBAR}
\alias{PH2XBAR}
\title{PH2XBAR}
\usage{
PH2XBAR(
    X
    ,PH1.info = list(
                  X = NULL
                  , mu = NULL
                  , sigma = NULL
                  , k = NULL
                  , n = NULL
                  , model = "ANOVA-based"
                )
    ,c.ii.info = list(
                  c.ii = NULL
                  , method = 'UC'
                  , ARL0 = 370
                  , p = 0.05
                  , eps = 0.1
                  , interval.c.ii.UC = c(1, 3.2)
                  , interval.c.ii.EPC = c(1, 10)
                  , UC.tol = .Machine$double.eps^0.25
                  , EPC.tol = .Machine$double.eps^0.25
                )
    #,alternative = '2-sided'
    ,c4.option = TRUE
    ,plot.option = TRUE
    ,maxsim = 100000
)
}
\arguments{
  \item{X}{ Phase II sample }
  \item{PH1.info}{ the information of Phase I chart. X - Phase I sample, and X can be null if other parameters are filled. mu - mean of Phase I sample. sigma - standard deviation of Phase I sample. k - Phase I batch size. n - Phase I subgroup size. model - 'ANOVA-based' - degrees of freedom equal to k - 1 and estimating the within-batch standard deviation; 'basic' - degree of freedo equal to k * (n - 1) and estimating the pooled standard deviation }
  \item{c.ii.info}{ the information of Phase II charting constant. c.ii - a nominal Phase II charting constant. If this is given, the function will not recompute the Phase II charting constant. method - 'UC' - using the unconditional method; 'EPC' - using the exceedance probability criterion; 'both' - using both of the unconditional method and the exceedance probability criterion. ARL0 - a nominal in-control average run length. p - exceedance probability. eps - a nominal allowable changes of the in-control average run length. interval.c.ii.UC - the searching range of charting constants based on the unconditional method. interval.c.ii.EPC - the searching range of charting constants based on the exceedance probability criterion. UC.tol - the tolerance level for UC. EPC.tol - the tolerance level for EPC}
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{maxsim}{ maximum simulations of computation }
}
\description{
Build a Phase II Shewhart control chart based on the unconditional method and/or the exceedance probability criterion}
\examples{
# Build a Phase I ANOVA-based Shewhart control chart
X <- PH2XBAR.data()
PH2XBAR(
    X
    ,PH1.info = list(
    				X = NULL
    				, mu = 245
    				, sigma = 2
    				, k = 30
    				, n = 5
    				, model = "ANOVA-based"
    			)
    ,c.ii.info = list(
    				c.ii = NULL
    				, method = 'both'
    				, ARL0 = 370
    				, p = 0.05
    				, eps = 0.1
    				, interval.c.ii.UC = c(1, 3.2)
    				, interval.c.ii.EPC = c(1, 10)
    				, UC.tol = .Machine$double.eps^0.25
            , EPC.tol = .Machine$double.eps^0.25
    			)
)
}
