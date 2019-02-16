\name{PH2.get.ARLb.EPC}
\alias{PH2.get.ARLb.EPC}
\title{PH2.get.ARLb.EPC}
 \usage{
PH2.get.ARLb.EPC(
	p, 
	k, 
	nu = k - 1, 
	c.ii, 
	ARL0 = NULL, 
	c4.option = TRUE, 
	interval = c(1, 1000), 
	u = runif(100000), 
	v = runif(100000),
	tol = .Machine$double.eps^0.25
)
}
\arguments{
  \item{p}{ exceedance probability }
  \item{k}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{c.ii}{ a nominal Phase II charting constant }
  \item{ARL0}{ a nominal in-control average run length }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{interval}{ the searching range of lower prediction bound of average run length }
  \item{u}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{v}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{tol}{ a tolerance level for the convergence }
}
\description{
Obtain the lower prediction bound of average run length
}
\examples{
PH2.get.ARLb.EPC(p = 0.05, k = 30, c.ii = 3)
}

