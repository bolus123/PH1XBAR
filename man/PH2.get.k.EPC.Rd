\name{PH2.get.k.EPC}
\alias{PH2.get.k.EPC}
\title{PH2.get.k.EPC}
 \usage{
PH2.get.k.EPC(
	p, 
	c.ii, 
	eps = 0.1, 
	ARL0 = 370, 
	c4.option = TRUE, 
	interval = c(1000, 3000), 
	model = 'ANOVA-based', 
	n = 10, 
	u = runif(100000), 
	v = runif(100000),
	tol = .Machine$double.eps^0.25
)
}
\arguments{
  \item{p}{ exceedance probability }
  \item{c.ii}{ a nominal Phase II charting constant }
  \item{eps}{ a nominal allowable changes of the in-control average run length }
  \item{ARL0}{ a nominal in-control average run length }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{interval}{ the searching range of Phase I batch size }
  \item{model}{ 'ANOVA-based' - degrees of freedom equal to k - 1 and estimating the within-batch standard deviation. 'basic' - degree of freedo equal to k * (n - 1) and estimating the pooled standard deviation }
  \item{n}{ subgroup size }
  \item{u}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{v}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{tol}{ a tolerance level for the convergence }

}
\description{
Obtain a Phase I batch size and a potential Phase I sample size by exceedance probability criterion.
}
\examples{
PH2.get.k.EPC(
	p = 0.1, 
	c.ii = 3,
	interval = c(1000, 500000)
)
}
