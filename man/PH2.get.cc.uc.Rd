\name{PH2.get.cc.uc}
\alias{PH2.get.cc.uc}
\title{PH2.get.cc.uc}
 \usage{
PH2.get.cc.uc(
	ARL0,
	k,
	nu = k - 1,
	c4.option = TRUE,
	interval = c(1, 10),
	u = runif(100000),
	v = runif(100000),
	tol = .Machine$double.eps^0.25
)
}
\arguments{
  \item{ARL0}{ a nominal in-control average run length }
  \item{k}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{interval}{ the searching range of Phase II charting constant }
  \item{u}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{v}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{tol}{ a tolerance level for the convergence }

}
\description{
Obtain a Phase II charting constant by the unconditional method.
}
\examples{
PH2.get.cc.uc(ARL0 = 370,	k = 30)
}
