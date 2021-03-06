\name{EX1gXK}
\alias{EX1gXK}
\title{
  A function for computing the bias adjusted point estimate for a
  statistic observed to cross the efficacy boundary. 
}
\description{
  A function for computing the bias adjusted point estimate for a
  statistic, on the Brownian scale, observed to cross the efficacy boundary. 
}
\usage{
EX1gXK(xk, b.eff, frac)
}
\arguments{
  \item{xk}{
    The observed value of the statistic, on the \dQuote{Brownian} scale.
  }
  \item{b.eff}{
    Efficacy boundary points at current and prior analyses
  }
  \item{frac}{
    Information fraction at current and prior analyses
  }
}
\value{
  Returns the expected value of X_1 given X_K, which is the bias adjusted point estimate
}
\references{
Emerson, S. S. (1993). Computation of the uniform minimum variance unibiased estimator of
a normal mean following a group sequential trialdiscrete sequential boundaries for
clinical trials. Computers and Biomedical Research 26 68--73.

Izmirlian, G. (2014). Estimation of the relative risk following group sequential procedure
based upon the weighted log-rank statistic. Statistics and its Interface 00 00--00
}
\author{
Grant Izmirlian <izmirlig@mail.nih.gov>
}
\note{
This works for the unweighted, proportional hazards case, but also works
in the case of the weighted log-rank statistic when we assume the chosen
weights are proportional to the true shape.
}
\seealso{
\code{\link{gsd.dens}}
}
\examples{
# if Z.K = U_K/V_K^0.5 is the log-rank statistic on the standard normal
# scale, then we obtain an estimate of the logged relative risk as follows
# Suppose we've stopped at analysis number K=4, and Z.K = 2.5
# suppose the end of trial variance of the log-rank statistic
# (specified in design and used to compute 'frac') is V.end = 100

K <- 4
Z.K <- 2.5
V.end <- 100

# Information fraction
frac <- c(0.15, 0.37, 0.64, 0.76)

# Efficacy Boundary 
gsb <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(spending=ObrienFleming, alpha=0.05))

# Efficacy boundary points
be <- gsb$table[,"b.e"] 

# Brownian scale
X.K <- Z.K*frac[K]

# expected value of X_1 given X_K
ex1gxk <- EX1gXK(X.K, be, frac)

# Crude estimate of logged relative risk
X.K/(frac[K]*V.end^0.5)

# Bias adjusted estimate of logged relative risk
ex1gxk/(frac[1]*V.end^0.5)
}
\keyword{Sequential Design}
\keyword{Bias Adjusted point estimate}
