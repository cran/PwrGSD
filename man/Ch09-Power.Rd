\name{Power}
\alias{Power}
\title{Extract the Power results}
\description{
  The function `Power' is used to summarize the `cpd.PwrGSD' object into  a
  dataframe containing power and type II error, summed over analysis
  times. The data frame is stacked by rows of `descr' and by `stat'
  (if there are multiple statistics being profiled per each component of
  `Elements'), for generating tables or performing other computations.
  
}
\usage{
Power(object, subset, nlook.ind = NULL)
}
\arguments{
  \item{object}{an object of class \code{cpd.PwrGSD}}
  \item{subset}{you may extract a subset via a logical expression in the
    variables of the index dataframe, \code{descr}}
  \item{nlook.ind}{(optional) a vector containing a subset of the
    indices of analysis times over which the sum is formed. Use this
    for example if you want to know the probability of stopping by the
    kth analysis under an unfavorable alternative.  Set nlook.ind to
    1:k}
}
\value{
  a dataframe, stacked by rows of `descr' and then by choices of `stat'
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\seealso{\code{\link{cpd.PwrGSD}} and \code{\link{PwrGSD}}}
\examples{
  ## See the `cpd.PwrGSD' example
}
\keyword{data}
