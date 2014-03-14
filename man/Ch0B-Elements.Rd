\name{Elements}
\alias{Elements}
\title{Create a subset of a "cpd.PwrGSD" object}
\description{
  Create a subset of a \code{cpd.PwrGSD} object
}
\usage{
  Elements(object, subset, na.action = na.pass)
}
\arguments{
  \item{object}{an object of class \code{cpd.PwrGSD}}
  \item{subset}{you may extract a subset via a logical expression in the
    variables of the index dataframe, \code{descr}}
  \item{na.action}{a method for handling \code{NA} values in the
    variables in the subset expression.}
}
\value{
  an object of class \code{cpd.PwrGSD}. See help on that topic for
  details.
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\seealso{\code{\link{cpd.PwrGSD}} and \code{\link{PwrGSD}}}
\examples{
  ## See the `cpd.PwrGSD' example
}
\keyword{data}

