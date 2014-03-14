\name{mystack}
\alias{mystack}
\title{Stack a dataset}
\description{
  Given a dataframe containing one or more variables named with a common
  prefix, this function creates a stacked dataset with one set of
  observed values of the variables (in order of occurence) per line.
}
\usage{
  mystack(object, fu.vars, create.idvar = FALSE)
}
\arguments{
  \item{object}{a dataframe containing one or more variables named with a common
    prefix}
  \item{fu.vars}{a list of the unique prefixes}
  \item{create.idvar}{Do you want to add an ID variable with a common
    value given to all records resulting from a given input record?
    Default is \code{FALSE}}
}
\value{
  A stacked dataframe
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\examples{
  ## none as yet
}
\keyword{data}

