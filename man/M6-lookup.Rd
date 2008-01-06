\name{lookup}
\alias{lookup}
\title{Lookup values for a piecewise constant function}
\description{
  Given the values and lefthand endpoints for intervals of constancy,
  lookup values of the function at arbitrary values of the independent
  variable.
}
\usage{
  lookup(xgrid, ygrid, x, y0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xgrid}{Lefthand endpoints of intervals of constancy}
  \item{ygrid}{Values on these intervals, of same length as \code{xgrid}}
  \item{x}{Input vector of arbitrary independent variables.}
  \item{y0}{Value to be returned for values of \code{x} that are smaller
    than \code{min(xgrid)}.}
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\examples{
  ## none as yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
