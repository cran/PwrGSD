\name{plot.cpd.PwrGSD}
\alias{plot.cpd.PwrGSD}
\title{Plot Method for cpd.PwrGSD objects}
\usage{
\method{plot}{cpd.PwrGSD}(x, formula, subset, na.action,\dots)
}

\description{Creates a trellis plot of type II error probability and
  power at each interim analysis, stacked, versus an effect size
  variable, conditioned upon levels of up to two factors.}

\arguments{
  \item{x}{an object of class \code{cpd.PwrGSD}}
  \item{formula}{a one sided formula of the form \code{~ effect | f1}
    or \code{~ effect | f1 * f2} where \code{effect}, \code{f1}, and
    \code{f2} are variables in the indexing dataframe \code{descr}, or
    the special variable \code{stat} which may be used when there are
    multiple test statistics per component of \code{Elements}.  See
    the example in the documentation for \code{cpd.PwrGSD}}.
  \item{subset}{the plot can be applied to a subset of rows of
    \code{descr} via a logical expression on its variables in
    combination with the special variable, \code{stat} when applicable.}
  \item{na.action}{a \code{na.action} method for handling \code{NA}
    values}
  \item{\dots}{other parameters to pass to the R function \code{coplot}
    usually not neccesary}
}

\value{
  Returns the object, \code{x}, invisibly
}

\note{This processes the \code{cpd.PwrGSD} object into a dataframe,
  stacked on interim looks and then passes the results to the R
  function \code{coplot}}
\references{
  Chambers, J. M. (1992)
  \emph{Data for models.}
  Chapter 3 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth \& Brooks/Cole.

  Cleveland, W. S. (1993) \emph{Visualizing Data.} New Jersey: Summit Press.
}
\author{Abovementioned \code{cpd.PwrGSD} processing done by Grant
  Izmirlian <izmirlian@nih.gov>}

\seealso{\code{\link{cpd.PwrGSD}} \code{\link{Power}} and \code{\link{Elements}}}
\examples{
  ## See the example in the 'cpd.PwrGSD' documentation
}
\keyword{survival}
\keyword{hplot}
