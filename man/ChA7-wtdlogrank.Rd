\name{wtdlogrank}
\alias{wtdlogrank}
\title{Weighted log-rank test}
\description{
  Computes a two sample weighted log-rank statistic with events weighted
  according to one of the available weighting function choices}
\usage{
  wtdlogrank(formula = formula(data), data = parent.frame(), WtFun = c("FH", "SFH", "Ramp"),
  param = c(0, 0), sided = c(2, 1), subset, na.action, w = FALSE)
}
\arguments{
  \item{formula}{a formula of the form \code{Surv(Time, Event) ~ arm} where \code{arm}
    is a dichotomous variable with values 0 and 1.}
  \item{data}{a dataframe}
  \item{WtFun}{a selection from the available list: \dQuote{FH} (Fleming-Harrington),
    \dQuote{SFH} (stopped Fleming-Harrington) or \dQuote{Ramp}.  See \code{param} in
    the following line.}
  \item{param}{Weight function parameters. Length and interpretation depends upon
    the selected value of \code{WtFun}:
    If \code{WtFun==FH} then \code{param} is a length 2 vector specifying the power
    of the pooled (across arms) kaplan meier estimate and its complement.
    If \code{WtFun==SFH} then \code{param} is a length 3 vector with first two components
    as in the \dQuote{FH} case, and third component the time (in the same units as the time
    to event) at which the \dQuote{FH} weight function is capped off at its current value.
    If \code{WtFun==SFH} then \code{param} is of length 1 specifying the time (same units
    as time to event) at which events begin to get equal weight. The \dQuote{Ramp} weight
    function is a linearly increasing deterministic weight function which is capped off at 1
    at the user specified time.
     }
  \item{sided}{One or Two sided test?  Set to 1 or 2}
  \item{subset}{Analysis can be applied to a subset of the dataframe based upon a logical expression
    in its variables}
  \item{na.action}{Method for handling \code{NA} values in the covariate, \code{arm}}
  \item{w}{currently no effect}
}
\value{
  An object of class \code{survtest} containing components
  \item{pn}{sample size}
  \item{wttyp}{internal representation of the \code{WtFun} argument}
  \item{par}{internal representation of the \code{param} argument}
  \item{time}{unique times of events accross all arms}
  \item{nrisk}{number at risk accross all arms at each event time}
  \item{nrisk1}{Number at risk in the experimental arm at each event time}
  \item{nevent}{Number of events accross all arms at each event time}
  \item{nevent1}{Number of events in the experimental arm at each event time}
  \item{wt}{Values of the weight function at each event time}
  \item{pntimes}{Number of event times}
  \item{stat}{The un-normalized weighted log-rank statistic, i.e. the summed
    weighted observed minus expected differences at each event time}
  \item{var}{Variance estimate for the above}
  \item{UQt}{Cumulative sum of increments in the sum resulting in \code{stat} above}
  \item{varQt}{Cumulative sum of increments in the sum resulting in \code{var} above}
  \item{var1t}{Cumulative sum of increments in the sum resulting in the variance of an
    unweighted version of the statistic}
  \item{pu0}{person units of follow-up time in the control arm}
  \item{pu1}{person units of follow-up time in the intervention arm}
  \item{n0}{events in the control arm}
  \item{n1}{events in the intervention arm}
  \item{n}{sample size, same as \code{pn}}
  \item{call}{the call that created the object}
}
\references{
  Harrington, D. P. and Fleming, T. R. (1982).
  A class of rank test procedures for censored survival data.
  \emph{Biometrika}
  \bold{69}, 553-566.}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\seealso{\code{\link{IntSurvDiff}}}
\examples{
  library(PwrGSD)
  data(lung)
  fit.wlr <- wtdlogrank(Surv(time, I(status==2))~I(sex==2), data=lung, WtFun="SFH", param=c(0,1,300))
}
\keyword{survival}
