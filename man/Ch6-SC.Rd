\name{SC}
\alias{SC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Stochastic Curtailment method of Boundary Construction}
\description{
  The function \code{SC} is used in calls to the functions
  \code{GrpSeqBnds} and \code{PwrGSD} as a possible setting for the
  arguments \code{EfficacyBoundary} and \code{FutilityBoundary}, in
  specification of the method whereby efficacy and or futility
  boundaries are to be constructed. The Stochastic Curtailment method is
  one of four currently availiable choices, the others being
  \code{LanDemets}, \code{Haybittle} and user specified.
}
\usage{
SC(be.end, prob, drift.end = NULL, from = NULL, to = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{be.end}{The value of the efficacy criterion in the scale of a
    standardized normal. This should be set to something further from
    the null than the single test $Z_alpha$. For example if the total type
    I error probability is 0.05 in a two sided test of the null than
    set be.end to 2.10 or larger (instead of 1.96).}
  \item{prob}{The criterion, a probability to be exceeded in
    order to stop. 0.90 or above is a good choice. See detail below.}
  \item{drift.end}{Required only if you are using \code{SC} to set the
    \code{FutilityBoundary}. In this case, set \code{drift.end} to the
    value of the drift function anticipated at the end of the
    trial. See detail below.}
  \item{from}{WARNING EXPERIMENTAL: you can actually construct boundaries via a
    hybrid of the 3 boundary construction methods, \code{LanDemets}, \code{SC},
    and 'user specified'. When using a hybrid boundry, set the argument
    \code{EfficacyBoundary} or \code{FutilityBoundary} respectively,
    to a list with components \code{LanDemets}, \code{SC}, or user
    specified numbers. In the former two cases, \code{from} and
    \code{to} are used in \code{LanDemets} and also in \code{SC}
    to stipulate how many interim analyses they are in effect.
    See the help for \code{GrpSeqBnds} and \code{PwrGSD}}
  \item{to}{See above.}
}
\details{
  When the stochastic curtailment procedure is used to construct the
  efficacy boundary, i.e. \code{EfficacyBoundary}=\code{SC(...)},
  the efficacy criterion is reached when the conditional probability,
  under the null hypothesis, that the last analysis results in
  statistical significance, given the present value of the statistic,
  exceeds 'prob'. In of itself, this doesn't produce a boundary on the
  scale of a standard normal, but it is easily converted to one as is
  done here. When this is used to construct a futility boundary, i.e.
  \code{FutilityBoundary}=\code{SC(...)}, the futility criterion is
  reached when the conditional probability, under the design alternative
  hypothesis, that the last analysis does not result in statistical
  significance, given the present value of the statistic, exceeds 'prob'.
  The design alternative corresponds to a drift function, which is the
  expected value of the statistic normalized to have variance equal to
  the information fraction at each interim analysis. For the unweighted
  log-rank statistic, the drift function is $(V_T)^(1/2)$ B f, where
  B is the logged relative risk, $V_T$ is the variance at the end of the
  trial and f is the information fraction.  If the two trial arms are
  balanced and the number at risk is roughly constant throughout the
  trial then $V_T = pi (1-pi) N_T$, where $pi$ is the constant proportion
  at risk in one of the trial arms and $N_T$ is the anticipated number of
  events.
}  

\value{
  An object of class \code{boundary.construction.method} which is really a list
  with the following components. The print method displays the original
  call.
  \item{type}{Gives the boundary construction method type, which is the character
    string "SC"}
  \item{be.end}{The numeric value passed to the argument 'be.end', which is the 
    value of the efficacy criterion in the scale of a standardized normal.}
  \item{prob}{The numeric value passed to the argument 'prob', which is the probability to be
    exceeded in order to stop.}
  \item{drift.end}{The numeric value passed to the argument 'drift.end', which is the
    value of the drift function at the end of the trial. See details.}
  \item{from}{The numeric value passed to the argument 'from'. See above.}
  \item{to}{The numeric value passed to the argument 'to'. See above.}
  \item{call}{returns the call}
}

\references{see references under \code{\link{PwrGSD}}}
\author{Grant Izmirlian}
\note{The print method returns the call by default} 
\seealso{\code{\link{LanDemets}}, \code{\link{GrpSeqBnds}}, \code{\link{PwrGSD}}}
\examples{
## example 1: what is the result of calling a Boundary Construction Method function
    ## A call to 'SC' just returns the call
    SC(be.end=2.10, prob=0.90)
    
    ## It does arguement checking...this results in an error
    \dontrun{
      SC(be.end=2.10)
    }
    
    ## but really its value is a list with the a component containing
    ## the boundary method type, "LanDemts", and components for each
    ## of the arguments.
    names(SC(be.end=2.10, prob=0.90))

    SC(be.end=2.10, prob=0.90, drift.end=2.34)$type
    SC(be.end=2.10, prob=0.90, drift.end=2.34)$be.end
    SC(be.end=2.10, prob=0.90, drift.end=2.34)$prob
    SC(be.end=2.10, prob=0.90, drift.end=2.34)$drift.end

## example 2: ...But the intended purpose of the spending functions its
    ## intended purpose is in constructing calls to 'GrpSeqBnds' and to 'PwrGSD':
     

    frac <- c(0.07614902,0.1135391,0.168252,0.2336901,0.3186155,
              0.4164776,0.5352199,0.670739,0.8246061,1)
    drift <- c(0.3836636,0.5117394,0.6918584,0.8657705,1.091984,
	      1.311094,1.538582,1.818346,2.081775,2.345386)

    test <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(alpha=0.05, spending=ObrienFleming),
                       FutilityBoundary=SC(be.end=2.10, prob=0.90, drift.end=drift[10]),
                       drift=drift)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{design}
