\name{GrpSeqBnds}
\alias{GrpSeqBnds}
\title{Computes efficacy and futility boundaries}
\description{
  This computes efficacy and futility boundaries for interim analysis
  and sequential designs.  Two sided symmetric efficacy boundaries can
  be computed by specifying half of the intended total type I error
  probability in the argument, \code{Alpha.Efficacy}.  Otherwise,
  especially in the case of efficacy and futility bounds only one sided
  boundaries are currently computed.  The computation allows for two
  different time scales--one must be the variance ratio, and the second
  can be a user chosen increasing scale beginning with 0 that takes the
  value 1 at the conclusion of the trial.
}
\usage{
GrpSeqBnds(EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming),
           FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming),
           NonBindingFutility = TRUE, frac, frac.ii = NULL, drift = NULL)
}

\arguments{
  \item{EfficacyBoundary}{This specifies the method used to construct the efficacy boundary. The available choices are:
    
    \samp{(i)  }\code{Lan-Demets(alpha=<total type I error>, spending=<spending function>)}. The Lan-Demets method is based upon
    a error probability spending approach. The spending function can be set to \code{ObrienFleming}, \code{Pocock}, or
    \code{Power(rho)}, where \code{rho} is the the power argument for the power spending function: rho=3 is roughly equivalent
    to the O'Brien-Fleming spending function and smaller powers result in a less conservative spending function.
    
    \samp{(ii)  }\code{Haybittle(alpha=<total type I error>, b.Haybittle=<user specified boundary point>)}. The Haybittle
    approach is the simplest, which sets the boundary points equal to \code{b.Haybittle}, a user specified value (try 3) for
    all analyses except the last, which is calculated so as to result in the total type I error, set with the argument
    \code{alpha}.
    
    \samp{(iii)  }\code{SC(be.end=<efficacy boundary point at trial end>, prob=<threshold for conditional type I error for
    efficacy stopping>)}. The stochastic curtailment method is based upon the conditional probability of type I error given the
    current value of the statistic. Under this method, a sequence of boundary points on the standard normal scale (as are
    boundary points under all other methods) is calculated so that the total probability of type I error is  maintained. This is
    done by considering the joint probabilities of continuing to the current analysis and then exceeding the threshold at the
    current analysis. A good value for the threshold value for the conditional type I error, \code{prob} is 0.90 or greater.
    
    \samp{(iv)  }User supplied boundary points in the form \code{c(b1, b2, b3, ..., b_m)}, where \code{m} is the number of looks.}
  \item{FutilityBoundary}{This specifies the method used to construct the futility boundary. The available choices are:
    
    \samp{(i)  }\code{Lan-Demets(alpha=<total type II error>, spending=<spending function>)}. The Lan-Demets method is based upon
    a error probability spending approach. The spending function can be set to \code{ObrienFleming}, \code{Pocock}, or
    \code{Power(rho)}, where \code{rho} is the the power argument for the power spending function: rho=3 is roughly equivalent
    to the O'Brien-Fleming spending function and smaller powers result in a less conservative spending function.
    
    \samp{(ii)  }\code{Haybittle(alpha=<total type I error>, b.Haybittle=<user specified boundary point>)}. The Haybittle
    approach is the simplest, which sets the boundary points equal to \code{b.Haybittle}, a user specified value (try 3) for
    all analyses except the last, which is calculated so as to result in the total type II error, set with the argument
    \code{alpha}.
    
    \samp{(iii)  }\code{SC(be.end=<efficacy boundary point at trial end>, prob=<threshold for conditional type II error for
    futility stopping>, drift.end=<projected drift at end of trial>)}. The stochastic curtailment method is based upon the
    conditional probability of type II error given the current value of the statistic. Under this method, a sequence of boundary
    points on the standard normal scale (as are boundary points under all other methods) is calculated so that the total
    probability of type II error, is maintained. This is done by considering the joint probabilities of continuing to the
    current analysis and then exceeding the threshold at the current analysis. A good value for the threshold value for the
    conditional type I error, \code{prob} is 0.90 or greater.
    
    \samp{(iv)  }User supplied boundary points in the form \code{c(b1, b2, b3, ..., b_m)}, where \code{m} is the number of looks.}
  \item{NonBindingFutility}{When using a futility boundary and this is set to 'TRUE', the efficacy boundary will be constructed   
    in the absence of the futility boundary, and then the futility boundary will be constructed given the resulting efficacy      
    boundary. This results in a more conservative efficacy boundary with true type I error less than the nominal level. This is   
    recommended due to the fact that futility crossings are viewed by DSMB's with much less gravity than an efficacy crossing     
    and as such, the consensus is that efficacy bounds should not be discounted towards the null hypothesis because of paths      
    which cross a futility boundary. Default value is 'TRUE'.}
  \item{frac}{The variance ratio. If the end of trial variance is unknown then normalize all previous variances by the current
    variance. In this case you must specify a second scale that is monotone increasing from 0 to 1 at the end of the
    trial. Required.}
  \item{frac.ii}{The second information scale that is used for type I and type II error probability spending. Optional (see
    above)}
  \item{drift}{The drift function of the underlying brownian motion, which is the expected value under the design alternative
    of the un-normalized weighted log-rank statistic, then normalized to have variance one when the variance ratio equals 1. See
    the examples below.}
}

\value{
  An object of class \code{boundaries} with components:
  "table"   "frac"    "frac.ii" "drift"   "call"
  \item{call}{The call that produced the returned results.}
  \item{frac}{The vector of variance ratios.}
  \item{frac.ii}{The vector of information ratios for type I and type II error probability spending, which differs from the
    above if the user sets the argument \code{frac.ii} to a second scale as mentioned above.}
  \item{drift}{The drift vector that is required as an argument when futility boundaries are calculated.}
  \item{table}{A matrix with components
    
    \samp{frac         }The information ratio for type I and type II error probability spending.
    
    \samp{b.f          }The calculated futility boundary (if requested).

    \samp{alpha.f      }The type II error probability spent at that analysis (if doing futility bounds).

    \samp{cum-alpha.f  }Cumulative sum of \code{alpha.f} (if doing futility bounds).
    
    \samp{b.e          }The calculated efficacy boundary.
    
    \samp{alpha.e      }The type I error probability spent at that analysis.

    \samp{cum-alpha.e  }Cumulative sum of \code{alpha.e}.
  }
}
\references{
  Gu, M.-G. and Lai, T.-L.  \dQuote{Determination of power and sample size in the design of clinical trials with failure-time
  endpoints and interim analyses.}
  Controlled Clinical Trials 20 (5): 423-438. 1999

  Izmirlian, G.  \dQuote{The PwrGSD package.}  NCI Div. of Cancer Prevention Technical
  Report. 2004

  Jennison, C. and Turnbull, B.W. (1999)
  Group Sequential Methods: Applications to Clinical Trials
  Chapman & Hall/Crc, Boca Raton FL

  Proschan, M.A., Lan, K.K.G., Wittes, J.T.
  (2006), corr 2nd printing (2008)
  Statistical Monitoring of Clinical Trials A Unified Approach
  Springer Verlag, New York
}
\author{Grant Izmirlian <izmirlian@nih.gov>}

\seealso{\code{\link{PwrGSD}}}
\examples{

    ## NOTE: In an unweighted analysis, the variance ratios and event ratios
    ## are the same, whereas in a weighted analysis, they are quite different.
    ##
    ## For example, in a trial with 7 or so years of accrual and maximum follow-up of 20 years
    ## using the stopped Fleming-Harrington weights, `WtFun' = "SFH", with paramaters
    ## `ppar' = c(0, 1, 10) we might get the following vector of variance ratios:

    frac    <- c(0.006995655, 0.01444565, 0.02682463, 0.04641363, 0.0585665,
                 0.07614902, 0.1135391, 0.168252, 0.2336901, 0.3186155, 0.4164776,
                 0.5352199, 0.670739, 0.8246061, 1)


    ## and the following vector of event ratios:

    frac.ii <- c(0.1494354, 0.1972965, 0.2625075, 0.3274323, 0.3519184, 0.40231,
                 0.4673037, 0.5579035, 0.6080742, 0.6982293, 0.7671917, 0.8195019,
                 0.9045182, 0.9515884, 1)

    ## and the following drift under a given alternative hypothesis
               
    drift <-   c(0.06214444, 0.1061856, 0.1731267, 0.2641265, 0.3105231, 0.3836636,
                 0.5117394, 0.6918584, 0.8657705, 1.091984, 1.311094, 1.538582,
                 1.818346, 2.081775, 2.345386)

    ## JUST ONE SIDED EFFICACY BOUNDARY
    ## In this call, we calculate a one sided efficacy boundary at each of 15 analyses
    ## which will occur at the given (known) variance ratios, and we use the variance
    ## ratio for type I error probability spending, with a total type I error probabilty
    ## of 0.05, using the Lan-Demets method with Obrien-Fleming spending (the default).

    gsb.all.just.eff <- GrpSeqBnds(frac=frac, 
                                   EfficacyBoundary=LanDemets(alpha=0.05, spending=ObrienFleming))

    ## ONE SIDED EFFICACY AND FUTILTY BOUNDARIES
    ## In this call, we calculate a one sided efficacy boundary at each of 15 analyses
    ## which will occur at the given (known) variance ratios, and we use the variance
    ## ratio for type I and type II error probability spending, with a total type I error
    ## probabilty of 0.05 and a total type II error probability of 0.10, using the Lan-Demets
    ## method with Obrien-Fleming spending (the default) for both efficacy and futilty.

    gsb.all.eff.fut <- GrpSeqBnds(frac=frac,  
                                  EfficacyBoundary=LanDemets(alpha=0.05, spending=ObrienFleming),
                                  FutilityBoundary=LanDemets(alpha=0.10, spending=ObrienFleming),
                                  drift=drift)

    ## Now suppose that we are performing the 7th interim analysis. We don't know what the variance
    ## will be at the end of the trial, so we normalize variances of the current and previous
    ## statistics by the variance of the current statistic.  This is equivalent to the following
    ## length 7 vector of variance ratios:

    frac7 <- frac[1:7]/frac[7]

    ## To proceed under the "unknown variance at end of trial" case, we must use a second
    ## scale for spending type I and II error probabilty. Unlike the above scale
    ## which is renormalized at each analysis to have value 1 at the current analysis, the
    ## alpha spending scale must be monotone increasing and attain the value 1 only at the
    ## end of the trial. A natural choice is the event ratio, which is known in advance if
    ## the trial is run until a required number of events is obtained, a so called
    ## maximum information trial:

    frac7.ii <- frac.ii[1:7]

    ## the first seven values of the drift function

    drift7 <- drift[1:7]/frac[7]^0.5

    gsb.1st7.eff.fut <- GrpSeqBnds(frac=frac7, frac.ii=frac7.ii,  
                                   EfficacyBoundary=LanDemets(alpha=0.05, spending=ObrienFleming),
                                   FutilityBoundary=LanDemets(alpha=0.10, spending=ObrienFleming),
				   drift=drift7)

    ## Of course there are other options not covered in these examples but this should get you
    ## started
}

\keyword{design}
\keyword{htest}
