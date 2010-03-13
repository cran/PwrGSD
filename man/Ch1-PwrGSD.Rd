\name{PwrGSD}
\alias{PwrGSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Power in a Group Sequential Design}
\description{

  Derives power in a two arm clinical trial under a group sequential design.  Allows for arbitrary number of interim analyses,
  arbitrary specification of arm-0/arm-1 time to event distributions (via survival or hazard), arm-0/arm-1 censoring
  distribution, provisions for two types of continuous time non-compliance according to arm-0/arm-1 rate followed by switch to
  new hazard rate.  Allows for analyses using (I) weighted log-rank statistic, with weighting function (a) a member of the
  Flemming-Harrington G-Rho class, or (b) a stopped version thereof, or (c) the ramp-plateau deterministic weights, or (II) the
  integrated survival distance (currently under method=="S" without futility only).  Stopping boundaries are computed via the
  Lan-Demets method, Haybittle method, or converted from the stochastic curtailment procedure.  The Lan-Demets boundaries can be
  constructed usign either O'Brien-Flemming, Pocock or Wang-Tsiatis power alpha-spending.  The C kernel is readily extensible,
  and further options will become available in the near future.}

\usage{
PwrGSD(EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming),
    FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming),
    NonBindingFutility = TRUE, sided = c("2>", "2<", "1>", "1<"),
    method = c("S", "A"), accru, accrat, tlook,
    tcut0 = NULL, h0 = NULL, s0 = NULL, tcut1 = NULL,
    rhaz = NULL, h1 = NULL, s1 = NULL, tcutc0 = NULL, hc0 = NULL,
    sc0 = NULL, tcutc1 = NULL, hc1 = NULL, sc1 = NULL, tcutd0A = NULL,
    hd0A = NULL, sd0A = NULL, tcutd0B = NULL, hd0B = NULL, sd0B = NULL,
    tcutd1A = NULL, hd1A = NULL, sd1A = NULL, tcutd1B = NULL,
    hd1B = NULL, sd1B = NULL, tcutx0A = NULL, hx0A = NULL, sx0A = NULL,
    tcutx0B = NULL, hx0B = NULL, sx0B = NULL, tcutx1A = NULL,
    hx1A = NULL, sx1A = NULL, tcutx1B = NULL, hx1B = NULL, sx1B = NULL,
    noncompliance = c("none", "crossover", "mixed", "user"),
    gradual = FALSE, WtFun = c("FH", "SFH", "Ramp"), ppar = cbind(c(0, 0)), 
    Spend.Info = c("Variance", "Events", "Hybrid(k)", "Calendar"), RR.Futility = NULL, 
    qProp.one.or.Q = c("one", "Q"), Nsim = NULL, detail = FALSE, StatType = c("WLR",
        "ISD"))
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
    boundary points under all other methods) is calculated so that the total probability of type I error is maintained. This is
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
  \item{sided}{Set to \dQuote{2>} (quoted) for two sided tests of the null hypothesis when a positive drift corresponds to
    efficacy. Set to \dQuote{2<} (quoted) for two sided tests of the null hypothesis when a negative drift corresponds to
    efficacy. Set to \dQuote{1>} or \dQuote{1<} for one sided tests of H0 when efficacy corresponds to a positive or negative drift,
    respectively. If \code{method}==\dQuote{S} then this must be of the same length as \code{StatType} because the
    interpretation of \code{sided} is different depending upon whether \code{StatType}==\dQuote{WLR} 
    (negative is benefit) or \code{StatType}==\dQuote{ISD} (positive is benefit)}
  \item{method}{Determines how to calculate the power.  Set to \dQuote{A} (Asymptotic method, the default) or \dQuote{S}
    (Simulation method)} 
  \item{accru}{The upper endpoint of the accrual period beginning with time 0.}
  \item{accrat}{The rate of accrual per unit of time.}
  \item{tlook}{The times of planned interim analyses.}
  \item{tcut0}{Left hand endpoints for intervals upon which the arm-0 specific mortality is constant. The last given component
    is the left hand endpoint of the interval having right hand endpoint infinity.}
  \item{h0}{A vector of the same length as \code{tcut0} which specifies the piecewise constant arm-0 mortality rate.} 
  \item{s0}{Alternatively, the arm-0 mortality distribution can be supplied via this argument, in terms of of the corresponding
    survival function values at the times given in the vector \code{tcut0}. If \code{s0} is supplied, then \code{h0}is derived
    internally, assuming the piecewise exponential distrubiton. If you specify \code{s0}, the first element must be 1, and
    correspondingly, the first component of \code{tcut0} will be the lower support point of the distribution. You must supply
    either \code{h0} or \code{s0} but not both.}  
  \item{tcut1}{Left hand endpoints for intervals upon which the arm-1 specific mortality is constant.  The last given component
    is the left hand endpoint of the interval having right hand endpoint infinity.}  
  \item{rhaz}{A vector of piecewise constant arm-1 versus arm-0 mortality rate ratios.  If \code{tcut1} and \code{tcut0} are not
    identical, then \code{tcut1}, \code{h0}, and \code{rhaz} are internally rederived at the union of the sequences \code{tcut0}
    and \code{tcut1}. In all cases the arm-1 mortality rate is then derived at the time cutpoints \code{tcut1} as \code{rhaz}
    times\code{h0}.}   
  \item{h1}{Alternatively, the arm-1 mortality distribution can be supplied via this argument by specifying the piecewise
    constant arm-1 mortality rate. See the comments above.}
  \item{s1}{Alternatively, the arm-1 mortality distribution can be supplied via this argument, in terms of of the corresponding
    survival function values at the times given in the vector \code{tcut1}. Comments regarding \code{s0} above apply here as
    well. You must supply exactly one of the following: \code{h1}, \code{rhaz}, or \code{s1}.}
  \item{tcutc0}{Left hand endpoints for intervals upon which the arm-0 specific censoring distribution hazard function is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity.}
  \item{hc0}{A vector of the same length as \code{tcutc0} which specifies the arm-0 censoring distribution in terms of a
    piecewise constant hazard function.} 
  \item{sc0}{Alternatively, the arm-0 censoring distribution can be supplied via this argument, in terms of of the corresponding
    survival function values at the times given in the vector \code{tcutc0}. See comments above. You must supply either
    \code{hc0} or \code{sc0} but not both.}  
  \item{tcutc1}{Left hand endpoints for intervals upon which the arm-1 specific censoring distribution hazard function is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity.}
  \item{hc1}{A vector of the same length as \code{tcutc1} which specifies the arm-1 censoring distribution in terms of a
    piecewise constant hazard function.}  
  \item{sc1}{Alternatively, the arm-1 censoring distribution can be supplied via this argument, in terms of of the
    corresponding survival function values at the times given in the vector \code{tcutc1}. See comments above. You must supply
    either \code{hc1} or \code{sc1} but not both.}
  \item{noncompliance}{(i) Seting \code{noncompliance} to \dQuote{none} for no non-compliance will automatically set the
    non-compliance arguments, below, to appropriate values for no compliance. This requires no additional user specification of 
    non-compliance parameters.  (ii) Setting \code{noncompliance} to \dQuote{crossover} will automatically set crossover values
    in the arm 0/1 specific \emph{post-cause-B-delay-mortality} for cross-over, i.e. simple interchange of the arm 0 and arm 1
    mortalities. The user is required to specify all parameters corresponding to the arm 0/1 specific \emph{cause-B-delay}
    distributions. The \emph{cause-A-delay} and \emph{post-cause-A-delay-mortality} are automatically set so as not to influence
    the calculations. Setting \code{noncompliance} to \dQuote{mixed} will set the arm 0/1 specific
    \emph{post-cause-B-delay-mortality} distributions for crossover as defined above. The user specifies the arm 0/1 specific 
    \emph{cause-B-delay} distribution as above, and in addition, all parameters related to the arm 0/1 specific
    \emph{cause-A-delay} distributions and corresponding arm 0/1 specific \emph{post-cause-A-delay-mortality}
    distributions. (iii) Setting \code{noncompliance} to \dQuote{user} requires the user to specify all non-compliance
    distribution parameters.}  
  \item{tcutd0A}{Left hand endpoints for intervals upon which the arm-0 specific \emph{cause-A delay} distribution hazard
    function is constant. The last given component is the left hand endpoint of the interval having right hand endpoint
    infinity. Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{hd0A}{A vector of the same length as \code{tcutd0A} containing peicewise constant hazard rates for the arm-0 
    \emph{cause-A delay} distribution.  Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{sd0A}{When required, the arm-0 \emph{cause-A-delay} distribution is alternately specified via a survival function. A
    vector of the same length as \code{tcutd0A}.}  
  \item{tcutd0B}{Left hand endpoints for intervals upon which the arm-0 specific \emph{cause-B delay} distribution hazard
    function is constant. The last given component is the left hand endpoint of the interval having right hand endpoint
    infinity. Always required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{hd0B}{A vector of the same length as \code{tcutd0B} containing peicewise constant hazard rates for the arm-0
    \emph{cause-B delay} distribution. Always required when \code{noncompliance} is set to any value other than 
    \dQuote{none}.}
  \item{sd0B}{When required, the arm-0 \emph{cause-B-delay} distribution is alternately specified via a
    survival function.  A vector of the same length as \code{tcutd0B}.}
  \item{tcutd1A}{Left hand endpoints for intervals upon which the arm-1 specific \emph{cause-A delay} distribution hazard
    function is constant. The last given component is the left hand endpoint of the interval having right hand endpoint
    infinity. Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{hd1A}{A vector of the same length as \code{tcutd1A} containing peicewise constant hazard rates for the arm-1
    \emph{cause-A delay} distribution. Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{sd1A}{When required, the arm-1 \emph{cause-A-delay} distribution is alternately specified via a survival function. A
    vector of the same length as \code{tcutd1A}.}
  \item{tcutd1B}{Left hand endpoints for intervals upon which the arm-1 specific \emph{cause-B delay} distribution hazard
    function is constant. The last given component is the left hand endpoint of the interval having right hand endpoint
    infinity. Always required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{hd1B}{A vector of the same length as \code{tcutd1B} containing peicewise constant hazard rates for the arm-1
    \emph{cause-B delay} distribution. Always required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{sd1B}{When required, the arm-1 \emph{cause-A-delay} distribution is alternately specified via a survival function. A
    vector of the same length as \code{tcutd1A}.}
  \item{tcutx0A}{Left hand endpoints for intervals upon which the arm-0 specific \emph{post-cause-A-delay-mortality} rate is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity. Required
    only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{hx0A}{A vector of the same length as \code{tcutx0A} containing the arm-0 \emph{post-cause-A-delay mortality} rates.
    Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{sx0A}{When required, the arm-0 \emph{post-cause-A-delay mortality} distribution is alternately specified via a
    survival function. A vector of the same length as \code{tcutx0A}.}
  \item{tcutx0B}{Left hand endpoints for intervals upon which the arm-0 specific \emph{post-cause-B-delay-mortality} rate is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity. Always
    required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{hx0B}{A vector of the same length as \code{tcutx0B} containing the arm-0 \emph{post-cause-B-delay mortality}
    rates. Always required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{sx0B}{When required, the arm-0 \emph{post-cause-B-delay mortality} distribution is alternately specified via a
    survival function.  A vector of the same length as \code{tcutx0B}.}
  \item{tcutx1A}{Left hand endpoints for intervals upon which the arm-1 specific \emph{post-cause-A-delay-mortality} rate is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity. Required
    only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{hx1A}{A vector of the same length as \code{tcutx1A} containing the arm-1 \emph{post-cause-A-delay mortality}
    rates. Required only when \code{noncompliance} is set to \dQuote{mixed} or \dQuote{user}.}
  \item{sx1A}{When required, the arm-1 \emph{post-cause-A-delay mortality} distribution is alternately specified via a
    survival function. A vector of the same length as \code{tcutx1A}.}
  \item{tcutx1B}{Left hand endpoints for intervals upon which the arm-1 specific \emph{post-cause-B-delay-mortality} rate is
    constant. The last given component is the left hand endpoint of the interval having right hand endpoint infinity. Always
    required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{hx1B}{A vector of the same length as \code{tcutx1B} containing the arm-1 \emph{post-cause-B-delay mortality}
    rates. Always required when \code{noncompliance} is set to any value other than \dQuote{none}.}
  \item{sx1B}{When required, the arm-1 \emph{post-cause-B-delay mortality} distribution is alternately specified via a
    survival function. A vector of the same length as \code{tcutx1B}.}
  \item{gradual}{Should the conversion to post-noncompliance mortality be gradual. Under the default behavior,
    \code{gradual}=\code{FALSE}, there is an immediate conversion to the post-noncompliance mortality rate function. If
    \code{gradual} is set to \code{TRUE} then this conversion is done \dQuote{gradually}. In truth, at the individual level what
    is done is that the new mortality rate function is a convex combination of the pre-noncompliance mortality and the
    post-noncompliance mortality, with the weighting in proportion to the time spent in compliance with the study arm protocal.}
  \item{WtFun}{Specifies the name of a weighting function (of time) for assigning relative weights to events according to the
    times at which they occur. The default, \dQuote{FH}, a two parameter weight function, specifies the
    \sQuote{Fleming-Harrington} \code{g-rho} family of weighting functions defined as the pooled arm survival function
    (Kaplan-Meier estimate) raised to the \code{g} times its complement raised to the \code{rho}. Note that
    \code{g}=\code{rho}=0 corresponds to the unweighted log-rank statistic. A second choice is the \dQuote{SFH} function, (for
    \sQuote{Stopped Fleming-Harrington}), meaning that the \dQuote{FH} weights are capped at their value at a user specified
    time, which has a total of 3 parameters.  A third choice is \code{Ramp(tcut)}. Under this choice, weights are assigned in a
    linearly manner from time 0 until a user specified cut-off time, \code{tcut}, after which events are weighted equally. It is
    possible to conduct computations on \code{nstat} candidate statistics within a single run. In this case, \code{WtFun} should
    be a character vector of length \code{nstat} having components set from among the available choices.}
  \item{ppar}{A vector containing all the weight function parameters, in the order determined by that of \dQuote{WtFun}.  For
    example, if \code{WtFun} is set to \code{c("FH","SFH","Ramp")} then \code{ppar} should be a vector of length six, with the
    \dQuote{FH} parameters in the first two elements, \dQuote{SFH} parameters in the next 3 elements, and \dQuote{Ramp}
    parameter in the last element.}
  \item{RR.Futility}{The relative risk corresponding to the alternative alternative hypothesis that is required in the
    construction of the futility boundary. Required if \code{Boundary.Futility} is set to a non-null value.}
  \item{Spend.Info}{When the test statistic is something other than the unweighted log-rank statistic, the variance
    information, i.e. the ratio of variance at interim analysis to variance at the end of trial, is something other than the
    ratio of events at interim analysis to the events at trial end.  The problem is that in practice one doesn't necessarily
    have a good idea what the end of trial variance should be.  In this case the user may wish to spend the type I and type II
    error probabilities according to a different time scale. Possible choices are \dQuote{Variance}, (default), which just uses
    the variance ratio scale, \dQuote{Events}, which uses the events ratio scale, \dQuote{Hybrid(k)}, which makes a linear
    transition from the \dQuote{Variance} scale to the \dQuote{Events} scale beginning with analysis number \code{k}.  The last
    choice, \dQuote{Calendar}, uses the calendar time scale}
  \item{qProp.one.or.Q}{If a futility boundary is specified, what assumption should be made about the drift function (the mean
    value of the weighted log-rank statistic at analysis \code{j} normalized by the square root of the variance function at
    analysis \code{k}).  In practice we don't presume to know the shape of the drift function. Set to \dQuote{one} or
    \dQuote{Q}.  The choice \dQuote{one} results in a more conservative boundary.}
  \item{Nsim}{If you specify \code{method}==\dQuote{S}, then you must specify the number of simulations.  1000 should be
    sufficient.}
  \item{detail}{If you specify \code{method}==\dQuote{S}, and want to see the full level of detail regarding arguments
    returned from the C level code, specify \code{detail}==TRUE}
  \item{StatType}{If you specify \code{method}==\dQuote{S}, then the available choices are \dQuote{WLR} (weighted log-rank)
    and \dQuote{ISD} (integrated survival difference).}  
}

\value{ Returns a value of class \code{PwrGSD} which has components listed below.  Note that the print method will display a
  summary table of estimated powers and type I errors as a \code{nstat} by 2 matrix.  The summary method returns the same
  object invisibly, but after computing the summary table mentioned above, and it is included in the returned value as a
  commponent \code{TBL}.  See examples below.

  \item{dPower}{A \code{length(tlook)} by \code{nstat} matrix containing in each column, an increment in power that resulted
    at that analysis time for the given statistic.}
  \item{dErrorI}{A \code{length(tlook)} by \code{nstat} matrix containing in each column, an increment in type I error that
    resulted at that analysis time for the given statistic.  Always sums to the total alpha specified in \code{alphatot}}
  \item{detail}{A list with components equal to the arguments of the C-call, which correspond in a natural way to the
    arguments specified in the R call, along with the computed results in \code{palpha0vec}, \code{palpha1vec}, \code{inffrac},
    and \code{mu}.  The first two are identical to \code{dErrorI} and \code{dPower}, explained above.  The last two are
    \code{length(tlook)} by \code{nstat} matrices. For each statistic specified in \code{par}, the corresponding columns of
    \code{pinffrac} and \code{mu} contain the information fraction and drift at each of the analysis times.}
  \item{call}{the call} 
} 

\references{ Gu, M.-G. and Lai, T.-L.  \dQuote{Determination of power and sample size in the design
    of clinical trials with failure-time endpoints and interim analyses.}  Controlled Clinical Trials 20 (5): 423-438. 1999
  
  Izmirlian, G.  \dQuote{The PwrGSD package.}  NCI Div. of Cancer Prevention Technical 
  Report. 2004

  Jennison, C. and Turnbull, B.W. (1999)
  Group Sequential Methods: Applications to Clinical Trials
  Chapman & Hall/Crc, Boca Raton FL

  Proschan, M.A., Lan, K.K.G., Wittes, J.T.
  (2006), corr 2nd printing (2008)
  Statistical Monitoring of Clinical Trials
  A Unified Approach
  Springer Verlag, New York
}
\author{Grant Izmirlian \email{izmirlian@nih.gov}}

\seealso{\code{\link{cpd.PwrGSD}}}
\examples{
library(PwrGSD)

test.example <-
  PwrGSD(EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming),
         FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming),
	 RR.Futility = 0.82, sided="1<",method="A",accru =7.73, accrat =9818.65,
         tlook =c(7.14, 8.14, 9.14, 10.14, 10.64, 11.15, 12.14, 13.14,
                  14.14, 15.14, 16.14, 17.14, 18.14, 19.14, 20.14),
         tcut0 =0:19, h0 =c(rep(3.73e-04, 2), rep(7.45e-04, 3),
                            rep(1.49e-03, 15)),
         tcut1 =0:19, rhaz =c(1, 0.9125, 0.8688, 0.7814, 0.6941,
                              0.6943, 0.6072, 0.5202, 0.4332, 0.6520,
                              0.6524, 0.6527, 0.6530, 0.6534, 0.6537,
                              0.6541, 0.6544, 0.6547, 0.6551, 0.6554),
         tcutc0 =0:19, hc0 =c(rep(1.05e-02, 2), rep(2.09e-02, 3),
                              rep(4.19e-02, 15)),
         tcutc1 =0:19, hc1 =c(rep(1.05e-02, 2), rep(2.09e-02, 3),
                              rep(4.19e-02, 15)),
         tcutd0B =c(0, 13), hd0B =c(0.04777, 0),
         tcutd1B =0:6, hd1B =c(0.1109, 0.1381, 0.1485, 0.1637, 0.2446,
                               0.2497, 0),
         noncompliance =crossover, gradual =TRUE,
         WtFun =c("FH", "SFH", "Ramp"),
         ppar =c(0, 1, 0, 1, 10, 10))
}
\keyword{survival}
\keyword{htest}
\keyword{design}
