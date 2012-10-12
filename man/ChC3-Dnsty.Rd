\name{Dnsty}
\alias{Dnsty}
\title{Preliminary version of function to help compute design adjusted p-values}
\description{Preliminary version of function to help compute design adjusted p-values}
\usage{
Dnsty(EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming),
      FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming),
      NonBindingFutility = TRUE, frac, frac.ii = NULL, drift = NULL)}
\arguments{
  \item{EfficacyBoundary}{This specifies the method used to construct the efficacy boundary. See the
    \code{\link{GrpSeqBnds}} documentation for details}
  \item{FutilityBoundary}{This specifies the method used to construct the futility boundary. See the
    \code{\link{GrpSeqBnds}} documentation for details}
  \item{NonBindingFutility}{If using a futility boundary set this to \code{TRUE} if a non-binding
    futility boundary is desired.  See the \code{\link{GrpSeqBnds}} documentation for details}
  \item{frac}{The variance ratio. See the \code{\link{GrpSeqBnds}} documentation for details}
  \item{frac.ii}{The second information scale that is used for type I and type II error probability
    spending. See the \code{\link{GrpSeqBnds}} documentation for details}
  \item{drift}{If using a futility boundary you must specify the drift function which is the expected
    value of the statistic, normalized to the brownian scale, at time \code{t} as time varies, under
    a specified alternative.}
  }
\author{Grant Izmirlian}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming), 
    FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming), 
    NonBindingFutility = TRUE, frac, frac.ii = NULL, drift = NULL) 
{
    .call. <- match.call()
    xc.km1 <- dF.km1 <- xc.k <- xr.k <- NULL
    is.frac.ii <- TRUE
    if (missing(frac.ii)) {
        frac.ii <- frac
        is.frac.ii <- FALSE
    }
    if (is.frac.ii && length(frac.ii) != length(frac)) 
        stop("Lengths of 'frac' and 'frac.ii' must agree")
    normut <- 8.20953615160139
    psimin <- 1e-15
    nlooks <- length(frac)
    call.PBS <- .call.
    call.PBS[[1]] <- as.name("ParseBoundarySelection")
    call.PBS$frac <- call.PBS$frac.ii <- NULL
    call.PBS$check.drift <- TRUE
    call.PBS$n.looks <- nlooks
    eval.PBS <- eval(call.PBS)
    do.efficacy <- TRUE
    do.futility <- !is.null(eval.PBS$frontend$FutilityBoundary)
    is.drift <- !missing(drift)
    nbf <- do.futility && NonBindingFutility
    if (!do.futility) 
        drift <- rep(0, nlooks)
    Alpha.Efficacy <- eval.PBS$backend$Alpha.Efficacy
    Alpha.Futility <- eval.PBS$backend$Alpha.Futility
    if (is.null(Alpha.Futility)) 
        Alpha.Futility <- 0
    nbnd.e <- eval.PBS$backend$nbnd.e
    nbnd.f <- eval.PBS$backend$nbnd.f
    nsf.e <- eval.PBS$backend$nsf.e
    nsf.f <- eval.PBS$backend$nsf.f
    rho.Efficacy <- eval.PBS$backend$rho.Efficacy
    rho.Futility <- eval.PBS$backend$rho.Futility
    b.Haybittle.e <- eval.PBS$backend$b.Haybittle.e
    b.Haybittle.f <- eval.PBS$backend$b.Haybittle.f
    drift.end <- eval.PBS$backend$drift.end
    be.end <- eval.PBS$backend$be.end
    prob.e <- eval.PBS$backend$prob.e
    prob.f <- eval.PBS$backend$prob.f
    my.Efficacy <- eval.PBS$backend$my.Efficacy
    my.Futility <- eval.PBS$backend$my.Futility
    is.myE <- !all(my.Efficacy == 0)
    is.myF <- !all(my.Futility == 0)
    nunq <- nlooks
    glegx <- glegx24
    glegw <- glegw24
    ngqnodes <- length(glegx)
    b.e <- alpha.e <- rep(0, nlooks)
    b.f <- alpha.f <- rep(0, nlooks)
    l <- 1
    l.act.e <- 1
    l.act.f <- 1
    fracold <- rep(0, 2)
    fracold.ii <- rep(0, 2)
    fracnew <- frac[1]
    fracnew.ii <- frac.ii[1]
    sc.drift.factor <- 1
    mu.end <- drift.end
    if (nbnd.e[1] != 2) 
        bold.e <- normut
    if (nbnd.e[1] == 2) 
        bold.e <- b.Haybittle.e[1]
    if (nbnd.f[1] != 2) 
        bold.f <- -normut
    if (nbnd.f[1] == 2) 
        bold.f <- b.Haybittle.f[1]
    y.e <- tmp.e <- rep(0, ngqnodes)
    y.f <- tmp.f <- rep(0, ngqnodes)
    x.e <- dF.e <- rep(0, ngqnodes)
    x.f <- dF.f <- rep(0, ngqnodes)
    if (nbnd.e[1] == 3) 
        my.Efficacy <- rep(0, nlooks)
    if (nbnd.f[1] == 3) 
        my.Futility <- rep(0, nlooks)
    while (l <= nlooks) {
        mu <- drift[l]
        if (nbnd.e[l] == 3) {
            my.Efficacy[l] <- .C("StCu2Bnds", pmu = double(2), 
                pfrac = as.double(fracnew.ii), pzcrit = as.double(be.end), 
                prho = as.double(prob.e[l]), pef = as.integer(0), 
                b = double(1), PACKAGE = "PwrGSD")$b
            is.myE <- TRUE
            nsf.e[l] <- 1
        }
        if (nbnd.f[l] == 3) {
            my.Futility[l] <- .C("StCu2Bnds", pmu = as.double(c(mu, 
                mu.end) * sc.drift.factor), pfrac = as.double(fracnew.ii), 
                pzcrit = as.double(be.end), prho = as.double(prob.f[l]), 
                pef = as.integer(1), b = double(1), PACKAGE = "PwrGSD")$b
            is.myF <- TRUE
            nsf.f[l] <- 1
        }
        bnew <- c(normut, -normut)
        if (nbnd.e[l] == 3 || nbnd.e[l] == 4) 
            bnew[1] <- my.Efficacy[l]
        if ((nbnd.f[l] == 3 || nbnd.f[l] == 4) && (1 - fracnew.ii >= 
            1e-06)) 
            bnew[2] <- my.Futility[l]
        ans <- .C("dnsty", nbf = as.integer(nbf), nbnd = as.integer(c(nbnd.e[l], 
            nbnd.f[l])), nsf = as.integer(c(nsf.e[l], nsf.f[l])), 
            rho = as.double(c(rho.Efficacy[l], rho.Futility[l])), 
            pnlook = as.integer(c(l.act.e, l.act.f)), palphtot = as.double(c(Alpha.Efficacy, 
                Alpha.Futility)), palpha = double(2), psimin = as.double(psimin), 
            dlact = integer(2), pfracold = as.double(fracold), 
            pfracnew = as.double(fracnew), pfracold.ii = as.double(fracold.ii), 
            pfracnew.ii = as.double(fracnew.ii), xc = as.double(c(x.e, 
                x.f)), y = as.double(c(y.e, y.f)), tmp = as.double(c(tmp.e, 
                tmp.f)), dF = as.double(c(dF.e, dF.f)), gqxw = as.double(c(glegx, 
                glegw)), pngqnodes = as.integer(ngqnodes), pxc = double(ngqnodes), 
            pxr = double(ngqnodes), bold = as.double(c(bold.e, 
                bold.f)), bnew = as.double(bnew), mybounds = as.integer(c(is.myE, 
                is.myF)), PACKAGE = "PwrGSD")
        xc.km1 <- cbind(xc.km1, ans$xc[1:ngqnodes])
        dF.km1 <- cbind(dF.km1, ans$dF[1:ngqnodes])
        xc.k <- cbind(xc.k, ans$pxc[1:ngqnodes])
        xr.k <- cbind(xr.k, ans$pxr[1:ngqnodes])
        dlact <- ans$dlact
        if (dlact[1] == 1) {
            if (nbnd.e[l] == 1 || nbnd.e[l] == 3 || nbnd.e[l] == 
                4) 
                b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], ans$bnew[1])
            if (nbnd.e[l] == 2) {
                if (1 - fracnew.ii >= 1e-06) {
                  b.e[l] <- bold.e <- b.Haybittle.e[l]
                  Alpha.Efficacy <- Alpha.Efficacy - ans$palpha[1]
                }
                else b.e[l] <- ans$bnew[1]
            }
            alpha.e[l] <- ans$palpha[1]
            x.e <- ans$x[1:ngqnodes]
            dF.e <- ans$dF[1:ngqnodes]
            fracold[1] <- ans$pfracnew
            fracold.ii[1] <- ans$pfracnew.ii
            l.act.e <- l.act.e + 1
        }
        else {
            b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], normut)
            alpha.e[l] <- ifelse(is.myE, ans$palpha[1], psimin)
            fracold.ii[1] <- ans$pfracnew.ii
        }
        if (do.futility && dlact[2] == 1) {
            if (nbnd.f[l] == 1 || nbnd.f[l] == 3 || nbnd.f[l] == 
                4) 
                b.f[l] <- bold.f <- ifelse((nbnd.f[l] == 3 || 
                  nbnd.f[l] == 4) && (1 - fracnew.ii >= 1e-06), 
                  bnew[2], ans$bnew[2])
            if (nbnd.f[l] == 2) {
                if (l < nlooks) {
                  b.f[l] <- bold.f <- b.Haybittle.f[l]
                  Alpha.Futility <- Alpha.Futility - ans$palpha[2]
                }
                else b.f[l] <- ans$bnew[2]
            }
            alpha.f[l] <- ans$palpha[2]
            x.f <- ans$x[ngqnodes + (1:ngqnodes)]
            dF.f <- ans$dF[ngqnodes + (1:ngqnodes)]
            fracold[2] <- ans$pfracnew
            fracold.ii[2] <- ans$pfracnew.ii
            l.act.f <- l.act.f + 1
        }
        if (do.futility && dlact[2] == 0) {
            b.f[l] <- bold.f <- ifelse(is.myF || (l < nlooks), 
                bnew[2], -normut)
            alpha.f[l] <- ifelse(is.myF, ans$palpha[2], psimin)
            fracold.ii[2] <- ans$pfracnew.ii
        }
        l <- l + 1
        fracnew <- frac[l]
        fracnew.ii <- frac.ii[l]
        mu <- drift[l]
    }
    out <- frac.ii
    nms <- "frac"
    if (do.futility) {
        out <- cbind(out, b.f, alpha.f, cumsum(alpha.f))
        nms <- c(nms, "b.f", "alpha.f", "cum-alpha.f")
    }
    out <- cbind(out, b.e, alpha.e, cumsum(alpha.e))
    nms <- c(nms, "b.e", "alpha.e", "cum-alpha.e")
    dimnames(out) <- list(1:nlooks, nms)
    ans <- list(table = out, frac = frac, frac.ii = frac.ii, 
        drift = drift, xc.km1 = xc.km1, dF.km1 = dF.km1, xc.k = xc.k, 
        xr.k = xr.k, call = .call.)
    class(ans) <- "boundaries"
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
