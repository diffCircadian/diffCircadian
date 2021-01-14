##' F test for sin function.
##'
##' Test the significance of circadian curve fitting using F test.
##' @title F test for detecting circadian pattern
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @return A list of amp, phase, offset, peak, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
##' \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
##' \item{tss}{Total sum of square}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma2}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' @author Caleb

##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' FTest(tt, yy)

FTest <- function(tt,yy, period = 24){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss
  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset
    
  fss <-tss-rss
  
  df1 <- 2
  df2 <- n-3
  fvalue <- fss/df1/(rss/df2)
  R2 <- 1-rss/tss
  pvalue <- pf(fvalue,df1,df2,lower.tail = F)
  
  return(list(amp=amp, phase=phase, offset=offset, rss=rss, tss=tss,R2=R2 ,df1 = df1, df2= df2, stat=fvalue,pvalue=pvalue))
}

