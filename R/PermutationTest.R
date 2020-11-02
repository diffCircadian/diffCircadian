##' Permutation test sin function
##'
##' Test the signficance of circadian curve fitting using Permutation test
##' @title Permutation test
##' @param tt time vector
##' @param yy expression vector
##' @param B number of permutations.Default is 1000.
##' @param period Period of the since curve. Default is 24.
##' @return A list of amp, phase, offset, peak, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{phase based on formula 1, phase is restricted within (0, period)}
##' \item{offset}{offset based on formula 1 or on formula 2}
##' \item{tss}{Total sum of square}
##' \item{rss}{residual sum of square, SSE/n is the MLE of the variance sigma2}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' \item{pvalue}{p-value of the permutation test}
##' \item{R2b}{Sequence of R squares for each permutaion}
##' @author Haocheng Ding
##' @export
##' @examples
##' set.seed(32608)
##' n <- 10
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' permTest(tt, yy)



permTest <- function(tt,yy,B=1000,period=24){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss
  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset
  
  R2 <- 1-rss/tss
  R2b <- rep(0,B)
  for(i in 1:B){
    tti <- sample(tt)
    rssi <- fitSinCurve(tti,yy,period)$rss
    R2b[i] <- 1-rssi/tss
  }
  return(list(amp=amp,phase=phase,offset=offset,R2=R2,rss=rss,tss=tss,pvalue=sum(R2b>=R2)/B,R2b=R2b))
}
#