##' Finite sample/ large sample Likelihood ratio test for circadian pattern detection
##'
##' Test the significance of circadian curve fitting using finite sample likelihood ratio test
##' @title LR Test for detecting circadian pattern.
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Test, finite sample "FN" or large sample "LS", default is "FN".
##' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
##' \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
##' \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
##' \item{sigma02}{Variance estimate under the null (intercept only)}
##' \item{sigmaA2}{Variance estimate under the alternative (since curve fitting)}
##' \item{l0}{Log likelihood under the null (intercept only)}
##' \item{l1}{Log likelihood under the alternative (since curve fitting)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
##' @author Caleb
##' @noRd
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' LRTest(tt, yy)

#model: y=A*sin(2*pi*x+B)+C
#y: a 1*n vector of data y
#A: estimated A^hat from fitCurve
#B: estimated B^hat from fitCurve
#C: estimated C^hat from fitCurve
#sigma0: sigma0^hat under H0
#sigmaA: sigmaA^hat under H1
#n: length of data y
#df0: df under H0
#df1: df under H1

LRTest <- function(tt,yy, period = 24,type=TRUE){
  fitCurveOut <- fitSinCurve(tt,yy,period=period)
  n <- length(yy)
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss

  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset
  
  sigma02 <- 1/(n)*sum((yy-mean(yy))^2)
  sigmaA2 <- 1/(n)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)
    
  l0 <- -n/2*log(2*pi*sigma02)-1/(2*sigma02)*sum((yy-mean(yy))^2)
  l1 <- -n/2*log(2*pi*sigmaA2)-1/(2*sigmaA2)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)
  
  dfdiff <- (n-1)-(n-3)
	LR_stat <- -2*(l0-l1)
	
	if(type==FALSE){
	  pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
	}
  else if(type==TRUE){
    r <- 2
    k <- 3
    LR_stat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(LR_stat,df1 = r, df2 = n-k, lower.tail = F)
  }
	  R2 <- 1-rss/tss
    res <- list(
	  amp = amp,
	  phase = phase,
		peakTime = (6 - phase) %% period, 
	  offset = offset, 
	  sigma02=sigma02, sigmaA2=sigmaA2, 
	  l0=l0, 
	  l1=l1, 
	  stat=LR_stat, 
	  pvalue=pvalue,R2=R2)
  return(res)
}