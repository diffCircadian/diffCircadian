##' Likelihood ratio test for sigma2
##'
##' Test diff sigma2 of circadian curve fitting using likelihood ratio test
##' @title LRTest_diff_sigma2
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{phase based on formula 1}
##' \item{offset}{offset based on formula 1 or on formula 2}
##' \item{A}{A based on formula 2}
##' \item{B}{B based on formula 2}
##' \item{SST}{Total sum of square}
##' \item{SSE}{Error sum of square}
##' \item{R2}{Pseudo R2 defined as (SST - SSE)/SST}
##' @author Caleb
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24) 
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24) 
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' LRTest_diff_sigma2(tt1, yy1, tt2, yy2)

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

LRTest_diff_sigma2 <- function(tt1, yy1, tt2, yy2, period = 24){
  this_opt_commonSigma <- opt_commonSigma(tt1, yy1, tt2, yy2, period=period)
  
  ## extract parameters	  
  A_1 <- this_opt_commonSigma$A_1
  B_1 <- this_opt_commonSigma$B_1
  sigma2_1 <- this_opt_commonSigma$sigma2_1
  offset_1 <- this_opt_commonSigma$offset_1
  n1 <- this_opt_commonSigma$n1
  
  A_2 <- this_opt_commonSigma$A_2
  B_2 <- this_opt_commonSigma$B_2
  sigma2_2 <- this_opt_commonSigma$sigma2_2
  offset_2 <- this_opt_commonSigma$offset_2
  n2 <- this_opt_commonSigma$n2
  
  sigma2_C <- this_opt_commonSigma$sigma2_C
  
  asin_1 <- sin(2*pi/period * tt1)
  acos_1 <- cos(2*pi/period * tt1)
  yhat_1 <- A_1 * asin_1 + B_1 * acos_1 + offset_1
  residual_1 <- yy1 - yhat_1
  
  asin_2 <- sin(2*pi/period * tt2)
  acos_2 <- cos(2*pi/period * tt2)
  yhat_2 <- A_2 * asin_2 + B_2 * acos_2 + offset_2
  residual_2 <- yy2 - yhat_2
    
	#1/(2*sigma2_C)*sum(residual_1^2) + 1/(2*sigma2_C)*sum(residual_2^2)
	#1/(2*sigma2_1)*sum(residual_1^2) + 1/(2*sigma2_2)*sum(residual_2^2)

  ## H0 (sigma2_C)
  l0_1 <- -n1/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_1^2)
  l0_2 <- -n2/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_2^2)
  l0 <- l0_1 + l0_2
  
  ## H1 (sigma2_1, sigma2_2)
  la_1 <- -n1/2*log(2*pi*sigma2_1) # - 1/(2*sigma2_1)*sum(residual_1^2)
  la_2 <- -n2/2*log(2*pi*sigma2_2) # - 1/(2*sigma2_2)*sum(residual_2^2)
  la <- la_1 + la_2
    
  dfdiff <- 1
  pvalue <- pchisq(-2*(l0-la),dfdiff,lower.tail = F)
  
  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_C=sigma2_C,
	  l0=l0, 
	  la=la, 
	  df = dfdiff, 
	  stat=-2*(l0-la), 
	  pvalue=pvalue)
  return(res)
}