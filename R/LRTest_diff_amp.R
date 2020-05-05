##' Likelihood ratio test for amplitute
##'
##' Test diff amp of circadian curve fitting using likelihood ratio test
##' @title LRTest_diff_amp
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
##' LRTest_diff_amp(tt1, yy1, tt2, yy2)

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


LRTest_diff_amp <- function(tt1, yy1, tt2, yy2, period = 24){
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))

	w <- 2*pi/period
	## fit Ha, individual curves
  par1 <- fitSinCurve(tt1,yy1)
	sum_diffy_1_sq <- par1$SSE
	theta1 <- n1/sum_diffy_1_sq

  par2 <- fitSinCurve(tt2,yy2)
	sum_diffy_2_sq <- par2$SSE
	theta2 <- n2/sum_diffy_2_sq
		
	l1_Ha <- 1/2 * n1 * log(theta1) - 1/2 * n1
	l2_Ha <- 1/2 * n2 * log(theta2) - 1/2 * n2
	
	la <- unlist(l1_Ha + l2_Ha)
	  
  beta0 <- c((par1$amp + par2$amp)/2, 
                par1$phase, par1$offset, theta1, 
                par2$phase, par2$offset, theta2
              )
	

  this_opt_commonAmp <- opt_commonAmp(tt1, yy1, tt2, yy2, period=period, parStart=beta0)
  
	## fit H0, individual curves
  ## extract parameters	  
  amp_c <- this_opt_commonAmp[1]
  phase_1 <- this_opt_commonAmp[2]
  offset_1 <- this_opt_commonAmp[3]
  theta_1 <- this_opt_commonAmp[4]
  phase_2 <- this_opt_commonAmp[5] 
  offset_2 <- this_opt_commonAmp[6]
  theta_2 <- this_opt_commonAmp[7] 

  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))
  
  yhat_1 <- amp_c * asin_1 + offset_1
  yhat_2 <- amp_c * asin_2 + offset_2
  
  diffy_1 <- yy1-yhat_1
  diffy_2 <- yy2-yhat_2
  
  l1a <- 1/2 * n1 * log(theta_1)
  l1b <- - 1/2 * sum(diffy_1^2) * theta_1
  
  l2a <- 1/2 * n2 * log(theta_2)
  l2b <- - 1/2 * sum(diffy_2^2) * theta_2
  
  l0 <- unlist(l1a + l1b + l2a + l2b)
     
  dfdiff <- 1
  pvalue <- pchisq(-2*(l0-la),dfdiff,lower.tail = F)
  
  res <- list(amp_c=amp_c, amp_1=par1$amp, amp_2=par2$amp,
	  l0=l0, 
	  la=la, 
	  df = dfdiff, 
	  stat=-2*(l0-la), 
	  pvalue=pvalue)
  return(res)
}