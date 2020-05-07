##' Wald test for differential amplitute
##'
##' Test differential amplitute of circadian curve fitting using Wald test
##' @title WaldTest_diff_amp
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{amplitute estimate of the 1st data}
##' \item{amp_2}{amplitute estimate of the 2nd data}
##' \item{amp_c}{amplitute estimate pooling all data together}
##' \item{l0}{log likelihood under the null (same variance between the two groups)}
##' \item{l1}{log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{degree of freedom for the LR test}
##' \item{stat}{the LR statistics}
##' \item{pvalue}{the p-value from the LR test}
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
##' WaldTest_diff_amp(tt1, yy1, tt2, yy2)


WaldTest_diff_amp <- function(tt1, yy1, tt2, yy2, period = 24){
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))

	w <- 2*pi/period
	## fit Ha, individual curves
  par1 <- fitSinCurve(tt1,yy1)
	sum_diffy_1_sq <- par1$rss
	theta1 <- n1/sum_diffy_1_sq

  par2 <- fitSinCurve(tt2,yy2)
	sum_diffy_2_sq <- par2$rss
	theta2 <- n2/sum_diffy_2_sq
			  
  beta0 <- c((par1$amp + par2$amp)/2, 
                par1$phase, par1$offset, theta1, 
                par2$phase, par2$offset, theta2
              )
	this_opt_commonAmp <- opt_commonAmp(tt1, yy1, tt2, yy2, period=period, parStart=beta0)
	
  beta_ha <- c(par1$amp, 
                par1$phase, 
								par1$offset, 
								theta1, 
								par2$amp,
								par2$phase, 
								par2$offset, 
								theta2
              )
	#
  beta_h0 <- c(this_opt_commonAmp[1], 
                this_opt_commonAmp[2], 
								this_opt_commonAmp[3], 
								this_opt_commonAmp[4], 
								this_opt_commonAmp[1],
								this_opt_commonAmp[5], 
								this_opt_commonAmp[6], 
								this_opt_commonAmp[7]
              )
  #
	I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)	
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)	
  stat <- t(beta_diff) %*% I8 %*% beta_diff
	
	if(F){
		beta_diff2 <- beta_diff[c(1,5)]
	  I2 <- solve(solve(I8)[c(1,5),c(1,5)])
	  stat <- t(beta_diff2) %*% I2 %*% beta_diff2		
	}
	   
  dfdiff <- 1
  pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  
  res <- list(amp_c=this_opt_commonAmp[1], amp_1=par1$amp, amp_2=par2$amp,
	  df = dfdiff, 
	  stat = stat, 
	  pvalue = pvalue)
  return(res)
}
