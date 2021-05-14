##' Finite sample/Large sample Wald test for differential sigma square.
##'
##' Test differential sigma square of circadian curve fitting using Wald test
##' @title Wald test for detecting differential sigma square
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{sigma2_1}{Variance estimate of the 1st data}
##' \item{sigma2_2}{Variance estimate of the 2nd data}
##' \item{sigma2_C}{Variance estimate pooling all data together}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
##' @author Caleb
##' @noRd
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
##' WaldTest_diff_sigma2(tt1, yy1, tt2, yy2)


WaldTest_diff_sigma2 <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	#period <- 24
	w <- 2*pi/period
	
	fit1 <- fitSinCurve(tt1, yy1, period = 24)
	fit2 <- fitSinCurve(tt2, yy2, period = 24)
	
	A1 <- fit1$amp
	A2 <- fit2$amp
	
	phase1 <- fit1$phase
	phase2 <- fit2$phase
	
	basal1 <- fit1$offset
	basal2 <- fit2$offset
	
  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss	
  sigma2_c <- 1/(n1 + n2) * (sigma2_1 * n1 + sigma2_2 * n2)
	
	theta1 <- 1/sigma2_1
	theta2 <- 1/sigma2_2
	thetaC <- 1/sigma2_c
	
	p1 <- c(A1, phase1, basal1, theta1)
	p2 <- c(A2, phase2, basal2, theta2)
	
	x_Ha <- c(p1, p2)
	x_H0 <- c(A1, phase1, basal1, thetaC, A2, phase2, basal2, thetaC)

  beta_ha <- x_Ha
  beta_h0 <- x_H0

  
  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)	
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)	
  #stat <- t(beta_diff) %*% I8 %*% beta_diff
  
  beta_diff2 <- beta_diff[c(4,8)]
  I2 <- solve(solve(I8)[c(4,8),c(4,8)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)
  
  dfdiff <- 1
  
  if(FN==FALSE){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }
  
  
  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_c=sigma2_c,
              stat=stat, 
              pvalue=pvalue)
  return(res)
}

