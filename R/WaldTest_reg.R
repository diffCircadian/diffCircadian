##' Large Sample Regular Wald test for circadian pattern detection
##'
##' Test the signficance of circadian curve fitting using Wald test
##' @title WaldTest
##' @param tt time vector
##' @param yy expression vector
##' @param period Period of the since curve. Default is 24.
##' @return A list of A, B, offset, df, stat, and pvalue
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{A}{A based on formula 2}
##' \item{B}{B based on formula 2}
##' \item{offset}{offset based on formula 1 or on formula 2}
##' \item{df}{degree of freedom for the Wald test}
##' \item{stat}{the Wald statistics}
##' \item{pvalue}{the p-value from the Wald test}
##' @author Caleb
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' WaldTest_reg(tt, yy)

WaldTest_reg <- function(tt, yy, period = 24){
  afit <- fitSinCurve(tt, yy)
  n <- length(tt)
  
  ## HA: 
  A <- afit$A 
  B <- afit$B 
  offset <- afit$offset 
  
  asin <- sin(2*pi/period * tt)
  acos <- cos(2*pi/period * tt)

  yhat <- A * asin + B * acos + offset
  residual <- yy - yhat
  
  sigmaA2 <- 1/n * sum(residual^2)
  invSigmaA2 <- 1/sigmaA2
  
  
  det_A_A <- - invSigmaA2 * sum(asin^2)
  det_A_B <- - invSigmaA2 * sum(asin * acos)
  det_A_offset <- - invSigmaA2 * sum(asin)
  det_A_sigma2 <- - invSigmaA2^2 * sum(residual * asin)
  
  det_B_B <- - invSigmaA2 * sum(acos^2) 
  det_B_offset <- - invSigmaA2 * sum(acos) 
  det_B_sigma2 <- - invSigmaA2^2 * sum(residual * acos)
  
  det_offset_offset <- - n * invSigmaA2
  det_offset_sigma2 <- - invSigmaA2^2 * sum(residual) 
  
  det_sigma2_sigma2 <- n/2*invSigmaA2^2 - invSigmaA2^3 * sum(residual^2) 
  
  r1 <- c(det_A_A,det_A_B,det_A_offset,det_A_sigma2)
  r2 <- c(det_A_B,det_B_B,det_B_offset,det_B_sigma2)
  r3 <- c(det_A_offset,det_B_offset,det_offset_offset,det_offset_sigma2)
  r4 <- c(det_A_sigma2,det_B_sigma2,det_offset_sigma2,det_sigma2_sigma2)
  I <- - rbind(r1,r2,r3,r4)
  
  I_test <- solve(solve(I)[1:2,1:2])
  Waldstat <- matrix(c(A, B), nrow = 1, ncol = 2) %*% 
    I_test %*% matrix(c(A, B), nrow = 2, ncol = 1)		
  
  df <- 2
  stat <- as.numeric(Waldstat)
  pvalue <- pchisq(stat,2,lower.tail = F)
  
	if(F){
		## for internal test purpose
	  res <- list(
		  A=A,B=B,offset=offset,invSigmaA2=invSigmaA2,
		  det_A_A=det_A_A, det_A_B=det_A_B, det_A_offset=det_A_offset, det_A_sigma2=det_A_sigma2,
		  det_B_B=det_B_B, det_B_offset=det_B_offset, det_B_sigma2=det_B_sigma2, 
		  det_offset_offset=det_offset_offset, det_offset_sigma2=det_offset_sigma2, 
		  det_sigma2_sigma2=det_sigma2_sigma2,
		  df=df, stat=stat, pvalue=pvalue)
	}
  res <- list(
	  A=A,B=B,offset=offset,
	  df=df, stat=stat, pvalue=pvalue
		)

  return(res)
  
}





















