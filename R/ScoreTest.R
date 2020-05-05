##' Score test for circadian pattern detection
##'
##' Test the signficance of circadian curve fitting using Score test
##' @title ScoreTest
##' @param tt time vector
##' @param yy expression vector
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
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' ScoreTest(tt, yy)

ScoreTest <- function(tt, yy, period = 24){
  afit <- fitSinCurve(tt, yy)
  n <- length(tt)
  
  ## HA: 
  if(F){
	  A <- afit$A 
	  B <- afit$B 
	  offset <- afit$offset   	
  }
  
  asin <- sin(2*pi/period * tt)
  acos <- cos(2*pi/period * tt)

  #yhat <- A * asin + B * acos + offset
  yhat <- mean(yy)
  residual <- yy - yhat
  sigma02 <- 1/(n)*sum(residual^2)
  invSigmaA0 <- 1/sigma02
  
  #sigmaA2 <- 1/n * sum(residual^2)
  #invSigmaA2 <- 1/sigmaA2
  
  s_A <- invSigmaA0 * sum(asin*residual)
  s_B <- invSigmaA0 * sum(acos*residual)
  s_offset <- invSigmaA0 * sum(residual)
  s_sigma2 <- invSigmaA0^2/2 * sum(residual^2) - n*invSigmaA0/2
  
  det_A_A <- - invSigmaA0 * sum(asin^2)
  det_A_B <- - invSigmaA0 * sum(asin * acos)
  det_A_offset <- - invSigmaA0 * sum(asin)
  det_A_sigma2 <- - invSigmaA0^2 * sum(residual * asin)
  
  det_B_B <- - invSigmaA0 * sum(acos^2) 
  det_B_offset <- - invSigmaA0 * sum(acos) 
  det_B_sigma2 <- - invSigmaA0^2 * sum(residual * acos)
  
  det_offset_offset <- - n * invSigmaA0
  det_offset_sigma2 <- - invSigmaA0^2 * sum(residual) 
  
  det_sigma2_sigma2 <- n/2*invSigmaA0^2 - invSigmaA0^3 * sum(residual^2) 
  
  r1 <- c(det_A_A,det_A_B,det_A_offset,det_A_sigma2)
  r2 <- c(det_A_B,det_B_B,det_B_offset,det_B_sigma2)
  r3 <- c(det_A_offset,det_B_offset,det_offset_offset,det_offset_sigma2)
  r4 <- c(det_A_sigma2,det_B_sigma2,det_offset_sigma2,det_sigma2_sigma2)
  I <- - rbind(r1,r2,r3,r4)
  
  vec_s <- matrix(c(s_A, s_B, s_offset, s_sigma2), nrow = 4)
  Scorestat <- t(vec_s) %*% solve(I) %*% vec_s
  
  vec_u <- matrix(c(s_A, s_B), nrow = 2)
  Scorestat2 <- t(vec_u) %*% solve(I)[1:2,1:2] %*% vec_u
  pvalue2 <- pchisq(Scorestat2,2,lower.tail = F)
    
  df <- 2
  stat <- Scorestat
  pvalue <- pchisq(Scorestat,2,lower.tail = F)
  
  res <- list(
	  invSigmaA0=invSigmaA0,
	  s_A=s_A, s_B=s_B, s_offset=s_offset, s_sigma2=s_sigma2,
	  det_A_A=det_A_A, det_A_B=det_A_B, det_A_offset=det_A_offset, det_A_sigma2=det_A_sigma2,
	  det_B_B=det_B_B, det_B_offset=det_B_offset, det_B_sigma2=det_B_sigma2, 
	  det_offset_offset=det_offset_offset, det_offset_sigma2=det_offset_sigma2, 
	  det_sigma2_sigma2=det_sigma2_sigma2,
	  df=df, stat=stat, pvalue=pvalue, 
	  stat2=Scorestat2, pvalue2=pvalue2 
	  )
  return(res)
  
}





















