##' Optimization for common offset exist
##'
##' Fit a joint sin curve for two conditions assuming common amplitute
##' @title fit the data based on joint sin curve assuming common amplitute
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param parStart initial value for optimzation purpose
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
##' @import minpack.lm
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
##' Phase2 <- 15
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' opt_commonOffset(tt1, yy1, tt2, yy2)


opt_commonOffset <- function(tt1, yy1, tt2, yy2, period = 24, 
                          parStart = NULL){
	
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	if(is.null(parStart)){
	  par1 <- fitSinCurve(tt1,yy1)
	  yhat1 <- par1$amp * sin(2*pi/24 * (tt1 + par1$phase)) + par1$offset
	  #theta1 <- 1/mean((yy1-yhat1)^2)
	  theta1 <- n1/par1$SSE
	  
	  par2 <- fitSinCurve(tt2,yy2)
	  yhat2 <- par2$amp * sin(2*pi/24 * (tt2 + par2$phase)) + par2$offset
	  #theta2 <- 1/mean((yy2-yhat2)^2)
	  theta2 <- n2/par2$SSE
	  
	  beta0 <- c(par1$amp, 
	             par1$phase, 
	             (par1$offset + par2$offset)/2,
	             theta1, 
	             par1$amp, 
	             par2$phase, 
	             theta2
	              )
	} else {
		beta0 <- parStart
	}
	
	w <- 2*pi/period
	
	# beta[1]: amp_1
	# beta[2]: phase_1
	# beta[3]: offset_c
	# beta[4]: theta_1
	# beta[5]: amp_2
	# beta[6]: phase_2
	# beta[7]: theta_2
	
	f <- function(beta){
	  amp_1 <- beta[1]
	  phase_1 <- beta[2]
	  offset_c <- beta[3]
	  theta_1 <- beta[4]
	  amp_2 <- beta[5]
	  phase_2 <- beta[6] 
	  theta_2 <- beta[7] 

	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_1 * asin_1 + offset_c
	  yhat_2 <- amp_2 * asin_2 + offset_c
	  
	  diffy_1 <- yy1-yhat_1
	  diffy_2 <- yy2-yhat_2
	  
	  l1a <- - 1/2 * n1 * log(theta_1)
	  l1b <- 1/2 * sum(diffy_1^2) * theta_1
	  
	  l2a <- - 1/2 * n2 * log(theta_2)
	  l2b <- 1/2 * sum(diffy_2^2) * theta_2
	  
	  ll <- l1a + l1b + l2a + l2b
	  as.numeric(ll)
	}
	#f(beta0)
	#optimx(beta0, fn=f)
	
	#beta <- beta0
	g <- function(beta){
	  amp_1 <- beta[1]
	  phase_1 <- beta[2]
	  offset_c <- beta[3]
	  theta_1 <- beta[4]
	  amp_2 <- beta[5]
	  phase_2 <- beta[6] 
	  theta_2 <- beta[7] 
	  
	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  acos_1 <- cos(w * (tt1 + phase_1))
	  acos_2 <- cos(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_1 * asin_1 + offset_c
	  yhat_2 <- amp_2 * asin_2 + offset_c
	  
	  diffy_1 <- yy1-yhat_1
	  diffy_2 <- yy2-yhat_2
	  
	  partial_amp_c_1 <- - theta_1 * sum(diffy_1 * asin_1) 
	  partial_phase_1 <- - amp_1 * w * theta_1 * sum(diffy_1 * acos_1) 
	  partial_offset_c <- - theta_1 * sum(diffy_1) - theta_2 * sum(diffy_2) 
	  partial_theta_1 <- - n1 / 2 / theta_1 + sum(diffy_1^2)/ 2 
	  
	  partial_amp_c_2 <- - theta_2 * sum(diffy_2 * asin_2) 
	  partial_phase_2 <- - amp_2 * w * theta_2 * sum(diffy_2 * acos_2) 
	  partial_theta_2 <- - n2 / 2 / theta_2 + sum(diffy_2^2)/ 2 
	  
	  c(partial_amp_c_1, 
	    partial_phase_1, 
	    partial_offset_c, 
	    partial_theta_1, 
	    partial_amp_c_2, 
	    partial_phase_2, 
	    partial_theta_2)
	}
	#g(beta0)
	#optimx(beta0, fn=f, gr = g, control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	
	h <- function(beta){
	  amp_1 <- beta[1]
	  phase_1 <- beta[2]
	  offset_c <- beta[3]
	  theta_1 <- beta[4]
	  amp_2 <- beta[5]
	  phase_2 <- beta[6] 
	  theta_2 <- beta[7] 
	  
	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  acos_1 <- cos(w * (tt1 + phase_1))
	  acos_2 <- cos(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_1 * asin_1 + offset_c
	  yhat_2 <- amp_2 * asin_2 + offset_c
	  
	  diffy_1 <- yy1-yhat_1
	  diffy_2 <- yy2-yhat_2
	  
	  h_amp1_amp1 <- theta_1 * sum(asin_1^2) 
	  h_amp1_phase1 <- theta_1 * w * sum((amp_1 * asin_1 - diffy_1) * acos_1)
	  h_amp1_offsetc <- theta_1 * sum(asin_1)
	  h_amp1_theta1 <- - sum(diffy_1 * asin_1)

	  h_phase1_phase1 <- theta_1 * amp_1 * w^2 * sum(amp_1 * acos_1^2 + diffy_1 * asin_1)
	  h_phase1_offsetc <- theta_1 * amp_1 * w * sum(acos_1)
	  h_phase1_theta1 <- - amp_1 * w * sum(diffy_1 * acos_1)
	  
	  h_offsetc_offsetc <- theta_1 * n1 + theta_2 * n2
	  h_offsetc_theta1 <- - sum(diffy_1)
	  h_offsetc_amp2 <- theta_2 * sum(asin_2)
	  h_offsetc_phase2 <- theta_2 * amp_2 * w * sum(acos_2)
	  h_offsetc_theta2 <- - sum(diffy_2)
	  
	  h_theta1_theta1 <- n1 / 2 / theta_1^2
	  
	  h_amp2_amp2 <- theta_2 * sum(asin_2^2) 
	  h_amp2_phase2 <- theta_2 * w * sum((amp_2 * asin_2 - diffy_2) * acos_2)
	  h_amp2_theta2 <- - sum(diffy_2 * asin_2)

	  h_phase2_phase2 <- theta_2 * amp_2 * w^2 * sum(amp_2 * acos_2^2 + diffy_2 * asin_2)
	  h_phase2_theta2 <- - amp_2 * w * sum(diffy_2 * acos_2)
	  
	  h_theta2_theta2 <- n2 / 2 / theta_2^2
	  
	  hmatrix <- matrix(0,nrow=7,ncol = 7)
	  hmatrix[1,1] <- h_amp1_amp1
	  hmatrix[1,2] <- hmatrix[2,1] <- h_amp1_phase1
	  hmatrix[1,3] <- hmatrix[3,1] <- h_amp1_offsetc
	  hmatrix[1,4] <- hmatrix[4,1] <- h_amp1_theta1
	  
	  hmatrix[2,2] <- h_phase1_phase1
	  hmatrix[2,3] <- hmatrix[3,2] <- h_phase1_offsetc
	  hmatrix[2,4] <- hmatrix[4,2] <- h_phase1_theta1
	  
	  hmatrix[3,3] <- h_offsetc_offsetc
	  hmatrix[3,4] <- hmatrix[4,3] <- h_offsetc_theta1
	  hmatrix[3,5] <- hmatrix[5,3] <- h_offsetc_amp2
	  hmatrix[3,6] <- hmatrix[6,3] <- h_offsetc_phase2
	  hmatrix[3,7] <- hmatrix[7,3] <- h_offsetc_theta2
	  
	  hmatrix[4,4] <- h_theta1_theta1

	  hmatrix[5,5] <- h_amp2_amp2
	  hmatrix[5,6] <- hmatrix[6,5] <- h_amp2_phase2
	  hmatrix[5,7] <- hmatrix[7,5] <- h_amp2_theta2
	  
	  hmatrix[6,6] <- h_phase2_phase2
	  hmatrix[6,7] <- hmatrix[7,6] <- h_phase2_theta2
	  
	  hmatrix[7,7] <- h_theta2_theta2
	  hmatrix
	}
	#h(beta0)
	if(F){
	  betaStar <- c(3.039525, 5.356788, 3.017614, 1.059177, 15.07628, 1.787642, 1.020192)
	  f(betaStar)
	  g(betaStar)
	  h(betaStar)
	  
	  f(beta0)
	  g(beta0)
	  h(beta0)
	  optimx(beta0, fn=f)
	  optimx(betaStar, fn=f)
	  
	  optimx(beta0, fn=f, gr = g, 
	         control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	  optimx(betaStar, fn=f, gr = g, 
	         control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	  optimx(beta0, fn=f, gr = g, hess=h, 
	         lower=c(-Inf, -Inf, -Inf, 0, -Inf, -Inf, 0), 
	         upper=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf), 
	         method = "L-BFGS-B", 
	         control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	  
	}
	anopt <- optimx(beta0, fn=f, gr = g, hess=h, 
	       lower=c(-Inf, -Inf, -Inf, 0, -Inf, -Inf, 0), 
	       upper=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf), 
	       method = "L-BFGS-B", 
	       control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,
					 starttests=FALSE, save.failures=TRUE, trace=0,
					 kkt=FALSE))

	res <- as.numeric(anopt[1:length(beta0)])
	if(any(is.na(res))){
	  bnopt <- optimx(beta0, fn=f, method=c("BFGS"))
	  #bnopt <- optimx(beta0, fn=f, method=c("nlm"))
	  res <- as.numeric(bnopt[1:length(beta0)])
	}
	res
	
}



