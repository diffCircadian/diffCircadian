##' Optimization for common amplitude exist
##'
##' Fit a joint sin curve for two conditions assuming common amplitude
##' @title fit the data based on joint sin curve assuming common amplitude
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param parStart initial value for optimzation purpose. This has the same order as the output vector.
##' @return A vector of 7 with the following order: amp_c, phase_1, offset_1, theta_1, phase_2, offset_2, theta_2
##' @author Caleb
##' @import optimx
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
##' opt_commonAmp(tt1, yy1, tt2, yy2)

opt_commonAmp <- function(tt1, yy1, tt2, yy2, period = 24, 
                          parStart = NULL){
	
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	## initialization
	if(is.null(parStart)){
	  par1 <- fitSinCurve(tt1,yy1)
	  yhat1 <- par1$amp * sin(2*pi/period * (tt1 + par1$phase)) + par1$offset
	  theta1 <- n1/par1$rss

	  par2 <- fitSinCurve(tt2,yy2)
	  yhat2 <- par2$amp * sin(2*pi/period * (tt2 + par2$phase)) + par2$offset
	  theta2 <- n2/par2$rss
	  
	  beta0 <- c((par1$amp + par2$amp)/2, 
	                par1$phase, par1$offset, theta1, 
	                par2$phase, par2$offset, theta2
	              )
	} else {
		beta0 <- parStart
	}
	
	w <- 2*pi/period
	
	# beta[1]: amp_c
	# beta[2]: phase_1
	# beta[3]: offset_1
	# beta[4]: theta_1
	# beta[5]: phase_2
	# beta[6]: offset_2
	# beta[7]: theta_2
	
	f <- function(beta){
	  amp_c <- beta[1]
	  phase_1 <- beta[2]
	  offset_1 <- beta[3]
	  theta_1 <- beta[4]
	  phase_2 <- beta[5] 
	  offset_2 <- beta[6]
	  theta_2 <- beta[7] 

	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_c * asin_1 + offset_1
	  yhat_2 <- amp_c * asin_2 + offset_2
	  
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
	  amp_c <- beta[1]
	  phase_1 <- beta[2]
	  offset_1 <- beta[3]
	  theta_1 <- beta[4]
	  phase_2 <- beta[5] 
	  offset_2 <- beta[6]
	  theta_2 <- beta[7] 
	  
	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  acos_1 <- cos(w * (tt1 + phase_1))
	  acos_2 <- cos(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_c * asin_1 + offset_1
	  yhat_2 <- amp_c * asin_2 + offset_2
	  
	  diffy_1 <- yy1-yhat_1
	  diffy_2 <- yy2-yhat_2
	  
	  partial_amp_c_1 <- - theta_1 * sum(diffy_1 * asin_1) 
	  partial_amp_c_2 <- - theta_2 * sum(diffy_2 * asin_2) 
	  partial_amp_c <- (partial_amp_c_1 + partial_amp_c_2)
	  
	  partial_phase_1 <- - amp_c * w * theta_1 * sum(diffy_1 * acos_1) 
	  partial_phase_2 <- - amp_c * w * theta_2 * sum(diffy_2 * acos_2) 
	  
	  partial_offset_1 <- - theta_1 * sum(diffy_1) 
	  partial_offset_2 <- - theta_2 * sum(diffy_2) 
	  
	  partial_theta_1 <- - n1 / 2 / theta_1 + sum(diffy_1^2)/ 2 
	  partial_theta_2 <- - n2 / 2 / theta_2 + sum(diffy_2^2)/ 2 
	  
	  c(partial_amp_c, 
	    partial_phase_1, 
	    partial_offset_1, 
	    partial_theta_1, 
	    partial_phase_2, 
	    partial_offset_2, 
	    partial_theta_2)
	}
	
	h <- function(beta){
	  amp_c <- beta[1]
	  phase_1 <- beta[2]
	  offset_1 <- beta[3]
	  theta_1 <- beta[4]
	  phase_2 <- beta[5] 
	  offset_2 <- beta[6]
	  theta_2 <- beta[7] 
	  
	  asin_1 <- sin(w * (tt1 + phase_1))
	  asin_2 <- sin(w * (tt2 + phase_2))
	  
	  acos_1 <- cos(w * (tt1 + phase_1))
	  acos_2 <- cos(w * (tt2 + phase_2))
	  
	  yhat_1 <- amp_c * asin_1 + offset_1
	  yhat_2 <- amp_c * asin_2 + offset_2
	  
	  diffy_1 <- yy1-yhat_1
	  diffy_2 <- yy2-yhat_2
	  
	  h_ampc_ampc <- theta_1 * sum(asin_1^2) + theta_2 * sum(asin_2^2)
	  h_ampc_phase1 <- theta_1 * w * sum((amp_c * asin_1 - diffy_1) * acos_1)
	  h_ampc_offset1 <- theta_1 * sum(asin_1)
	  h_ampc_theta1 <- - sum(diffy_1 * asin_1)
	  h_ampc_phase2 <- theta_2 * w * sum((amp_c * asin_2 - diffy_2) * acos_2)
	  h_ampc_offset2 <- theta_2 * sum(asin_2)
	  h_ampc_theta2 <- - sum(diffy_2 * asin_2)
	  
	  h_phase1_phase1 <- theta_1 * amp_c * w^2 * sum(amp_c * acos_1^2 + diffy_1 * asin_1)
	  h_phase1_offset1 <- theta_1 * amp_c * w * sum(acos_1)
	  h_phase1_theta1 <- - amp_c * w * sum(diffy_1 * acos_1)
	  
	  h_offset1_offset1 <- theta_1 * n1
	  h_offset1_theta1 <- - sum(diffy_1)
	  
	  h_theta1_theta1 <- n1 / 2 / theta_1^2
	  
	  h_phase2_phase2 <- theta_2 * amp_c * w^2 * sum(amp_c * acos_2^2 + diffy_2 * asin_2)
	  h_phase2_offset2 <- theta_2 * amp_c * w * sum(acos_2)
	  h_phase2_theta2 <- - amp_c * w * sum(diffy_2 * acos_2)
	  
	  h_offset2_offset2 <- theta_2 * n2
	  h_offset2_theta2 <- - sum(diffy_2)
	  
	  h_theta2_theta2 <- n2 / 2 / theta_2^2
	  
	  hmatrix <- matrix(0,nrow=7,ncol = 7)
	  hmatrix[1,1] <- h_ampc_ampc
	  hmatrix[1,2] <- hmatrix[2,1] <- h_ampc_phase1
	  hmatrix[1,3] <- hmatrix[3,1] <- h_ampc_offset1
	  hmatrix[1,4] <- hmatrix[4,1] <- h_ampc_theta1
	  hmatrix[1,5] <- hmatrix[5,1] <- h_ampc_phase2
	  hmatrix[1,6] <- hmatrix[6,1] <- h_ampc_offset2
	  hmatrix[1,7] <- hmatrix[7,1] <- h_ampc_theta2
	  
	  hmatrix[2,2] <- h_phase1_phase1
	  hmatrix[2,3] <- hmatrix[3,2] <- h_phase1_offset1
	  hmatrix[2,4] <- hmatrix[4,2] <- h_phase1_theta1
	  
	  hmatrix[3,3] <- h_offset1_offset1
	  hmatrix[3,4] <- hmatrix[4,3] <- h_offset1_theta1
	
	  hmatrix[4,4] <- h_theta1_theta1

	  hmatrix[5,5] <- h_phase2_phase2
	  hmatrix[5,6] <- hmatrix[6,5] <- h_phase2_offset2
	  hmatrix[5,7] <- hmatrix[7,5] <- h_phase2_theta2
	  
	  hmatrix[6,6] <- h_offset2_offset2
	  hmatrix[6,7] <- hmatrix[7,6] <- h_offset2_theta2
	  
	  hmatrix[7,7] <- h_theta2_theta2
	  hmatrix
	}

	anopt <- optimx(beta0, fn=f, gr = g, hess=h, 
	       lower=c(-Inf, -Inf, -Inf, 0, -Inf, -Inf, 0), 
	       upper=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf), 
	       method = "L-BFGS-B", 
	       control=list(follow.on = FALSE,usenumDeriv=FALSE,
					 starttests=FALSE, save.failures=FALSE, trace=0, 
					 dowarn=FALSE, 
					 kkt=FALSE))
	
	res <- as.numeric(anopt[1:length(beta0)])
	if(any(is.na(res))){
	  bnopt <- suppressWarnings(optimx(beta0, fn=f, method=c("BFGS")))
	  #bnopt <- optimx(beta0, fn=f, method=c("nlm"))
	  res <- as.numeric(bnopt[1:length(beta0)])
	}
	res
}

