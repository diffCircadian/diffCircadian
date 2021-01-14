##' Fisher information matrix when two conditions exist
##'
##' Obtain the Fisher information matrix when two conditions exist
##' @title Fisher information matrix when two conditions exist
##' @param beta parameter vector of 8 with the following order: amp_1, phase_1, offset_1, theta_1, amp_2, phase_2, offset_2, theta_2
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @return The Fisher information matrix, this is a 8*8 matrix, with the same order as the input beta parameter.
##' @author Caleb
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
##' beta <- c(Amp1,Phase1,Offset1,1,Amp2,Phase2,Offset2,2)
##' fisherInformation2(beta, tt1, yy1, tt2, yy2)

fisherInformation2 <- function(beta, tt1, yy1, tt2, yy2, period = 24){
	
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))

	stopifnot(length(beta) == 8)
	
	w <- 2*pi/period
	
  amp_1 <- beta[1]
  phase_1 <- beta[2]
  offset_1 <- beta[3]
  theta_1 <- beta[4]
  amp_2 <- beta[5]
  phase_2 <- beta[6] 
  offset_2 <- beta[7]
  theta_2 <- beta[8] 

  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))
  
  acos_1 <- cos(w * (tt1 + phase_1))
  acos_2 <- cos(w * (tt2 + phase_2))
  
  yhat_1 <- amp_1 * asin_1 + offset_1
  yhat_2 <- amp_2 * asin_2 + offset_2
  
  diffy_1 <- yy1-yhat_1
  diffy_2 <- yy2-yhat_2

  h_amp1_amp1 <- theta_1 * sum(asin_1^2) 
  h_amp1_phase1 <- theta_1 * w * sum((amp_1 * asin_1 - diffy_1) * acos_1)
  h_amp1_offset1 <- theta_1 * sum(asin_1)
  h_amp1_theta1 <- - sum(diffy_1 * asin_1)
  
  h_phase1_phase1 <- theta_1 * amp_1 * w^2 * sum(amp_1 * acos_1^2 + diffy_1 * asin_1)
  h_phase1_offset1 <- theta_1 * amp_1 * w * sum(acos_1)
  h_phase1_theta1 <- - amp_1 * w * sum(diffy_1 * acos_1)
  
  h_offset1_offset1 <- theta_1 * n1
  h_offset1_theta1 <- - sum(diffy_1)
  
  h_theta1_theta1 <- n1 / 2 / theta_1^2
  
  h_amp2_amp2 <- theta_2 * sum(asin_2^2) 
  h_amp2_phase2 <- theta_2 * w * sum((amp_2 * asin_2 - diffy_2) * acos_2)
  h_amp2_offset2 <- theta_2 * sum(asin_2)
  h_amp2_theta2 <- - sum(diffy_2 * asin_2)
  
  h_phase2_phase2 <- theta_2 * amp_2 * w^2 * sum(amp_2 * acos_2^2 + diffy_2 * asin_2)
  h_phase2_offset2 <- theta_2 * amp_2 * w * sum(acos_2)
  h_phase2_theta2 <- - amp_2 * w * sum(diffy_2 * acos_2)
  
  h_offset2_offset2 <- theta_2 * n2
  h_offset2_theta2 <- - sum(diffy_2)
  
  h_theta2_theta2 <- n2 / 2 / theta_2^2
	  	  	  
  hmatrix <- matrix(0,nrow=8,ncol = 8)
	
  hmatrix[1,1] <- h_amp1_amp1
  hmatrix[1,2] <- hmatrix[2,1] <- h_amp1_phase1
  hmatrix[1,3] <- hmatrix[3,1] <- h_amp1_offset1
  hmatrix[1,4] <- hmatrix[4,1] <- h_amp1_theta1  
  hmatrix[2,2] <- h_phase1_phase1
  hmatrix[2,3] <- hmatrix[3,2] <- h_phase1_offset1
  hmatrix[2,4] <- hmatrix[4,2] <- h_phase1_theta1 
  hmatrix[3,3] <- h_offset1_offset1
  hmatrix[3,4] <- hmatrix[4,3] <- h_offset1_theta1
  hmatrix[4,4] <- h_theta1_theta1

  hmatrix[5,5] <- h_amp2_amp2
  hmatrix[5,6] <- hmatrix[6,5] <- h_amp2_phase2
  hmatrix[5,7] <- hmatrix[7,5] <- h_amp2_offset2
  hmatrix[5,8] <- hmatrix[8,5] <- h_amp2_theta2  
  hmatrix[6,6] <- h_phase2_phase2
  hmatrix[6,7] <- hmatrix[7,6] <- h_phase2_offset2
  hmatrix[6,8] <- hmatrix[8,6] <- h_phase2_theta2 
  hmatrix[7,7] <- h_offset2_offset2
  hmatrix[7,8] <- hmatrix[8,7] <- h_offset2_theta2
  hmatrix[8,8] <- h_theta2_theta2

  hmatrix
}

