##' Finite sample/Large sample Wald test for differential phase.
##'
##' Test differential phase of circadian curve fitting using Wald test
##' @title Wald test for detecting differential phase
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{phase_1}{Phase estimate of the 1st data, phase is restricted in (0, period)}
##' \item{phase_2}{Phase estimate of the 2nd data, phase is restricted in (0, period)}
##' \item{phase_c}{Phase estimate pooling all data together, phase is restricted in (0, period)}
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
##' WaldTest_diff_phase(tt1, yy1, tt2, yy2)


WaldTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	#period <- 24
	w <- 2*pi/period
	
	fit1 <- fitSinCurve(tt1, yy1, period = period)
	fit2 <- fitSinCurve(tt2, yy2, period = period)
	
	A1 <- fit1$amp
	A2 <- fit2$amp
	
	phase1 <- fit1$phase
	phase2 <- fit2$phase

	if(phase2 - phase1 > period/2){
		phase2 <- phase2 - period
	} else if(phase1 - phase2 > period/2){
		phase1 <- phase1 - period
	}
	
	basal1 <- fit1$offset
	basal2 <- fit2$offset
	
  sigma2_1 <- 1/n1 * fit1$rss
  sigma2_2 <- 1/n2 * fit2$rss	
	
	theta1 <- 1/sigma2_1
	theta2 <- 1/sigma2_2
	
	p1 <- c(A1, phase1, basal1, theta1)
	p2 <- c(A2, phase2, basal2, theta2)
	
	x_Ha <- c(p1, p2)
			
	eval_f_list <- function(x) {		
		p1 <- x[1:4]
		p2 <- x[5:8]
			
		A1 <- p1[1]
		phase1 <- p1[2]
		basel1 <- p1[3]
		theta1 <- p1[4]				
		asin1 <- sin(w * (tt1 + phase1) )
		acos1 <- cos(w * (tt1 + phase1) )
		yhat1 <- A1 * asin1 + basel1
		
		A2 <- p2[1]
		phase2 <- p2[2]
		basel2 <- p2[3]
		theta2 <- p2[4]
		asin2 <- sin(w * (tt2 + phase2) )
		acos2 <- cos(w * (tt2 + phase2) )
		yhat2 <- A2 * asin2 + basel2
						
		ll1_a <- log(theta1)/2
		ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
		ll1 <- ll1_a - ll1_b
		
		ll2_a <- log(theta2)/2
		ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
		ll2 <- ll2_a - ll2_b
		
		partial_A1 <- - theta1 * sum((yy1 - yhat1) * asin1)
		partial_phase1 <- - theta1 * A1 * w * sum((yy1 - yhat1) * acos1)
		partial_C1 <- - theta1 * sum(yy1 - yhat1) 
		partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

		partial_A2 <- - theta2 * sum((yy2 - yhat2) * asin2)
		partial_phase2 <- - theta2 * A2 * w * sum((yy2 - yhat2) * acos2)
		partial_C2 <- - theta2 * sum(yy2 - yhat2) 
		partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2
		
		
	    return( list( "objective" = -sum(ll1) - sum(ll2),
	                  "gradient"  = c(partial_A1, partial_phase1, partial_C1, partial_theta1, 
						  				partial_A2, partial_phase2, partial_C2, partial_theta2)	
			 		) 
			 )
	}
				
	# Equality constraints
	eval_g_eq <- function(x)
	{
		p1 <- x[1:4]
		p2 <- x[5:8]
			
		phase1 <- p1[2]
		phase2 <- p2[2]
		
		phase1 - phase2
	}
	
	# Equality constraints
	eval_g_eq_jac <- function(x)
	{				
		c(0, 1, 0, 0,
			0, -1, 0, 0)
	}
	
	
	# Lower and upper bounds
	lb <- c(0,-Inf,-Inf,0, 0, -Inf,-Inf, 0)
	ub <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
	#initial values
	
	## Error in is.nloptr(ret) : 
#  If you want to use equality constraints, then you should use one of these algorithms NLOPT_LD_AUGLAG, NLOPT_LN_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LN_AUGLAG_EQ, NLOPT_GN_ISRES, NLOPT_LD_SLSQP

  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  "local_opts" = local_opts
	opts <- list( "algorithm"= "NLOPT_LD_SLSQP",
	              "xtol_rel"= 1.0e-15,
	              "maxeval"= 160000,
	              "local_opts" = local_opts,
	              "print_level" = 0
				  #"check_derivatives"=TRUE
				  )

	res <- nloptr ( x0 = x_Ha,
	                eval_f = eval_f_list,
	                #eval_grad_f=eval_g,					
	                lb = lb,
	                ub = ub,
	                #eval_g_ineq = eval_g_ineq,
	                eval_g_eq = eval_g_eq,
									eval_jac_g_eq = eval_g_eq_jac,
	                opts = opts)

	#
	#x_Ha
	x_H0 <- res$solution

  beta_ha <- x_Ha
  beta_h0 <- x_H0

  #
  I8 <- fisherInformation2(beta_ha, tt1, yy1, tt2, yy2, period=period)	
  beta_diff <- matrix(beta_ha - beta_h0,ncol=1)	
  #stat <- t(beta_diff) %*% I8 %*% beta_diff
  
  beta_diff2 <- beta_diff[c(2,6)]
  I2 <- solve(solve(I8)[c(2,6),c(2,6)])
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
  
  
	phase_c <- x_H0[2]
	phase_c2 <- x_H0[6]
  
  res <- list(phase_1=phase1, phase_2=phase2, phase_c=phase_c, 
              #df = dfdiff, 
              stat = stat, 
              pvalue = pvalue)
  return(res)
}

