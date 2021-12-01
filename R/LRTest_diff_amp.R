##' Finite sample/Large sample Likelihood ratio test for differential amplitude.
##'
##' Test differential amplitude of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential amplitudes. 
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample. 
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{Amplitude estimate of the 1st data}
##' \item{amp_2}{Amplitude estimate of the 2nd data}
##' \item{amp_c}{Amplitude estimate pooling all data together}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
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
##' LRTest_diff_amp(tt1, yy1, tt2, yy2)


LRTest_diff_amp<- function(tt1, yy1, tt2, yy2, period = 24,FN=TRUE){

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
	
	E1 <- A1 * cos(w * phase1)
	F1 <- A1 * sin(w * phase1)

	E2 <- A2 * cos(w * phase2)
	F2 <- A2 * sin(w * phase2)
	
	basal1 <- fit1$offset
	basal2 <- fit2$offset
	
    sigma2_1 <- 1/n1 * fit1$rss
    sigma2_2 <- 1/n2 * fit2$rss	
	
	theta1 <- 1/sigma2_1
	theta2 <- 1/sigma2_2
	
	p1 <- c(E1, F1, basal1, theta1)
	p2 <- c(E2, F2, basal2, theta2)
	
	x_Ha <- c(p1, p2)
		
	asin1 <- sin(w * tt1)
	acos1 <- cos(w * tt1)
	asin2 <- sin(w * tt2)
	acos2 <- cos(w * tt2)
	
	eval_f_list <- function(x,asin1,acos1,asin2,acos2) {		
		p1 <- x[1:4]
		p2 <- x[5:8]
			
		E1 <- p1[1]
		F1 <- p1[2]
		basel1 <- p1[3]
		theta1 <- p1[4]
		yhat1 <- E1 * asin1 + F1 * acos1 + basel1
		
		E2 <- p2[1]
		F2 <- p2[2]
		basel2 <- p2[3]
		theta2 <- p2[4]
		yhat2 <- E2 * asin2 + F2 * acos2 + basel2
						
		ll1_a <- log(theta1)/2
		ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
		ll1 <- ll1_a - ll1_b
		
		ll2_a <- log(theta2)/2
		ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
		ll2 <- ll2_a - ll2_b
		
		partial_E1 <- - theta1 * sum((yy1 - yhat1) * asin1)
		partial_F1 <- - theta1 * sum((yy1 - yhat1) * acos1)
		partial_C1 <- - theta1 * sum(yy1 - yhat1) 
		partial_theta1 <-  sum((yy1 - yhat1)^2)/2 - n1/2/theta1

		partial_E2 <- - theta2 * sum((yy2 - yhat2) * asin2)
		partial_F2 <- - theta2 * sum((yy2 - yhat2) * acos2)
		partial_C2 <- - theta2 * sum(yy2 - yhat2) 
		partial_theta2 <-  sum((yy2 - yhat2)^2)/2 - n2/2/theta2
		
		
	    return( list( "objective" = -sum(ll1) - sum(ll2),
	                  "gradient"  = c(partial_E1, partial_F1, partial_C1, partial_theta1, 
						  				partial_E2, partial_F2, partial_C2, partial_theta2)	
			 		) 
			 )
	}
				
	# Equality constraints
	eval_g_eq <- function(x,asin1,acos1,asin2,acos2)
	{
		p1 <- x[1:4]
		p2 <- x[5:8]
			
		E1 <- p1[1]
		F1 <- p1[2]
		#basel1 <- p1[3]
		theta1 <- p1[4]
		#yhat1 <- E1 * asin1 + F1 * acos1 + basel1
		
		E2 <- p2[1]
		F2 <- p2[2]
		#basel2 <- p2[3]
		theta2 <- p2[4]
		#yhat2 <- E2 * asin2 + F2 * acos2 + basel2
		
		A2_1 <- (E1^2 + F1^2)
		A2_2 <- (E2^2 + F2^2)
		A2_1 - A2_2
	}
	
	# Equality constraints
	eval_g_eq_jac <- function(x,asin1,acos1,asin2,acos2)
	{
		p1 <- x[1:4]
		p2 <- x[5:8]
			
		E1 <- p1[1]
		F1 <- p1[2]
		#basel1 <- p1[3]
		theta1 <- p1[4]
		#yhat1 <- E1 * asin1 + F1 * acos1 + basel1
		
		E2 <- p2[1]
		F2 <- p2[2]
		#basel2 <- p2[3]
		theta2 <- p2[4]
		#yhat2 <- E2 * asin2 + F2 * acos2 + basel2
		
		A2_1 <- (E1^2 + F1^2)
		A2_2 <- (E2^2 + F2^2)
		A2_1 * theta1 - A2_2 * theta2
		
		c(2 * E1, 2 * F1, 0, 0,
			- 2 * E2, - 2 * F2, 0, 0)
	}
	
	
	# Lower and upper bounds
	lb <- c(-Inf,-Inf,-Inf,0, -Inf, -Inf,-Inf, 0)
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
	                opts = opts,
					asin1=asin1,
					acos1=acos1,
					asin2=asin2,
					acos2=acos2)

	#
	#x_Ha
	x_H0 <- res$solution
		
	l0 <- - eval_f_list(x_H0,asin1,acos1,asin2,acos2)$objective
	la <- - eval_f_list(x_Ha,asin1,acos1,asin2,acos2)$objective
	
	LR_stat <- -2*(l0-la)
	
    dfdiff <- 1
    if(!FN){
      pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
    } else if(FN){      
      r <- 1
      k <- 6
      n <- n1+n2
      Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
      pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
    } else{
    	stop("FN has to be TRUE or FALSE")
    }
  	
	amp_c <- sqrt(x_H0[1]^2 + x_H0[2]^2)
	amp_c2 <- sqrt(x_H0[5]^2 + x_H0[6]^2)
  
    res <- list(amp_1=A1, amp_2=A2, amp_c=amp_c, 
                l0=l0, 
                la=la, 
                #df = dfdiff, 
                stat=-2*(l0-la), 
                pvalue=pvalue)
    return(res)
}

