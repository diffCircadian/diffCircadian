##' Finite sample/Large sample Wald test for differential amplitude.
##'
##' Test differential amplitude of circadian curve fitting using Wald test
##' @title Wald test for detecting differential amplitude
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of Wald test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{Amplitude estimate of the 1st data}
##' \item{amp_2}{Amplitude estimate of the 2nd data}
##' \item{amp_c}{Amplitude estimate pooling all data together}
##' \item{df}{Degree of freedom for the Wald test}
##' \item{stat}{Wald statistics}
##' \item{pvalue}{P-value from the Wald test}
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


WaldTest_diff_amp <- function(tt1, yy1, tt2, yy2, period = 24,type="FN"){
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
  #stat <- t(beta_diff) %*% I8 %*% beta_diff
  
  beta_diff2 <- beta_diff[c(1,5)]
  I2 <- solve(solve(I8)[c(1,5),c(1,5)])
  stat <- as.numeric(t(beta_diff2) %*% I2 %*% beta_diff2)
  
  dfdiff <- 1
  
  if(type=="LS"){
    pvalue <- pchisq(stat,dfdiff,lower.tail = F)
  }
  else if(type=="FN"){
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- stat*(n-k)/n/r
    pvalue <- pf(Fstat,df1=r,df2=n-k,lower.tail = F)
  }
  
  
  res <- list(amp_1=par1$amp, amp_2=par2$amp, amp_c=this_opt_commonAmp[1], 
              df = dfdiff, 
              stat = stat, 
              pvalue = pvalue)
  return(res)
}

