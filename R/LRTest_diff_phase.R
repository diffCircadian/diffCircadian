##' Finite sample/Large sample Likelihood ratio test for differential phase.
##'
##' Test differential phase of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential phase.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{phase_1}{Phase estimate of the 1st data, phase is restricted in (0, period)}
##' \item{phase_2}{Phase estimate of the 2nd data, phase is restricted in (0, period)}
##' \item{phase_c}{Phase estimate pooling all data together, phase is restricted in (0, period)}
##' \item{l0}{Log likelihood under the null (same variance between the two groups)}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{Degree of freedom for the LR test}
##' \item{stat}{LR statistics}
##' \item{pvalue}{P-value from the LR test}
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
##' LRTest_diff_phase(tt1, yy1, tt2, yy2)


LRTest_diff_phase <- function(tt1, yy1, tt2, yy2, period = 24,type="FN"){
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
  
  l1_Ha <- 1/2 * n1 * log(theta1) - 1/2 * n1
  l2_Ha <- 1/2 * n2 * log(theta2) - 1/2 * n2
  
  la <- unlist(l1_Ha + l2_Ha)
  
  beta0 <- c(par1$amp, 
             (par1$phase + par2$phase)/2,
             par1$offset,
             theta1, 
             par2$amp, 
             par2$offset, 
             theta2
  )
  
  #opt_commonAmp(tt1, yy1, tt2, yy2, period=period)
  #this_opt_commonAmp <- unlist(opt_commonAmp(tt1, yy1, tt2, yy2, period=period, parStart=beta0))
  this_opt_commonPhase <- opt_commonPhase(tt1, yy1, tt2, yy2, period=period, parStart=beta0)
  
  ## fit H0, individual curves
  ## extract parameters	  
  amp_1 <- this_opt_commonPhase[1]
  phase_c <- this_opt_commonPhase[2]
  offset_1 <- this_opt_commonPhase[3]
  theta_1 <- this_opt_commonPhase[4]
  amp_2 <- this_opt_commonPhase[5]
  offset_2 <- this_opt_commonPhase[6] 
  theta_2 <- this_opt_commonPhase[7] 
  
  asin_1 <- sin(w * (tt1 + phase_c))
  asin_2 <- sin(w * (tt2 + phase_c))
  
  yhat_1 <- amp_1 * asin_1 + offset_1
  yhat_2 <- amp_2 * asin_2 + offset_2
  
  diffy_1 <- yy1-yhat_1
  diffy_2 <- yy2-yhat_2
  
  l1a <- 1/2 * n1 * log(theta_1)
  l1b <- - 1/2 * sum(diffy_1^2) * theta_1
  
  l2a <- 1/2 * n2 * log(theta_2)
  l2b <- - 1/2 * sum(diffy_2^2) * theta_2
  
  l0 <- unlist(l1a + l1b + l2a + l2b)
  
  dfdiff <- 1
  
  if(type=="LS"){
    pvalue <- pchisq(-2*(l0-la),dfdiff,lower.tail = F)
  }
  else if(type=="FN"){
    LR_stat <- -2*(l0-la)
    r <- 1
    k <- 6
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  }
  
  
  res <- list(phase_1=par1$phase, phase_2=par2$phase, phase_c=phase_c, 
              l0=l0, 
              la=la, 
              df = dfdiff, 
              stat=-2*(l0-la), 
              pvalue=pvalue)
  return(res)
}

