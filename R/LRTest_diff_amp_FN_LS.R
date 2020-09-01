##' Finite sample/Large sample Likelihood ratio test for differential amplitude.
##'
##' Test differential amplitude of circadian curve fitting using likelihood ratio test
##' @title Likelihood ratio test for detecting differential amplitudes. 
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is finite sample. 
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp_1}{amplitude estimate of the 1st data}
##' \item{amp_2}{amplitude estimate of the 2nd data}
##' \item{amp_c}{amplitude estimate pooling all data together}
##' \item{l0}{log likelihood under the null (same variance between the two groups)}
##' \item{l1}{log likelihood under the alternative (different variance between the two groups)}
##' \item{df}{degree of freedom for the LR test}
##' \item{stat}{the LR statistics}
##' \item{pvalue}{the p-value from the LR test}
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
##' LRTest_diff_amp(tt1, yy1, tt2, yy2)


LRTest_diff_amp<- function(tt1, yy1, tt2, yy2, period = 24,type="FN"){
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
  
  beta0 <- c((par1$amp + par2$amp)/2, 
             par1$phase, par1$offset, theta1, 
             par2$phase, par2$offset, theta2
  )
  
  
  this_opt_commonAmp <- opt_commonAmp(tt1, yy1, tt2, yy2, period=period, parStart=beta0)
  
  ## fit H0, individual curves
  ## extract parameters	  
  amp_c <- this_opt_commonAmp[1]
  phase_1 <- this_opt_commonAmp[2]
  offset_1 <- this_opt_commonAmp[3]
  theta_1 <- this_opt_commonAmp[4]
  phase_2 <- this_opt_commonAmp[5] 
  offset_2 <- this_opt_commonAmp[6]
  theta_2 <- this_opt_commonAmp[7] 
  
  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))
  
  yhat_1 <- amp_c * asin_1 + offset_1
  yhat_2 <- amp_c * asin_2 + offset_2
  
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
  
  
  res <- list(amp_1=par1$amp, amp_2=par2$amp, amp_c=amp_c, 
              l0=l0, 
              la=la, 
              df = dfdiff, 
              stat=-2*(l0-la), 
              pvalue=pvalue)
  return(res)
}

