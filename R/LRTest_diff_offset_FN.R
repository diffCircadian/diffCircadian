##' Finite sample Likelihood ratio test for differential offset.
##'
##' Test differential offset of circadian curve fitting using likelihood ratio test
##' @title LRTest_diff_offset_FN
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{offset_1}{offset estimate of the 1st data}
##' \item{offset_2}{offset estimate of the 2nd data}
##' \item{offset_c}{offset estimate pooling all data together}
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
##' LRTest_diff_offset_FN(tt1, yy1, tt2, yy2)


LRTest_diff_offset_FN <- function(tt1, yy1, tt2, yy2, period = 24){
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
             par1$phase, 
             (par1$offset + par2$offset)/2,
             theta1, 
             par2$amp, 
             par2$phase, 
             theta2
  )
  
  #opt_commonAmp(tt1, yy1, tt2, yy2, period=period)
  #this_opt_commonAmp <- unlist(opt_commonAmp(tt1, yy1, tt2, yy2, period=period, parStart=beta0))
  this_opt_commonOffset <- opt_commonOffset(tt1, yy1, tt2, yy2, period=period, parStart=beta0)
  
  ## fit H0, individual curves
  ## extract parameters	  
  amp_1 <- this_opt_commonOffset[1]
  phase_1 <- this_opt_commonOffset[2]
  offset_c <- this_opt_commonOffset[3]
  theta_1 <- this_opt_commonOffset[4]
  amp_2 <- this_opt_commonOffset[5]
  phase_2 <- this_opt_commonOffset[6] 
  theta_2 <- this_opt_commonOffset[7] 
  
  asin_1 <- sin(w * (tt1 + phase_1))
  asin_2 <- sin(w * (tt2 + phase_2))
  
  yhat_1 <- amp_1 * asin_1 + offset_c
  yhat_2 <- amp_2 * asin_2 + offset_c
  
  diffy_1 <- yy1-yhat_1
  diffy_2 <- yy2-yhat_2
  
  l1a <- 1/2 * n1 * log(theta_1)
  l1b <- - 1/2 * sum(diffy_1^2) * theta_1
  
  l2a <- 1/2 * n2 * log(theta_2)
  l2b <- - 1/2 * sum(diffy_2^2) * theta_2
  
  l0 <- unlist(l1a + l1b + l2a + l2b)
  
  dfdiff <- 1
  #pvalue <- pchisq(-2*(l0-la),dfdiff,lower.tail = F)
  LR_stat <- -2*(l0-la)
  r <- 1
  k <- 8
  n <- n1+n2
  Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
  pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  
  
  res <- list(offset_1=par1$offset, offset_2=par2$offset, offset_c=offset_c, 
              l0=l0, 
              la=la, 
              df = dfdiff, 
              stat=-2*(l0-la), 
              pvalue=pvalue)
  return(res)
}

