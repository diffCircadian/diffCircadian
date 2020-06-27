##' Finite sample/Large sample Likelihood ratio test for differential sigma2
##'
##' Test differential sigma2 of circadian curve fitting using likelihood ratio test
##' @title LRTest_diff_sigma2_FN_LS
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param type Type of likelihood ratio test to use, "FN" or "LS". Default is large sample.
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{sigma2_1}{variance estimate of the 1st data}
##' \item{sigma2_2}{variance estimate of the 2nd data}
##' \item{sigma2_C}{variance estimate pooling all data together}
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
##' LRTest_diff_sigma2_FN_LS(tt1, yy1, tt2, yy2)


LRTest_diff_sigma2_FN_LS <- function(tt1, yy1, tt2, yy2, period = 24,type="LS"){
  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))
  
  this_opt_commonSigma <- opt_commonSigma(tt1, yy1, tt2, yy2, period=period)
  
  ## extract parameters	  
  sigma2_1 <- this_opt_commonSigma$sigma2_1
  sigma2_2 <- this_opt_commonSigma$sigma2_2
  sigma2_C <- this_opt_commonSigma$sigma2_C
  
  ## H0 (sigma2_C)
  l0_1 <- -n1/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_1^2)
  l0_2 <- -n2/2*log(2*pi*sigma2_C) # - 1/(2*sigma2_C)*sum(residual_2^2)
  l0 <- l0_1 + l0_2
  
  ## H1 (sigma2_1, sigma2_2)
  la_1 <- -n1/2*log(2*pi*sigma2_1) # - 1/(2*sigma2_1)*sum(residual_1^2)
  la_2 <- -n2/2*log(2*pi*sigma2_2) # - 1/(2*sigma2_2)*sum(residual_2^2)
  la <- la_1 + la_2
  
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
  
  
  res <- list(sigma2_1=sigma2_1, sigma2_2=sigma2_2, sigma2_c=sigma2_C,
              l0=l0, 
              la=la, 
              df = dfdiff, 
              stat=-2*(l0-la), 
              pvalue=pvalue)
  return(res)
}