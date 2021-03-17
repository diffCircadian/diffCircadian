##' Likelihood-based tests for circadian pattern detection. 
##'
##' Test the significance of circadian curve fitting using likelihood-based tests.
##' @title Likelihood-based Tests for Detecting Circadian Pattern.
##' @param tt Time vector
##' @param yy Expression vector
##' @param period Period of the since curve. Default is 24.
##' @param method Testing methods can be "Wald" or "LR". Default is "LR".
##' @param FN Type of Test, finite sample if TRUE or large sample if FALSE. Default is TRUE. 
##' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue. 
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
##' \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{sigma02}{Variance estimate under the null (intercept only).}
##' \item{sigmaA2}{Variance estimate under the alternative (since curve fitting).}
##' \item{l0}{Log likelihood under the null (intercept only).}
##' \item{l1}{Log likelihood under the alternative (since curve fitting).}
##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##' @author Zhiguang Huo, Haocheng Ding
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' LR_rhythmicity(tt, yy, period=24, method="LR", FN=TRUE)
LR_rhythmicity <- function(tt,yy,period=24,method="LR",FN=TRUE){
  if(method=="Wald"){
    WaldTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else if(method=="LR"){
    LRTest(tt=tt, yy=yy, period = period, type=FN)
  }
  else(("Please check your input! Method only supports 'Wald','LR','F' or 'Permutation'."))
}