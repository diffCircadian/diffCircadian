##' Permutation test for differential R2
##'
##' Test differential R2 of circadian curve fitting using Permutation test
##' @title PermTest_diff_R2
##' @param tt1 time matrix for miltiple genes of condition 1
##' @param yy1 expression matrix for miltiple genes of condition 1
##' @param tt2 time matrix for miltiple genes of condition 2
##' @param yy2 expression matrix for miltiple genes of condition 2
##' @param B number of permutation for each gene.
##' @param period Period of the since curve. Default is 24.
##' @return pvalues, see details below. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{pvalue}{the p-value from the Permutation test}
##' @author Haocheng Ding
##' @export
##' @examples
##' K=100
##' n=10
##' tt1 <- matrix(0,K,n)
##' yy1 <- matrix(0,K,n)
##' tt2 <- matrix(0,K,n)
##' yy2 <- matrix(0,K,n)
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' for(k in 1:K){
##' set.seed(k)
##' tt1[k,] <- runif(n,0,24) 
##' yy1[k,] <- Amp1 * sin(2*pi/24 * (tt1[k,] + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2[k,] <- runif(n,0,24) 
##' yy2[k,] <- Amp2 * sin(2*pi/24 * (tt2[k,] + Phase2)) + Offset2 + rnorm(n,0,5)
##' }
##' PermTest_diff_R2(tt1, yy1, tt2, yy2, B=1000)


library(foreach)
library(doParallel)
PermTest_diff_R2 <- function(tt1,yy1,tt2,yy2,B=1000,period=24){
  numCores <- detectCores()
  registerDoParallel(numCores)
  
  K <- nrow(tt1)
  stopifnot(K == nrow(tt2))
  stopifnot(nrow(yy1) == nrow(yy2))
  
  delta_R2 <- rep(0,K)
  pvalue <- rep(0,K)
  
  
  for (k in 1:K){
    par1 <- fitSinCurve(tt1[k,],yy1[k,])
    par2 <- fitSinCurve(tt2[k,],yy2[k,])
    R2_1 <- 1-par1$rss/par1$tss
    R2_2 <- 1-par2$rss/par2$tss
    delta_R2[k] <- abs(R2_1-R2_2)
  }
  p <- foreach (k = 1:K,.combine=c) %dopar% {
    perm_delta_R2 <- rep(0,B)
    for(b in 1:B){
      set.seed(b)
      
      perm_yy1 <- sample(c(yy1[k,],yy2[k,]),length(yy1[k,]),replace = T)
      perm_yy2 <- sample(c(yy1[k,],yy2[k,]),length(yy2[k,]),replace = T)
      
      perm_par1 <- fitSinCurve(tt1[k,],perm_yy1)
      perm_par2 <- fitSinCurve(tt2[k,],perm_yy2)
      perm_R2_1 <- 1-perm_par1$rss/perm_par1$tss
      perm_R2_2 <- 1-perm_par2$rss/perm_par2$tss
      
      perm_delta_R2[b] <- abs(perm_R2_1-perm_R2_2)
    }
    perm_delta_R2
  }
  
  for(k in 1:K){
    pvalue[k] <- sum(p>delta_R2[k])/(K*B)
  }
  return(pvalue)
}