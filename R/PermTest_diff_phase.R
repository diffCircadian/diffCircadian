##' Permutation test for differential phase
##'
##' Test differential phase of circadian curve fitting using Permutation test
##' @title PermTest_diff_phase
##' @import foreach
##' @import doParallel
##' @import parallel
##' @param tt1 time vector for miltiple genes of condition 1
##' @param yy1 expression matrix for miltiple genes of condition 1
##' @param tt2 time vector for miltiple genes of condition 2
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
##' tt1 <- seq(2,2*n,2)
##' yy1 <- matrix(0,K,n)
##' tt2 <- seq(2,2*n,2)
##' yy2 <- matrix(0,K,n)
##' for(k in 1:K){
##' set.seed(k)
##' tt1 <- runif(n,0,24)
##' tt2 <- runif(n,0,24)
##' Amp1 <- runif(1,1,2)
##' Amp2 <- runif(1,3,4)
##' Phase1 <- runif(1,0,2)
##' Phase2 <- runif(1,5,7)
##' Offset1 <- runif(1,0,3)
##' Offset2 <- runif(1,5,8)
##' yy1[k,] <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' yy2[k,] <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,3)
##' }
##' PermTest_diff_phase(tt1, yy1, tt2, yy2, B=1000)


PermTest_diff_phase <- function(tt1,yy1,tt2,yy2,B=1000,period=24){
  numCores <- detectCores()
  registerDoParallel(numCores)
  
  K <- nrow(yy1)
  stopifnot(K == nrow(yy2))
  
  delta_phase <- rep(0,K)
  phase1 <- rep(0,K)
  phase2 <- rep(0,K)
  pvalue <- rep(0,K)
  
  
  for (k in 1:K){
    par1 <- fitSinCurve(tt1,yy1[k,])
    par2 <- fitSinCurve(tt2,yy2[k,])
    phase1[k] <- par1$phase
    phase2[k] <- par2$phase
    
    delta_phase[k] <- abs(phase1[k]-phase2[k])
    if(delta_phase[k]>12){
      delta_phase[k] <- 24-delta_phase[k]
    }
  }
  
  p <- foreach (k = 1:K,.combine=c) %dopar% {
    perm_delta_phase <- rep(0,B)
    for(b in 1:B){
      set.seed(b)
      
      n1 <- length(tt1)
      n2 <- length(tt2)
      tt1_ind <- rep(1,n1)
      tt2_ind <- rep(2,n2)
      index <- c(tt1_ind,tt2_ind)
      tt12 <- c(tt1,tt2)
      index <- sample(index)
      ttmat <- cbind(tt12,index)
      #tt12 index==1
      perm_tt1 <- ttmat[which(index==1),1]
      perm_tt2 <- ttmat[which(index==2),1]
      
      perm_phase1 <- fitSinCurve(perm_tt1,yy1[k,])$phase
      perm_phase2 <- fitSinCurve(perm_tt2,yy2[k,])$phase
      
      perm_delta_phase[b] <- abs(perm_phase1-perm_phase2)
      if(perm_delta_phase[b]>12){
        perm_delta_phase[b] <- 24-perm_delta_phase[b]
      }
    }
    perm_delta_phase
  }
  
  for(k in 1:K){
    pvalue[k] <- sum(p>delta_phase[k])/(K*B)
  }
  return(list(pvalue=pvalue,phase1=phase1,phase2=phase2,deltaPhase=delta_phase))
}