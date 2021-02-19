##' Likelihood-based tests for differential circadian pattern detection.
##'
##' Test differential amplitude/phase/basal level/rhythmicity of circadian curve fitting using likelihood-based tests.
##' @title Likelihood-based Tests for Detecting Differential Circadian Pattern
##' @param tt1 Time vector of condition 1.
##' @param yy1 Expression vector of condition 1.
##' @param tt2 Time vector of condition 2.
##' @param yy2 Expression vector of condition 2.
##' @param period Period of the since curve. Default is 24.
##' @param method Test used to detect differential circadian pattern. It can be chosen either "LR" or "Wald". Default is LR.
##' @param FN "TRUE" if using finite sample likelihood-based tests and "FALSE" if using general large sample likelihood-based tests. Default is "TRUE". 
<<<<<<< HEAD
##' @param type Test differential circadian pattern in differential "amplitude", "phase", "basal", "fit" or "all". Default is "all".   
=======
##' @param type Test differential circadian pattern in differential "amplitude", "phase", "offset" (basal level), "rhythmicity" or "all". Default is "all".   
>>>>>>> parent of 8655d70 (update LR_diff basal)
##' @return A list, see details below. 
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset.}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset.}
##' \item{amp_1}{Amplitude estimate of the 1st data.}
##' \item{amp_2}{Amplitude estimate of the 2nd data.}
##' \item{amp_c}{Amplitude estimate pooling all data together.}
##' \item{phase_1}{Phase estimate of the 1st data, phase is restricted in (0, period).}
##' \item{phase_2}{Phase estimate of the 2nd data, phase is restricted in (0, period).}
##' \item{phase_c}{Phase estimate pooling all data together, phase is restricted in (0, period).}
##' \item{offset_1}{Basal level estimate of the 1st data.}
##' \item{offset_2}{Basal level estimate of the 2nd data.}
##' \item{offset_c}{Basal level estimate pooling all data together.}
##' \item{sigma2_1}{Variance estimate of the 1st data.}
##' \item{sigma2_2}{Variance estimate of the 2nd data.}
##' \item{sigma2_C}{Variance estimate pooling all data together.}

##' \item{l0}{Log likelihood under the null (same variance between the two groups).}
##' \item{l1}{Log likelihood under the alternative (different variance between the two groups).}
##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' @author Zhiguang Huo, Haocheng Ding
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
##' LR_diff(tt1, yy1, tt2, yy2)

LR_diff <- function(tt1,yy1,tt2,yy2,period=24,method="LR",FN=TRUE,type="all"){
  # LR 
  if (method=="LR"){
    # LR + all
    if(type=="all"){
      list(LRTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      LRTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      LRTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      LRTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN))
    }
    # LR + amp
    else if(type=="amplitude"){
      LRTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="phase"){
      LRTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="offset"){
      LRTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else if(type=="fit"){
      LRTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else{
    	stop("Please check your input! type = 'all','amplitude','phase','offset' or 'fit' and test = 'LR' or 'Wald'")
    }
  }
  
  
  
  
  
  # Wald 
  else if (method=="Wald"){
    # Wald + all
    if(type=="all"){
      list(WaldTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      WaldTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      WaldTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN),
      WaldTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN))
    }
    # Wald + amp
    else if(type=="amplitude"){
      WaldTest_diff_amp(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    # Wald + phase
    else if(type=="phase"){
      WaldTest_diff_phase(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    # Wald + offset
    else if(type=="basal"){
      WaldTest_diff_offset(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    # Wald + sigma2
    else if(type=="fit"){
      WaldTest_diff_sigma2(tt1=tt1,yy1=yy1,tt2=tt2,yy2=yy2,period=period,FN=FN)
    }
    else("Please check your input! type = 'all','amplitude','phase','offset' or 'rhythmicity' and test = 'LR' or 'Wald'")
  }
  else("Please check your input! type = 'all','amplitude','phase','basal' or 'fit' and test = 'LR' or 'Wald'")
  
}