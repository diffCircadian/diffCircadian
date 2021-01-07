##' Optimization for common sigma square exist
##'
##' Fit a joint sin curve for two conditions assuming common sigmas
##' @title Fit the data based on joint sin curve assuming common sigma square
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param parStart Initial value for optimzation purpose. This is a 2 element vector. Only sigma2_1 and sigma2_2 are needed.
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{sigma2_1}{Variance estimate of the 1st data}
##' \item{sigma2_2}{Variance estimate of the 2nd data}
##' \item{sigma2_C}{Variance estimate pooling all data together}
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
##' opt_commonSigma(tt1, yy1, tt2, yy2)

opt_commonSigma <- function(tt1, yy1, tt2, yy2, period = 24, parStart = NULL){
	
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	if(is.null(parStart)){
		fit1 <- fitSinCurve(tt1, yy1, period = 24)
		fit2 <- fitSinCurve(tt2, yy2, period = 24)
	  sigma2_1 <- 1/n1 * fit1$rss
	  sigma2_2 <- 1/n2 * fit2$rss		
	} else {
		sigma2_1 <- parStart[1]
		sigma2_2 <- parStart[2]
	}
		
  sigma2_C <- 1/(n1 + n2) * (sigma2_1 * n1 + sigma2_2 * n2)

	res <- list(sigma2_1 = sigma2_1, sigma2_2 = sigma2_2, sigma2_C = sigma2_C)
	res
}



