##' Optimization for common sigmas exist
##'
##' Fit a joint sin curve for two conditions assuming common sigmas
##' @title fit the data based on joint sin curve assuming common sigmas
##' @param tt1 time vector of condition 1
##' @param yy1 expression vector of condition 1
##' @param tt2 time vector of condition 2
##' @param yy2 expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param parStart initial value for optimzation purpose
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1}
##' \item{phase}{phase based on formula 1}
##' \item{offset}{offset based on formula 1 or on formula 2}
##' \item{A}{A based on formula 2}
##' \item{B}{B based on formula 2}
##' \item{SST}{Total sum of square}
##' \item{SSE}{Error sum of square}
##' \item{R2}{Pseudo R2 defined as (SST - SSE)/SST}
##' @author Caleb
##' @import minpack.lm
##' @import optimx
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

opt_commonSigma <- function(tt1, yy1, tt2, yy2, period = 24, parStart = list(amp=3,phase=0, offset=0)){
	
	stopifnot(length(tt1) == length(yy1))
	stopifnot(length(tt2) == length(yy2))
	
	getPred <- function(parS, tt) {			
		parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
	}

	residFun <- function(p, yy, tt) yy - getPred(p,tt)

	nls.out1 <- nls.lm(par=parStart, fn = residFun, yy = yy1, tt = tt1)
	nls.out2 <- nls.lm(par=parStart, fn = residFun, yy = yy2, tt = tt2)
	
	apar1 <- nls.out1$par	
	amp0_1 <- apar1$amp
	asign_1 <- sign(amp0_1)
	## restrict amp > 0
	amp_1 <- amp0_1 * asign_1	
	phase0_1 <- apar1$phase
	#phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period 
	phase_1 <- (phase0_1 + ifelse(asign_1==1,0,period/2)) %% period 
	offset_1 <- apar1$offset	
	peak_1 <- (period/2 * sign(amp0_1) - period/4 - phase_1) %%period
	if(peak_1 > period/4*3) peak_1 = peak_1 - period	
	A_1 <- amp0_1 * cos(2*pi/period * phase0_1)
	B_1 <- amp0_1 * sin(2*pi/period * phase0_1)
			
	apar2 <- nls.out2$par	
	amp0_2 <- apar2$amp
	asign_2 <- sign(amp0_2)
	## restrict amp > 0
	amp_2 <- amp0_2 * asign_2	
	phase0_2 <- apar2$phase
	#phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period 
	phase_2 <- (phase0_2 + ifelse(asign_2==1,0,period/2)) %% period 
	offset_2 <- apar2$offset	
	peak_2 <- (period/2 * sign(amp0_2) - period/4 - phase_2) %%period
	if(peak_2 > period/4*3) peak_2 = peak_2 - period	
	A_2 <- amp0_2 * cos(2*pi/period * phase0_2)
	B_2 <- amp0_2 * sin(2*pi/period * phase0_2)

    asin_1 <- sin(2*pi/period * tt1)
    acos_1 <- cos(2*pi/period * tt1)
    yhat_1 <- A_1 * asin_1 + B_1 * acos_1 + offset_1
    residual_1 <- yy1 - yhat_1

    asin_2 <- sin(2*pi/period * tt2)
    acos_2 <- cos(2*pi/period * tt2)
    yhat_2 <- A_2 * asin_2 + B_2 * acos_2 + offset_2
    residual_2 <- yy2 - yhat_2
  	
	n1 <- length(yy1)
	n2 <- length(yy2)
	
    sigma2_1 <- 1/n1 * sum(residual_1^2) 
    sigma2_2 <- 1/n2 * sum(residual_2^2)
    sigma2_C <- 1/(n1 + n2) * (sum(residual_1^2) + sum(residual_2^2))


	res <- list(amp_1=amp_1, phase_1=phase_1, offset_1=offset_1, peak_1=peak_1, A_1=A_1, B_1=B_1, sigma2_1 = sigma2_1, n1=n1, 
				amp_2=amp_2, phase_2=phase_2, offset_2=offset_2, peak_2=peak_2, A_2=A_2, B_2=B_2, sigma2_2 = sigma2_2, n2=n2, 
				sigma2_C = sigma2_C	
		)
	res
}



