##' Fit sin function
##'
##' Fit a sine curve where tt is time, and yy is expression value.
##' @title Fit Data Based on Sine Curve
##' @param tt Time vector.
##' @param yy Expression vector.
##' @param period Period of the sine curve. Default is 24.
##' @param parStart Initial value for optimization purpose.
##' @return A list of amp, phase, offset, peak, A, B, SST, SSE, R2. 
##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period \times tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{amp}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, period).}
##' \item{offset}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{A}{A based on formula 2.}
##' \item{B}{B based on formula 2.}
##' \item{tss}{Total sum of square.}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma square.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##' @author Caleb
##' @import minpack.lm
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24) 
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
##' fitSinCurve(tt, yy)

fitSinCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){
	
	getPred <- function(parS, tt) {			
		parS$amp * sin(2*pi/period * (tt + parS$phase)) + parS$offset
	}

	residFun <- function(p, yy, tt) yy - getPred(p,tt)

	nls.out <- nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)
	
	apar <- nls.out$par
	
	amp0 <- apar$amp
	asign <- sign(amp0)
	## restrict amp > 0
	amp <- amp0 * asign
	
	phase0 <- apar$phase
	#phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period 
	phase <- (phase0 + ifelse(asign==1,0,period/2)) %% period 
	offset <- apar$offset
	
	peak <- (period/2 * sign(amp0) - period/4 - phase) %%period
	if(peak > period/4*3) peak = peak - period
	
	A <- amp0 * cos(2*pi/period * phase0)
	B <- amp0 * sin(2*pi/period * phase0)
	
	rss <- sum(nls.out$fvec^2)
	tss <- sum((yy - mean(yy))^2)
	R2 <- 1 - rss/tss
	
	if(F){
		amp <- apar$amp
		phase <- apar$phase
		offset <- apar$offset		
	}
	
	res <- list(amp=amp, phase=phase, offset=offset, peak=peak, A=A, B=B, tss=tss, rss=rss, R2=R2)
	res
}



