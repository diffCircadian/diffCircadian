library(optimx)
library(diffCircadian)

set.seed(32608)
n <- 50
tt <- runif(n,0,24) 
Amp <- 2
Phase <- 6
Offset <- 3
yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
fitSinCurve_opt3(tt, yy)


fitSinCurve_opt3 <- function(tt, yy, period = 24, parStart = c(amp=3,phase=0, offset=0)){
	## amp: vec[1]; phase: vec[2]; offset: vec[3]; theta: vec[4]. 
	n <- length(tt)

	f <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]
		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		sum(diffy^2)
	}
	#f(parStart)
	#optimx(parStart, fn=f)
	#unlist(fitSinCurve(tt, yy))

  #parStart <- c(-1.915421, -5.463607, 3.028328)
	g <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]
	
		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		acos <- cos(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		g_amp <- - sum(diffy*asin)
		g_phase <- - amp * w * sum(diffy*acos)
		g_offset <- - sum(diffy)

					
		c(g_amp, g_phase, g_offset)	
	}
	#g(parStart)
	#optimx(parStart, fn=f, gr = g)
	#optimx(parStart, fn=f, gr = g, control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	
	
	#unlist(fitSinCurve(tt, yy))

	h <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]

		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		acos <- cos(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		h_amp_amp <- sum(asin^2)
		h_amp_phase <- w * sum((amp * asin - diffy) * acos)
		h_amp_offset <- sum(asin)

		h_phase_phase <- amp * w^2 * sum(amp * acos^2 + diffy * asin)
		h_phase_offset <- amp * w * sum(acos)

					
		h_offset_offset <- n
		
		matrix(c(h_amp_amp, h_amp_phase, h_amp_offset, 
			h_amp_phase, h_phase_phase, h_phase_offset, 
			h_amp_offset, h_phase_offset, h_offset_offset 
			)	
			,ncol=3)
	}
	anopt <- optimx(parStart, fn=f, gr = g, hess=h, lower=c(-Inf, -Inf, -Inf), upper=c(Inf,Inf,Inf), method = "L-BFGS-B", , control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	anopt
}

