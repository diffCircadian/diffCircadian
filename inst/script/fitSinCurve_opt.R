library(optimx)
library(diffCircadian)

set.seed(32608)
n <- 50
tt <- runif(n,0,24) 
Amp <- 2
Phase <- 6
Offset <- 3
yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
fitSinCurve_opt(tt, yy)


fitSinCurve_opt <- function(tt, yy, period = 24, parStart = c(amp=3,phase=0, offset=0, theta = 1)){

	## amp: vec[1]; phase: vec[2]; offset: vec[3]; theta: vec[4]. 
	n <- length(tt)

	f <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]
		theta <- vec[4]
		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		-n/2 * log(theta) + theta/2 * sum(diffy^2)
	}
	#f(parStart)
	#optimx(parStart, fn=f)
	#unlist(fitSinCurve(tt, yy))

	g <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]
		theta <- vec[4]
		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		acos <- cos(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		g_amp <- -theta * sum(diffy*asin)
		g_phase <- -theta * amp * w * sum(diffy*acos)
		g_offset <- -theta * sum(diffy)
		g_theta <- -n/2/theta + 1/2*sum(diffy^2)
			
		c(g_amp, g_phase, g_offset, g_theta)	
	}
	#g(parStart)
	#optimx(parStart, fn=f, gr = g)
	#unlist(fitSinCurve(tt, yy))

	h <- function(vec){
		amp <- vec[1]
		phase <- vec[2]
		offset <- vec[3]
		theta <- vec[4]
		w <- 2*pi/period
		
		asin <- sin(w * (tt + phase))
		acos <- cos(w * (tt + phase))
		yhat <- amp * asin + offset
		diffy <- yy-yhat
		
		h_amp_amp <- theta * sum(asin^2)
		h_amp_phase <- theta * w * sum((amp * asin - diffy) * acos)
		h_amp_offset <- theta * sum(asin)
		h_amp_theta <- - sum(diffy * asin)
		
		h_phase_phase <- theta * amp * w^2 * sum(amp * acos^2 + diffy * asin)
		h_phase_offset <- theta * amp * w * sum(acos)
		h_phase_theta <- - amp * w * sum(diffy * acos)
			
		h_offset_offset <- theta * n
		h_offset_theta <- - sum(diffy)

		h_theta_theta <- n / 2 / theta^2

		matrix(c(h_amp_amp, h_amp_phase, h_amp_offset, h_amp_theta,
			h_amp_phase, h_phase_phase, h_phase_offset, h_phase_theta,
			h_amp_offset, h_phase_offset, h_offset_offset, h_offset_theta,
			h_amp_theta, h_phase_theta, h_offset_theta, h_theta_theta
			)	
			,ncol=4)
	}
	#h(parStart)
	#optimx(parStart, fn=f, gr = g, hess=h)
	#optimx(parStart, fn=f, gr = g, hess=h, lower=c(-Inf, -Inf, -Inf, 0), upper=c(Inf,Inf,Inf,Inf), method = "L-BFGS-B", control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	#optimx(parStart, fn=f)
	#unlist(fitSinCurve(tt, yy))

	anopt <- optimx(parStart, fn=f, gr = g, hess=h, lower=c(-Inf, -Inf, -Inf, 0), upper=c(Inf,Inf,Inf,Inf), method = "L-BFGS-B", control=list(follow.on = TRUE,usenumDeriv=FALSE,kkttol=10,starttests=FALSE, save.failures=TRUE, trace=0))
	anopt
}



