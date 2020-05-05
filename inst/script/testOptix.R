library(optimx)

## e^(x+y) + x^4 + y^2
f <- function(x){
  exp(x[1] + x[2]) + x[1]^4 + x[2]^2
}
f(c(1,2))

startx <- c(1,2)
f(startx)
optimx(startx, fn=f)

g <- function(x){
  p1 <- exp(x[1] + x[2]) + 4*x[1]^3 
  p2 <- exp(x[1] + x[2]) + 2*x[2] 
  c(p1, p2)
}
startx <- c(1,2)


g <- function(x){
  p1 <- exp(x[1] + x[2]) + 4*x[1]^3 
  p2 <- exp(x[1] + x[2]) + 2*x[2] 
  c(p1, p2)
}

g(startx)
optimx(startx, fn=f, gr = g)

h <- function(x){
  h11 <- exp(x[1] + x[2]) + 12*x[1]^2 
  h12 <- h21 <- exp(x[1] + x[2]) 
  h22 <- exp(x[1] + x[2]) + 2
  matrix(c(h11,h21,h12,h22),ncol=2)
}

h(startx)
optimx(startx, fn=f, gr = g, hess=h, lower=c(0, 1), upper=c(10,10), method = "L-BFGS-B")
optimx(startx, fn=f, gr = g, hess=h, lower=c(-Inf, -Inf), upper=c(10,10), method = "L-BFGS-B")
optimx(startx, fn=f, gr = g, hess=h, lower=1, upper=5, method = "L-BFGS-B")
optimx(startx, fn=f, gr = g, hess=h, lower=-Inf, upper=Inf, method = "BFGS")
