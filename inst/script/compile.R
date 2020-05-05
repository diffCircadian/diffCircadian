library(devtools)

WD_root <- "/Users/zhuo/Desktop/diffCircadian"

setwd(WD_root)
document()
install()

library(diffCircadian)
help(package = "diffCircadian")

?fitSinCurve
?Ftest
?LRtest

set.seed(32608)
n <- 50
tt1 <- runif(n,0,24) 
Amp1 <- 2
Phase1 <- 6
Offset1 <- 3
yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
tt2 <- runif(n,0,24) 
Amp2 <- 3
Phase2 <- 15
Offset2 <- 2
yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
opt_commonAmp(tt1, yy1, tt2, yy2)


set.seed(32608)
n <- 50
tt1 <- runif(n,0,24) 
Amp1 <- 2
Phase1 <- 6
Offset1 <- 3
yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
tt2 <- runif(n,0,24) 
Amp2 <- 3
Phase2 <- 15
Offset2 <- 2
yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
opt_commonSigma(tt1, yy1, tt2, yy2)
LRTest_diff_sigma2(tt1, yy1, tt2, yy2)
