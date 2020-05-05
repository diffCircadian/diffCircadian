library(diffCircadian)

set.seed(32608)

n <- 50

tt1 <- runif(n,0,24) 
Amp1 <- 3
Phase1 <- 6
Offset1 <- 3
yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
tt2 <- runif(n,0,24) 
Amp2 <- 7
Phase2 <- 15
Offset2 <- 2
yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
atest <- LRTest_diff_amp(tt1, yy1, tt2, yy2)	
atest

