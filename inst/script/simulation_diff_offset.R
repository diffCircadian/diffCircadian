library(diffCircadian)

set.seed(32608)
B <- 1000
pvalues <- rep(NA, B)
for(b in 1:B){
	if(b %% 100 == 0)	print(b)
	n <- 50
	set.seed(b)
	tt1 <- runif(n,0,24) 
	Amp1 <- 3
	Phase1 <- 6
	Offset1 <- 3
	yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
	tt2 <- runif(n,0,24) 
	Amp2 <- 2
	Phase2 <- 15
	Offset2 <- 3
	yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
	atest <- LRTest_diff_offset(tt1, yy1, tt2, yy2)	
	pvalues[b] <- atest$pvalue
}

hist(pvalues)


