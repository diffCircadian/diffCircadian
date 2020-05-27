library(diffCircadian)

set.seed(32608)
B <- 10000
n <- 200

pvalues <- rep(NA, B)
pvalues_f <- rep(NA, B)
for(b in 1:B){
	if(b %% 100 == 0)	print(b)
	set.seed(b)
	tt <- runif(n,0,24) 
	Amp <- 0
	Phase <- 0
	Offset <- 3
	yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
	pvalues[b] <- WaldTest(tt, yy)$pvalue
	pvalues_f[b] <- WaldTest_finiteN(tt, yy)$pvalue
}

mean(pvalues<0.05)
mean(pvalues_f<0.05)

hist(pvalues)



set.seed(32608)
B <- 10000
n <- 10

pvalues <- rep(NA, B)
pvalue_L <- rep(NA, B)
for(b in 1:B){
	if(b %% 100 == 0)	print(b)
	set.seed(b)
	tt <- runif(n,0,24) 
	Amp <- 0
	Phase <- 0
	Offset <- 3
	yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)
	pvalues[b] <- LRTest(tt, yy)$pvalue
	pvalue_L[b] <- LRTest_finiteN(tt, yy)$pvalue
	
}

mean(pvalues<0.05)
mean(pvalue_L<0.05)

hist(pvalues)
