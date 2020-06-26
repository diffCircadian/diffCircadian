if(F){
	library(devtools)
	install_github("diffCircadian/diffCircadian", auth_token = "9ca607bbfcca3128b7a87d30b175f33383004d90") 	
}

library(diffCircadian)
library(foreach)
library(doParallel)

WD <- "/ufrc/zhuo/zhuo/research/Students/Haocheng/circadianAnalysis/result/simulation/raw"

B <- 10
KK <- 1000
ns <- c(10,20,50,100,200)
numCores <- 32

n <- ns[1]

for(n in ns){
	folder_N <- paste0("N",n)
	thisFolder_N <- file.path(WD, folder_N)
	dir.create(thisFolder_N, re=TRUE)	
	setwd(thisFolder_N)

	b <- 1
	for(b in 1:B){
	#result <- foreach(b = 1:B) %dopar% {
		#
		afile <- paste0("astat_", b, ".csv")	
		bfile <- file.path(thisFolder_N, afile)
		
		if(file.exists(bfile)){
			cat(afile, "exists", "\n")
		} else {

			testStat <- NULL
			
			cl<-makeCluster(numCores)
			registerDoParallel(cl)

			testStat <- foreach(bb = 1:KK, .combine=rbind) %dopar% {
			#for(bb in 1:KK){
				library(diffCircadian)

				bbKK <- bb + KK * (b - 1)
				if(bb%%100==0) print(bb)
				
				set.seed(bbKK)
				tt1 <- runif(n,0,24) 
				Amp1 <- 3
				Phase1 <- 6
				Offset1 <- 3
				yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
				tt2 <- runif(n,0,24) 
				Amp2 <- 3
				Phase2 <- 6
				Offset2 <- 3
				yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)

				atest_LRTest_diff_amp_reg <- LRTest_diff_amp(tt1, yy1, tt2, yy2)				
				atest_LRTest_diff_phase_reg <- LRTest_diff_phase(tt1, yy1, tt2, yy2)	
				atest_LRTest_diff_offset_reg <- LRTest_diff_offset(tt1, yy1, tt2, yy2)	
				atest_LRTest_diff_sigma2_reg <- LRTest_diff_sigma2(tt1, yy1, tt2, yy2)	

				atest_WaldTest_diff_amp_reg <- WaldTest_diff_amp(tt1, yy1, tt2, yy2)	
				atest_WaldTest_diff_phase_reg <- WaldTest_diff_phase(tt1, yy1, tt2, yy2)	
				atest_WaldTest_diff_offset_reg <- WaldTest_diff_offset(tt1, yy1, tt2, yy2)	
				atest_WaldTest_diff_sigma2_reg <- WaldTest_diff_sigma2(tt1, yy1, tt2, yy2)	
				
				atibble <- data.frame(N=n*2,seed = bbKK, Method = c(rep("LR", 4), rep("Wald", 4)), 
						Type = rep(c("amp", "phase", "offset", "sigma2"), 2),
						stat = c(atest_LRTest_diff_amp_reg$stat, 
							atest_LRTest_diff_phase_reg$stat, 
							atest_LRTest_diff_offset_reg$stat, 
							atest_LRTest_diff_sigma2_reg$stat, 
							atest_WaldTest_diff_amp_reg$stat, 
							atest_WaldTest_diff_phase_reg$stat, 
							atest_WaldTest_diff_offset_reg$stat,
							atest_WaldTest_diff_sigma2_reg$stat
							)
				)
				atibble
			}
			stopCluster(cl)				
			write.csv(testStat, bfile, row.names = FALSE)
		} ## end of ifelse
	} ## end of loop for b
} ## end of loop for n


