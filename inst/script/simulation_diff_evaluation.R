if(F){
	library(devtools)
	install_github("diffCircadian/diffCircadian", auth_token = "9ca607bbfcca3128b7a87d30b175f33383004d90") 	
}

library(tidyverse)

WD <- "/ufrc/zhuo/zhuo/research/Students/Haocheng/circadianAnalysis/result/simulation/raw"

B <- 10
KK <- 1000
ns <- c(10,20,50,100,200)
numCores <- 32

n <- ns[1]

testStat_all <- NULL

for(n in ns){
	folder_N <- paste0("N",n)
	thisFolder_N <- file.path(WD, folder_N)
	setwd(thisFolder_N)

	b <- 1
	
	for(b in 1:B){
		afile <- paste0("astat_", b, ".csv")	
		bfile <- file.path(thisFolder_N, afile)
		
		testStat <- read_csv(bfile)
		testStat_all <- rbind(testStat_all, testStat)
	} ## end of loop for b
}

testStat_all_Wald <- testStat_all %>% 
	filter(Method=="Wald") %>% 
	mutate(Fstat_k8r1 = stat*(N-8)/N/1, 
			Fstat_k7r1 = stat*(N-7)/N/1, 
			Fstat_k6r1 = stat*(N-6)/N/1, 
	)

testStat_all_LR <- testStat_all %>% 
	filter(Method=="LR") %>% 
	mutate(Fstat_k8r1 = (exp(stat/N) - 1) * (N-8) / 1, 
			Fstat_k7r1 = (exp(stat/N) - 1) * (N-7) / 1, 
			Fstat_k6r1 = (exp(stat/N) - 1) * (N-6) / 1 
	)
testStat_all2 <- rbind(testStat_all_Wald, testStat_all_LR)

testStat_all3 <- testStat_all2 %>% 
	mutate(pvalue_reg = pchisq(stat,1,lower.tail=FALSE), 
				pvalue_k8r1 = pf(Fstat_k8r1,df1 = 1, df2 = N-8, lower.tail = F),
				pvalue_k7r1 = pf(Fstat_k7r1,df1 = 1, df2 = N-7, lower.tail = F),
				pvalue_k6r1 = pf(Fstat_k6r1,df1 = 1, df2 = N-6, lower.tail = F)			
	)


#
table_summary <- testStat_all3 %>% 
	group_by(N, Method, Type) %>%
	summarize(ErrorI_reg = mean(pvalue_reg < 0.05), 
						ErrorI_k8r1 = mean(pvalue_k8r1 < 0.05), 
						ErrorI_k7r1 = mean(pvalue_k7r1 < 0.05), 
						ErrorI_k6r1 = mean(pvalue_k6r1 < 0.05)
	)
#
setwd("~")
write.csv(table_summary, "table_summary.csv", row.names=FALSE)

table_summary %>% as.data.frame()

table_summary %>% filter(Type=="sigma2", Method=="LR")
table_summary %>% filter(Type=="sigma2", Method=="Wald")

table_summary %>% filter(Type=="amp", Method=="LR")
table_summary %>% filter(Type=="amp", Method=="Wald")

table_summary %>% filter(Type=="offset", Method=="LR")
table_summary %>% filter(Type=="offset", Method=="Wald")

table_summary %>% filter(Type=="phase", Method=="LR")
table_summary %>% filter(Type=="phase", Method=="Wald")

