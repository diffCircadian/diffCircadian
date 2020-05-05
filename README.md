# diffCircadian
Omics data circadian and differential circadian analysis

## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("diffCircadian/diffCircadian") 
```

## Citation

* To be updated


## Full tutorial

* To be updated

## Short tutorial for circadian pattern detection

```{R}
library(diffCircadian)

set.seed(32608)
n <- 50
tt <- runif(n,0,24) 
Amp <- 2
Phase <- 6
Offset <- 3
yy <- Amp * sin(2*pi/24 * (tt + Phase)) + Offset + rnorm(n,0,1)

WaldTest(tt, yy) ## Wald test
LRTest(tt, yy) ## LR test
FTest(tt, yy)  ## F test
```

## Short tutorial for differential circadian analysis

```{R}
set.seed(32608)
n <- 50
tt1 <- runif(n,0,24) 
Amp1 <- 2
Phase1 <- 6
Offset1 <- 3
yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
tt2 <- runif(n,0,24) 
Amp2 <- 3
Phase2 <- 5
Offset2 <- 2
yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
## Differential pattern fitting
LRTest_diff_sigma2(tt1, yy1, tt2, yy2)

## Differential amplitute
LRTest_diff_amp(tt1, yy1, tt2, yy2)

## Differential phase 
LRTest_diff_phase(tt1, yy1, tt2, yy2)

## Differential offset 
LRTest_diff_offset(tt1, yy1, tt2, yy2)


```
