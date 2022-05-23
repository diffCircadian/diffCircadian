# diffCircadian
Omics data circadian and differential circadian analysis

## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("diffCircadian/diffCircadian") 
```

## Citation

* Ding, H., Meng, L., Liu, A.C., Gumz, M.L., Bryant, A.J., Mcclung, C.A., Tseng, G.C., Esser, K.A. and Huo, Z., 2021. Likelihood-based Tests for Detecting Circadian Rhythmicity and Differential Circadian Patterns in Transcriptomic Applications. *Briefings in Bioinformatcs* 

* The manuscript can be found here: [https://www.biorxiv.org/content/10.1101/2021.02.23.432538v1](https://academic.oup.com/bib/article-abstract/22/6/bbab224/6297167)


## Full tutorial

http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
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

LR <- LR_rhythmicity(tt, yy) ## LR test
LR

## curve fitting visualization

tt0 <- seq(0,24,0.1) 
yy0 <- LR$amp * sin(2*pi/24 * (tt0 + LR$phase)) + LR$offset
plot(tt,yy,pch=20)
lines(tt0,yy0,col="red",lwd=2)

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
```


```{R}
## Differential pattern fitting
LR_diff(tt1, yy1, tt2, yy2, type="fit")

## Differential amplitute
LR_diff(tt1, yy1, tt2, yy2, type="amplitude")

## Differential phase 
LR_diff(tt1, yy1, tt2, yy2, type="phase")

## Differential offset (basal level)
LR_diff(tt1, yy1, tt2, yy2, type="basal")
```

