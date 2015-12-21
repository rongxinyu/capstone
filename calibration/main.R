rm(list=ls())
load("spxOptionMetrics.rData")
library(stinepack)
source("BlackScholes.R")
source("plotIvols.R")
source("optionMetricsToIvols.R")
source("sviFit0.R")
source("svi.R")
source('hybrid_bss.R')
source('bssFit.R')
source('bss.R')

set.seed(9081)
spxData050915$strike_price <- spxData050915$strike_price/1000
spxOptData <- generateOptionMetricsIvols(spxData050915)
spxOptData <- spxOptData[233:264, ]

paths<-100
n<-100

svifit0 <- sviFit(spxOptData)
bssfit0 <- bssFit(spxOptData, paths, n, kappa=1)
plotIvols_BSS(ivolData=spxOptData, paths=1e3, n=400, bssMatrix=bssfit0, kappa=1, sviMatrix=svifit0)