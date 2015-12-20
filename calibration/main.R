rm(list=ls())


load("spxOptionMetrics.rData")
#head(spxData050915)


library(stinepack)


source("BlackScholes.R")
source("plotSVI.R")
source("sviJW.R")
source("plotIvols.R")
source("sviarbitrage.R")
source("optionMetricsToIvols.R")
source("sviRoots.R")
source("sviFit0.R")
source("sviSqrtFit.R")
source("computeImpliedVols.R")
source("sviFitQuarticRoots.R")
source("sviVolSurface.R")
source("svi.R")
source('hybrid_bss.R')
source('bssFit.R')
source('bss.R')
set.seed(9081)
spxData050915$strike_price <- spxData050915$strike_price/1000 
spxIvols050915  <- generateOptionMetricsIvols(spxData050915)





#sviGuess <- list(a = mean(midVar,na.rm=T), b = 0.1, sig = 0.1, rho = - 0.7, m = 0); 

spxOptData <- spxIvols050915
head(spxOptData)

paths<-1e3
n<-400

# Fit each slice individually
bssfit0 <- bssFit(spxOptData, paths, n, kappa=1)
svifit0<- sviFit(spxOptData)
#params <- list(S0=1, xi=0.235^2, eta=1.9, alpha=-0.43, rho=-0.9)
# bssMatrix<- data.frame(params)
# for (i in 2:8){
#   bssMatrix<- rbind(bssMatrix, params)
# }   
# bssMatrix

# 
plotIvols_BSS(ivolData = spxOptData, paths =1e3 , n =400 ,bssMatrix=bssfit0,
              kappa = 1,sviMatrix = svifit0)


