rm(list=ls())

library(stinepack)
# 
# download.file(url="http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip",
#               destfile="oxfordRvData.zip")
# unzip(zipfile="oxfordRvData.zip")
# OxfordManRVRaw.data <- read.csv("OxfordManRealizedVolatilityIndices.csv", stringsAsFactors = F)
# save(OxfordManRVRaw.data, file='OxfordManRVRaw.rData')

#load("spxOptionMetrics.rData")
load('OxfordManRVRaw.rData')

source("BlackScholes.R")
source("optionMetricsToIvols.R")
source("sviFit0.R")
source("svi.R")
source('sviarbitrage.R')
source('sviRoots.R')
source('sviSqrtFit.R')
source("computeImpliedVols.R")
source("sviFitQuarticRoots.R")
source("sviVolSurface.R")
source('plotIvols_svibss.R')
source('hybrid_bss.R')
source('bssFit.R')
source('bss.R')
source('bssHistFit.R')

download.file(url="http://mfe.baruch.cuny.edu/wp-content/uploads/2014/09/spxOptionMetrics.rData_.zip", destfile="spxOptionMetrics.rData.zip")
unzip(zipfile="spxOptionMetrics.rData.zip")


load("spxOptionMetrics.rData")


set.seed(9081)
spxData110915$strike_price <- spxData110915$strike_price/1000
spxOptData <- generateOptionMetricsIvols(spxData110915)
# spxOptData <- spxOptData[spxOptData$Texp>0.5 & spxOptData$Texp<0.78,  ]

symbol<- 'SPX2.rk'
endDate<- '2011-09-15'
# rvData<- OxfordManRV.clean(rv.data = OxfordManRVRaw.data)
# rvParam<- OxfordManRV.param(rv.list = rvData, endDate = '2011-09-15', symbol = symbol)


paths<-1e3*5
n<-400


svifitSqrt <- sviSqrtFit(spxOptData)
svifitQR<- sviFitQR(ivolData = spxOptData, sviGuess = svifitSqrt )

# ### Estimate the bss params based on pure fitting optimization  
# #bssfit0 <- bssFit(spxOptData, paths, n, kappa=1)
# ###
# 

t0<- proc.time()

bssMatrix<- bss_histfit(ivolData = spxOptData, OxfordManRVData = OxfordManRVRaw.data, 
            symbol = symbol, EndDate = endDate, paths = paths, n = n)

proc.time()-t0

# 
# bssparam<- list(S0=1, xi= 0.235^2, eta= rvParam$eta, alpha= rvParam$h-.5, rho=-.9)
plotIvols_svibss(ivolData=spxOptData, paths=paths, n=400, 
              bssMatrix=bssMatrix,  kappa=1, sviMatrix=svifitQR)
