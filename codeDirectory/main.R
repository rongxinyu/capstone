rm(list=ls())

library(stinepack)
# 
# download.file(url="http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip",
#               destfile="oxfordRvData.zip")
# unzip(zipfile="oxfordRvData.zip")
# OxfordManRVRaw.data <- read.csv("OxfordManRealizedVolatilityIndices.csv", stringsAsFactors = F)
# save(OxfordManRVRaw.data, file='OxfordManRVRaw.rData')

#load("spxOptionMetrics.rData")
# load('OxfordManRVRaw.rData')

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


# download.file(url="http://mfe.baruch.cuny.edu/wp-content/uploads/2014/09/spxOptionMetrics.rData_.zip", destfile="spxOptionMetrics.rData.zip")
# unzip(zipfile="spxOptionMetrics.rData.zip")


load("spxOptionMetrics.rData")


spxData110915$strike_price <- spxData110915$strike_price/1000
spxOptData <- generateOptionMetricsIvols(spxData110915)
# spxOptData <- spxOptData[spxOptData$Texp>0.04 ,  ]

symbol<- 'SPX2.rk'
endDate<- '2011-09-15'
# rvData<- OxfordManRV.clean(rv.data = OxfordManRVRaw.data)
# rvParam<- OxfordManRV.param(rv.list = rvData, endDate = '2011-09-15', symbol = symbol)


paths<-1e3*5
n<-400


# alpha<- rep(-0.43, 14)
# S0<- rep(1,14)
# eta<- rep(1.9,14)
# rho<- rep(-.9, 14)
# xi<- rep(0.235^2,14)
# 
# 
# bssmatrix<- data.frame(S0=S0, alpha= alpha, eta=eta, xi=xi, rho=rho)



bssmatrix<- bssSqFit(spxOptData, paths = 5e3, n = 200)
save(bssmatrix, file= 'spx20110915bss.rData')

# svifitSqrt <- sviSqrtFit(spxOptData)
# svifitQR<- sviFitQR(ivolData = spxOptData, sviGuess = svifitSqrt ) 

load('spx20110915bss.rData')

plotIvols_svibss(ivolData = spxOptData, paths= 1e5, n= 200, bssplotscale = .6, 
                 kappa = 1, bssMatrix = bssmatrix$best  )


# 
# # ### Estimate the bss params based on pure fitting optimization  
# #bssfit0 <- bssFit(spxOptData, paths, n, kappa=1)
# ###
# 

# t0<- proc.time()
# 
# bssMatrix<- bss_histfit(ivolData = spxOptData, OxfordManRVData = OxfordManRVRaw.data, 
#             symbol = symbol, EndDate = endDate, paths = paths, n = n)
# 
# proc.time()-t0

# 
# bssparam<- list(S0=1, xi= 0.235^2, eta= rvParam$eta, alpha= rvParam$h-.5, rho=-.9)
# plotIvols_svibss(ivolData=spxOptData, paths=paths, n=400, 
#               bssMatrix=bssMatrix,  kappa=1, sviMatrix=svifitQR)
# 


# params <- list(S0=1, xi=0.235^2, eta=1.9, alpha=-0.43, rho=-0.9)
# 
# finalPrices<- hybridScheme(params)(1e4, 400, 1, 0.041)
# 
# k<-5
# clr<- rainbow(k)
# finalP<- list()
# 
# 
# for (i in 1:k){
#   
#   param<- params
#   param$xi<- (sqrt(param$xi)+0.04*i)^2
#   #print(i)
#   
#   finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
# }
# save(finalP, file='finalP_xi.rData')
# 
# curve(ImpliedVol(x, finalPrices), from=-0.4, to=0.2, col='black',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0.1, 0.sm5) ,lwd=2)
# print(ImpliedVol(0, finalPrices))
# 
# for (i in 1:k){
#   curve(ImpliedVol(x, finalP[[i]]),from=-0.15, to=0.1, add=TRUE, col=clr[i], lwd=2)
#   print(ImpliedVol(0, finalP[[i]]))
# }
# 
# legend( 'topright',
#         legend=c ('xi= 0.235^2', 
#                   'xi= 0.245^2',
#                   'xi= 0.255^2',
#                   'xi= 0.265^2',
#                   'xi= 0.275^2',
#                   'xi= 0.285^2'), 
#         col = c('black', clr), lty=1)
# 
# 
# title( 'Sensitivity on Xi')



