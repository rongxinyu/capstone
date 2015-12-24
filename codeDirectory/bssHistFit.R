
library(stinepack)
source('OxfordManRV.R')
source('hybrid_bss.R')

bss_histfit<- function (ivolData, OxfordManRVData, symbol, EndDate, 
                        paths, n , 
                        rho_guess=-.3, rho_up=-.1, rho_low=-.9, kappa=1){
  rv.data<- OxfordManRV.clean(rv.data = OxfordManRVData)
  rv.param<- OxfordManRV.param(rv.list = rv.data, symbol = symbol, endDate = EndDate)
  h<- rv.param$h
  alpha<- h-.5
  eta<- rv.param$eta
  ## now we compute alpha and eta.
  ## The following is to calibrate xi and rho. 
  ## Note that xi is the atmvol. We optimize rho. 
  
  rho_xi<- rho_optim (ivolData, paths, n, alpha, eta, rho_guess, rho_up, rho_low, kappa)
  ## rho_optim returns a list(rho, atmvol, texp, atmskew)
  
  slices<- length(rho_xi$texp)
  bssMatrix<- data.frame(texp=rho_xi$texp, S0=rep(1, slices), alpha=rep(alpha, slices), 
                          eta=rep(eta,slices), xi=rho_xi$atmvol^2, rho=rho_xi$rho) 


  return (bssMatrix)
}

rho_optim<- function (ivolData, paths, n, alpha, eta, rho_guess, rho_up, rho_low, kappa=1){
  
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  slices<- 1:nSlices
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  rho<-numeric(nSlices)
  par(mfrow = c(rows, columns), mex = 0.5)
  

  for (slice in slices) {
    t <- expDates[slice]
    texp <- ivolData$Texp
    bidVol <- bidVols[texp == t]
    askVol <- askVols[texp == t]
    midVol <- (bidVol + askVol)/2
    f <- (ivolData$Fwd[texp == t])[1]
    k <- log(ivolData$Strike[texp == t]/f)
    include <- !is.na(bidVol)
    kmin <- min(k[include])
    kmax <- max(k[include])
#     ybottom <- 0.8 * min(bidVol[include])
#     ytop <- 1.2 * max(askVol[include], na.rm = T)
#     xrange <- c(kmin, kmax)
#     yrange <- c(0.1, ytop)
#     plot(k, bidVol, col = "red", pch = 18, cex = 0.5, xlim = xrange, 
#          ylim = yrange, main = paste("T =", format(t, digits = 2, 
#                                                    nsmall = 2)), xlab = "Log-Strike", ylab = "Implied Vol.")
#     par(new = T)
#     plot(k, askVol, col = "blue", pch = 18, cex = 0.5, xlim = xrange, 
#          ylim = yrange, main = NA, xlab = NA, ylab = NA)
#     
#     
    
#     ### plot svi fit
#     if ((!is.null(sviMatrix))) {
#       vol <- function(k) {
#         sqrt(svi(sviMatrix[slice, ], k)/t)
#       }
#       par(new = T)
#       curve(vol(x), from = kmin, to = kmax, col = "orange", 
#             lwd = 2, add = T)
#     }
#     
#     
#     ### plot bss fit
#     if (! is.null(bssMatrix)){
#       finalPrices<- hybridScheme(bssMatrix)(paths,n,kappa, t)
#       curve(ImpliedVol(x, finalPrices), from = kmin, to = kmax, col='dark green', add=TRUE)
#     }
#     
    
    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[slice] <- volInterp(0)
    atmSkew[slice] <- (volInterp(0.01) - volInterp(-0.01))/0.02
    
    xi<- atmVol[slice]^2
    
    minimum_k<- optimize(volInterp, c(1,-1))$minimum
    
    guess=rho_guess
    if(slice>1) {guess= rho[slice-1]}
    
    bssparam<- list(S0=1, xi=xi, eta=eta, alpha=alpha, rho=guess )
    finalPrices<- hybridScheme(bssparam)(paths, n, kappa, t)
    vol0<- ImpliedVol(minimum_k, finalPrices)
    vol1<- ImpliedVol(minimum_k+0.01, finalPrices)
    
    if (vol1>vol0){
      while(vol1>vol0 & bssparam$rho>rho_low & bssparam$rho<rho_up){
        bssparam$rho<- bssparam$rho-0.01
        finalPrices<- hybridScheme(bssparam)(paths, n, kappa,t)
        vol0<- ImpliedVol(minimum_k, finalPrices)
        vol1<- ImpliedVol(minimum_k+0.01,finalPrices)
      }
    }else{
      while(vol1<= vol0 & bssparam$rho>rho_low & bssparam$rho<rho_up){
        bssparam$rho<- bssparam$rho+0.01
        finalPrices<- hybridScheme(bssparam) (paths, n, kappa, t)
        vol0<- ImpliedVol(minimum_k, finalPrices)
        vol1<- ImpliedVol(minimum_k+0.01, finalPrices)
      }
    }
    
    rho[slice]<- bssparam$rho
    
    print(rho[slice])
  }
  
#   par(mfrow = c(1, 1), mex = 1)
#   par(new = F)
  return(list(texp = expDates, atmvol = atmVol, atmskew = atmSkew, rho=rho))
  
  

}