# Fit BSS independently to each slice

source('bssIv.R')
bssSqFit<- function (ivolData, paths=5e3, n=400 , kappa=1, kscale= .4,
                     guess=list(S0=1, eta=2 , xi=0.24^2, alpha=-0.44, rho=-.9 ) ,slices= NULL){
  bidVols <- as.numeric(ivolData$Bid)
  askVols <- as.numeric(ivolData$Ask)
  expDates <- unique(ivolData$Texp)
  nSlices <- length(expDates)
  if (!is.null(slices)) {
    nSlices <- length(slices)
  }
  else {
    slices <- 1:nSlices
  }
  colnum <- sqrt(nSlices * 2)
  rows <- round(colnum/2, 0)
  columns <- round(colnum, 0)
  while (rows * columns < nSlices) {
    rows <- rows + 1
  }
  atmVol <- numeric(nSlices)
  atmSkew <- numeric(nSlices)
  
  S0<- rep(1, nSlices)
  alpha<- numeric(nSlices)
  eta<-   numeric(nSlices)
  rho<-   numeric(nSlices)
  xi<-    numeric(nSlices)
  
  TEXP<- expDates[slices]
  alpha1<- numeric(nSlices)
  eta1<-   numeric(nSlices)
  rho1<-   numeric(nSlices)
  xi1<-    numeric(nSlices)
  
  
  
#   par(mfrow = c(rows, columns), mex = 0.5)
  

  i<-1 
  
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
#       finalPrices<- hybridScheme(bssMatrix[slice,])(paths,n,kappa, t)
#       print(bssMatrix[slice,])
#       curve(ImpliedVol(x, finalPrices, t), from = kmin, to = kmax, col='dark green', add=TRUE)
#     }
#     

    kIn <- k[!is.na(midVol)]
    volIn <- midVol[!is.na(midVol)]
    volInterp <- function(xout) {
      stinterp(x = kIn, y = volIn, xout)$y
    }
    atmVol[i] <- volInterp(0)
#    atmSkew[slice] <- (volInterp(0.01) - volInterp(-0.01))/0.02
    N<-n;
    if (t<0.1) { N<- max(N,400)}
    xi_obj<- function (param_xi, param, k_in, midVol_in){
#       k_low<- kmin*.3
#       k_high<- kmax*.3
      
      param$xi<-param_xi
      
#       idx<- which(kIn> k_low & kIn<k_high)
#       midVol_in<- volIn[idx]
#       k_in<- kIn[idx]
      
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
    }
    
    eta_obj<- function (param_eta, param, k_in, midVol_in){
#       k_low<- kmin*.3
#       k_high<- kmax*.3
      
      param$eta<- param_eta
      
#       idx<- which(kIn> k_low & kIn<k_high)
#       midVol_in<- volIn[idx]
#       k_in<- kIn[idx]
      
      
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
      
    }
    
    rho_obj<-function(param_rho, param, k_in, midVol_in){
#       k_low<- kmin*.3
#       k_high<- kmax*.3
      
      param$rho<- param_rho
      
#       idx<- which(kIn> k_low & kIn<k_high)
#       midVol_in<- volIn[idx]
#       k_in<- kIn[idx]
      
      bssVol<- bssiv(bssparams = param, paths = paths, n = N, kappa = kappa, k = k_in, t = t )
      
      res<- sum((bssVol- midVol_in)^2) *1e4
      
      return (res)
      
    }
    
    bssparam<- guess
    
   # kscale<- 0.3
    k_low<- kmin*kscale
    k_high<- kmax*kscale
    idx<- which(kIn> k_low & kIn<k_high)
    midVol_in<- volIn[idx]
    k_in<- kIn[idx]
    
    ## First, optimize xi
    set.seed(9081)
    lower<-(atmVol[i]-0.1)^2
    upper<- (atmVol[i]+0.1)^2
    fit_xi<- optimize(f = xi_obj, param= bssparam, k_in= k_in, midVol_in= midVol_in, 
                      lower = lower, upper = upper, maximum = F)
    print(fit_xi)
    if(fit_xi$minimum<= lower + 0.0001 || fit_xi$minimum >= upper-0.0001){
      message("Warning! The boundary of Xi is hit.")
    }
    bssparam$xi<- fit_xi$minimum
    fit_xi$param<- bssparam
    
    ## Second, optimize eta, with optimized xi
    set.seed(9081)
    lower<- 1.6
    upper<- 2.8
    fit_eta<- optimize(f = eta_obj, param= bssparam,k_in=k_in, midVol_in=midVol_in, 
                       lower = lower, upper = upper, maximum = F)
    print(fit_eta)
    if(fit_eta$minimum<=lower+0.0001 || fit_eta$minimum >= upper-0.0001){
      message('Warning! The boundary of Eta is hit.')
    }
    bssparam$eta<- fit_eta$minimum
    fit_eta$param<- bssparam
    
    ## Last, optimize rho, with optimized xi and eta
    set.seed(9081)
    lower<- -0.99
    upper<--.5
    
    fit_rho<- optimize(f= rho_obj, param= bssparam, k_in= k_in, midVol_in=midVol_in, 
                       lower = lower, upper= upper, maximum = F)
    print(fit_rho)
    if(fit_rho$minimum< lower+0.0001 || fit_rho$minimum > upper-0.0001){
      message('Warning! The boundary of Rho is hit.')
    }
    bssparam$rho<- fit_rho$minimum
    fit_rho$param<- bssparam
    
    if (!(fit_xi$objective>= fit_eta$objective/3 && fit_eta$objective>= fit_rho$objective/3)){
      message ('The optimization kernal shows suspicious contradicts with t= ', t)
    }
    
    alpha[i]<- bssparam$alpha
    rho[i]  <- bssparam$rho
    eta[i]  <- bssparam$eta
    xi[i]   <- bssparam$xi
    
    bestfit<- fit_xi
    if (fit_eta$object< bestfit$object) {bestfit<- fit_eta}
    if (fit_rho$object< bestfit$object) {bestfit<- fit_rho}
    
    alpha1[i]<- bestfit$param$alpha
    rho1[i]<-   bestfit$param$rho
    eta1[i]<-   bestfit$param$eta
    xi1[i]<-    bestfit$param$xi
    
    i<- i+1
  }
#   par(mfrow = c(1, 1), mex = 1)
#   par(new = F)
#  return(list(expiries = expDates, atmVol = atmVol, atmSkew = atmSkew))
  res<- list(normal= data.frame(texp= TEXP, S0=S0, xi= xi, alpha= alpha, eta= eta, rho=rho),
             best= data.frame(texp= TEXP, S0= S0, xi= xi1, alpha= alpha1, eta= eta1, rho=rho1))
  return (res)
  
}



bssFit <- function(ivolData, paths, n, kappa=1){
  
  bidVols <- as.numeric(ivolData$Bid);
  askVols <- as.numeric(ivolData$Ask);
  expDates <- unique(ivolData$Texp);
  nSlices <- length(expDates);
  
  slices <- 1:nSlices;
  bssMatrix <- array(dim=c(nSlices,5));

  ######################################
  for (slice in slices){
    t <- expDates[slice];
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0);
    pick <- !is.na(pick);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
    k <- log(ivolData$Strike[texp==t]/f)[pick]; # Plot vs log-strike
    
    bssGuess <- list(S0=1)
    # bssGuess <- list(S0=1, xi = mean(midVar,na.rm=T), eta= 1.9,  rho = -0.5, alpha=-0.43); 
    # bssGuess$xi<- 0.235^2
    
    # Define objective function
    obj <- function(bssparams){
      # tmp1 <- as.list(bssparams);
      tmp1 <- list(S0=1.0, xi = mean(midVar,na.rm=T), eta= bssparams,  rho = -0.2064355, alpha=-0.43); 
      names(tmp1) <- c("S0","xi","eta","rho","alpha")
      bssVar<- bss(tmp1, k, paths, n, kappa, t);
      tmp <- sum((midVar-bssVar)^2, na.rm=T);
      return(tmp*1e4);
    };
   
    # fit <- optim(bssGuess, obj);
    fit <- optimize(obj, c(1, 3))

    # tmp1 <- as.list(fit$par)
    tmp1 <- list(S0=1.0, xi = mean(midVar,na.rm=T), eta= fit$minimum,  rho = -0.2064355, alpha=-0.43); 
    
    # bssMatrix[slice,] <- tmp1 * c(1,1,1,1,1);
    bssMatrix[slice, ] <- unlist(tmp1);
    print(bssMatrix)

  }# end of for loop
  ######################################
  
  bssMatrix <- as.data.frame(bssMatrix);
  colnames(bssMatrix) <- c("S0","xi","eta","rho","alpha");
  return(bssMatrix);
}# end of function
