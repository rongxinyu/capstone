# Fit BSS independently to each slice

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
