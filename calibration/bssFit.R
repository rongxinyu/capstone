# Fit BSS independently to each slice

bssFit <- function(ivolData, paths, n, kappa=1){
  
  bidVols <- as.numeric(ivolData$Bid);
  askVols <- as.numeric(ivolData$Ask);
  expDates <- unique(ivolData$Texp);
  nSlices <- length(expDates);
  slices <- 1:nSlices;
  bssMatrix <- array(dim=c(nSlices,5));
 #print(nSlices)
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
    
    bssGuess <- list(S0=1, xi = mean(midVar,na.rm=T), eta= 1.9,  rho = - 0.9, alpha=-0.43); 
    bssGuess$xi<- 0.235^2
    
    # Define objective function
    obj <- function(bssparams){
      tmp1 <- as.list(bssparams);
      names(tmp1) <- c("S0","xi","eta","rho","alpha")
      #sviVar <- svi(tmp1,k);
      bssVar<- bss(tmp1, paths, n, kappa, t)(k);
      #minVar <- tmp1$a+tmp1$b*tmp1$sig*sqrt(abs(1-tmp1$rho^2));
      #negVarPenalty <- min(100,exp(-1/minVar));
      tmp <- sum((midVar-bssVar)^2,na.rm=T);
      return(tmp*100000);
    };
   
    #fit <- optim(bssGuess,obj);
 
    #if(abs(fit$par["rho"])>.999){
      fit <- optim(bssGuess,obj,method="Nelder-Mead", hessian = FALSE
                 # lower=bssGuess,
                  # upper=bssGuess
                 );
    print(t)
                   #upper=c(1,0.5^2,10,+.999,0.0000001));
      #}

    bssMatrix[slice,] <- fit$par*c(1,1,1,1,1);##########

  }# end of for loop
  ######################################
  bssMatrix <- as.data.frame(bssMatrix);
  colnames(bssMatrix) <- c("S0","xi","eta","rho","alpha");
#   bssGuess <- list(S0=1, xi = mean(midVar,na.rm=T), eta= 1.9,  rho = - 0.9, alpha=-0.43); 
#   bssMatrix<- data.frame(bssGuess)
#   for (i in 2:8){
#     bssMatrix<- rbind(bssMatrix, bssGuess)
#   }   
  return(bssMatrix);
}# end of function
