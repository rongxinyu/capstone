rm(list=ls())
library(MASS)
library(gsl)
library(ggplot2)
Myconvolve <- function(x, y, length){
  return(sum(x[1:length] * rev(y[1:length])))
}

bstar <- function(k, alpha){
  result <- ((k^(alpha+1) - (k-1)^(alpha+1)) / (alpha+1)) ^ (1/alpha)
  return(result)
}

BSImpliedVolCall <- function(S0, K, T, r, C)
{
  nK <- length(K);
  sigmaL <- rep(1e-10,nK);
  CL <- BSFormula(S0, K, T, r, sigmaL);
  sigmaH <- rep(10,nK);
  CH <- BSFormula(S0, K, T, r, sigmaH);
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    CM <- BSFormula(S0, K, T, r, sigma);
    CL <- CL + (CM < C)*(CM-CL);
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL);
    CH <- CH + (CM >= C)*(CM-CH);
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH);
  }
  return(sigma);
  
}

MyImpliedVolCall <- function(S0, K, T, r, C){
  sigmaL <- 1e-3
  sigmaH <- 10
  while(sigmaH-sigmaL > 1e-10){
    sigma <- (sigmaL + sigmaH)/2
    if(BSFormula(S0, K, T, r, sigma)==C) {break}
    else if(BSFormula(S0, K, T, r, sigma)>C) {sigmaH <- sigma}
    else{
      sigmaL <- sigma
    }
  }
  return(sigma)
}

BSFormula <- function(S0, K, T, r, sigma)
{
  x <- log(S0/K)+r*T;
  sig <- sigma*sqrt(T);
  d1 <- x/sig+sig/2;
  d2 <- d1 - sig;
  pv <- exp(-r*T);
  return( S0*pnorm(d1) - pv*K*pnorm(d2));
}

hybridScheme <- function(params){
  
  S0 <- params$S0
  xi <- params$xi
  eta <- params$eta
  alpha <- params$alpha
  rho <- params$rho
  
  #n = 100
  #kappa = 1
  covMatrix <- function(n, kappa){
    ### generate the covariance matrix
    ### @param n: the distance between two points is 1/n
    ### @param kappa: how many power low terms to include around zero
    
    sigma <- matrix(0, nrow=kappa+1, ncol=kappa+1)
    sigma[1,1] <- 1/n
    for(j in 2:(kappa+1)){
      sigma[1, j] <- ((j-1)^(alpha+1)-(j-2)^(alpha+1)) / ((alpha+1)*n^(alpha+1))
      sigma[j, 1] <- ((j-1)^(alpha+1)-(j-2)^(alpha+1)) / ((alpha+1)*n^(alpha+1))
    }
    for(j in 2:(kappa+1)){
      for(k in 2:(kappa+1)){
        if(j==k){
          sigma[j,k] <- ((j-1)^(2*alpha+1)-(j-2)^(2*alpha+1)) / ((2*alpha+1)*n^(2*alpha+1))
        }
        else{
          #print(((j-1)(k-1))^(alpha+1) * f21hyper(1,2*(alpha+1),alpha+2,(k-1)/(k-j))) #- 
          #print(j)
          #print(k)
          #print((k-1)/(k-j))
          #print(hyperg_2F1(1,2*(alpha+1),alpha+2,(k-1)/(k-j)))
          #sigma[j,k] <- ((j-1)(k-1))^(alpha+1) * f21hyper(1,2*(alpha+1),alpha+2,(k-1)/(k-j)) #- 
          #               ((j-2)(k-2))^(alpha+1) * f21hyper(1,2*(alpha+1),alpha+2,(k-2)/(k-j))) /
          #               (j-k) * (alpha+1) * n^(2*alpha+1)  
          
        }
      }
    }
    #print(sigma)
    return(sigma)
  }
  
  Simulation <- function(n, kappa, T){
    W <- mvrnorm(floor(n*T), mu=rep(0, kappa+1), Sigma=covMatrix(n, kappa))
    Wperp <- rnorm(floor(n*T), sd=sqrt(1/n))
    Z <- rho * W[,1] + sqrt(1-rho*rho)*Wperp
    Gamma <- sapply(seq(1:floor(n*T)), function(x){(bstar(x, alpha)/n)^alpha}) 
    Gamma[1:kappa] <- 0
    Y2 <- sapply(seq(1,floor(n*T)), function(x){Myconvolve(Gamma, W[,1], x)})
    Y1 <- rep(0, floor(n*T))
    
    for(i in 1:floor(n*T)){
      Y1[i] <-0
      for (k in 1:min(i,kappa)){
        Y1[i] = Y1[i] + W[i+1-k,k+1]
      }
    }
    Y <- Y1+Y2 ## The simulated series of main interst
    v <- xi*exp(eta*sqrt(2*alpha+1)*Y - eta*eta/2*sapply(seq(1:floor(n*T)),function(x){(x/n)^(2*alpha+1)}))
    v <- c(xi, v[1:length(v)-1])
    S <- S0 * exp(sum(v^0.5*Z) - 1/2*sum(v)/n)
    return(S)
  }
  
  MC <- function(N, n, kappa, T){
    return(sapply(rep(T, N), function(x){Simulation(n, kappa, x)}))
  }
  return(MC)
}

impvol <- function(k, st, T){
  payoff <- (st > exp(k)) * (st - exp(k))
  return(BSImpliedVolCall(1, exp(k), T, 0, mean(payoff)))
}
params <- list(S0=1, xi=0.235^2, eta=1.9, alpha=-0.43, rho=-0.9)
T <- 0.1
n <- 1000 # step length is 1/n
paths <- 10000
set.seed(9081)
finalPrices <- hybridScheme(params)(paths, n, 1, T)
summary(finalPrices)

vol <- function(k){sapply(k, function(x){impvol(x, finalPrices, T)})}
payoff <- function(k){sapply(k, function(x){mean((finalPrices > exp(x)) * (finalPrices - exp(x)))})}
curve(vol(x), from=-0.5, to=0.5, xlab="Log strike", ylab="Implied Vol")
