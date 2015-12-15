rm(list=ls())
library(MASS)
library(parallel)

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
    if(kappa==0){
      return (sigma)
    }
    else{
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
  }
  
  Simulation <- function(n, kappa, T, W, Z, Gamma, tseq){
    #print(W)
    
    Y2 <- convolve(Gamma, rev(W[,1]), type="open")[1:floor(n*T)]
    Y1 <- rep(0, floor(n*T))
    
    for(i in 1:floor(n*T)){
      Y1[i] <-0
      if (kappa!=0 ){
      for (k in 1:min(i,kappa)){
        Y1[i] = Y1[i] + W[i+1-k,k+1]
      }
      }
    }
    Y <- Y1+Y2 ## The simulated series of main interst
    v <- xi*exp(eta*sqrt(2*alpha+1)*Y - tseq)
    v <- c(xi, v[1:length(v)-1])
    S <- S0 * exp(sum(v^0.5*Z) - 1/2*sum(v)/n)
    return(S)
  }
  
  MC <- function(N, n, kappa, T){
    Gamma <- sapply(seq(1:floor(n*T)), function(x){(bstar(x, alpha)/n)^alpha}) 
    Gamma[1:kappa] <- 0
    tseq <- eta*eta/2*sapply(seq(1:floor(n*T)),function(x){(x/n)^(2*alpha+1)})
    steps <- floor(n*T)
    W <- mvrnorm(steps*N, mu=rep(0, kappa+1), Sigma=covMatrix(n, kappa))
    Wperp <- rnorm(steps*N, sd=sqrt(1/n))
    Z <- rho * W[,1] + sqrt(1-rho*rho)*Wperp
    return(sapply(seq(1:N), function(loopNum){Simulation(
      n, kappa, T, W[(1+(loopNum-1)*steps):(loopNum*steps),],Z[(1+(loopNum-1)*steps):(loopNum*steps)], Gamma, tseq)}))
  }
  return(MC)
}

impvol <- function(k, st, T){
  payoff <- (st > exp(k)) * (st - exp(k))
  return(BSImpliedVolCall(1, exp(k), T, 0, mean(payoff)))
}




params <- list(S0=1, xi=0.235^2, eta=1.9, alpha=-0.43, rho=-0.9)
T <- 1
n <- 20 # step length is 1/n
paths <- 1e5
set.seed(9081)

Rprof()
finalPrices <- hybridScheme(params)(paths, n, 1, T)
Rprof(NULL)
summaryRprof()

summary(finalPrices)



vol <- function(k, finalP){sapply(k, function(x){impvol(x, finalP, T)})}
payoff <- function(k, finalP){sapply(k, function(x){mean((finalP > exp(x)) * (finalP - exp(x)))})}


# 
# # param2<- list(S0=1, xi=0.7, eta=1.9, alpha= -0.43, rho= -0.9)
# finalPrices2<- hybridScheme(params)(paths, n, 0, T)
# 
# save(finalPrices, file='finalPT1kappa1.rData')
# save(finalPrices2,file ='finalPT1kappa0.rData')
# 
# 
# curve(vol(x, finalPrices), from=-0.3, to=0.2, col='blue',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0, 0.55))
# curve(vol(x, finalPrices2),from=-0.18, to=0.08, add=TRUE, col='red')
# legend('topright', legend= c('kappa= 0', 'kappa= 1'), col= c('red', 'blue'), lty=1)
# title(main = 'T=1')
# 
# 

# 
# k<-5
# clr<- rainbow(k)
# finalP<- list()
# 
# 
# for (i in 1:k){
#   
#   param<- params
#   param$xi<- param$xi+0.01*i
#   print(i)
#   
#   finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
# }
# save(finalP, file='finalP_xi.rData')
# 
# curve(vol(x, finalPrices), from=-0.15, to=0.15, col='black',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0.1, 0.55) ,lwd=2)
# 
# for (i in 1:k){
#   curve(vol(x, finalP[[i]]),from=-0.15, to=0.15, add=TRUE, col=clr[i], lwd=2)
# }
# 
# 
# 
# legend( 'topright',
#         legend=c ('xi= 0.235^2', 
#                   'xi= 0.235^2+ 0.01',
#                   'xi= 0.235^2+ 0.02',
#                   'xi= 0.235^2+ 0.03',
#                   'xi= 0.235^2+ 0.04',
#                   'xi= 0.235^2+ 0.05'), 
#         col = c('black', clr), lty=1)
# 
# 
# title( 'Sensitivity on Xi')




# 
# k<-3
# clr<- rainbow(k)
# finalP<- list()
# 
# 
# for (i in 1:k){
#   
#   param<- params
#   param$eta<- param$eta+0.3*i
#   print(i)
#   
#   finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
# }
# save(finalP, file='finalP_eta.rData')
# 
# curve(vol(x, finalPrices), from=-0.15, to=0.15, col='black',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0.1, 0.55) ,lwd=2)
# 
# for (i in 1:k){
#   curve(vol(x, finalP[[i]]),from=-0.15, to=0.15, add=TRUE, col=clr[i], lwd=2)
# }
# 
# 
# 
# legend( 'topright',
#         legend=c ('eta= 1.9', 
#                   'eta= 2.2',
#                   'eta= 2.5',
#                   'eta= 2.8'), 
#         col = c('black', clr), lty=1)
# 
# 
# title( 'Sensitivity on Eta')
# 
# 
# 
# 
# 
# 
# 
# k<-4
# clr<- rainbow(k)
# finalP<- list()
# 
# 
# for (i in 1:k){
#   
#   param<- params
#   param$alpha<- param$alpha+0.02*i
#   print(i)
#   
#   finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
# }
# save(finalP, file='finalP_alpha.rData')
# 
# curve(vol(x, finalPrices), from=-0.15, to=0.1, col='black',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0.1, 0.55) ,lwd=2)
# 
# for (i in 1:k){
#   curve(vol(x, finalP[[i]]),from=-0.18, to=0.1, add=TRUE, col=clr[i], lwd=2)
# }
# 
# 
# 
# legend( 'topright',
#         legend=c ('alpha= -0.43', 
#                   'alpha= -0.41',
#                   'alpha= -0.39',
#                   'alpha= -0.37',
#                   'alpha= -0.35'), 
#         col = c('black', clr), lty=1)
# 
# 
# title( 'Sensitivity on Alpha')









# 
# k<-4
# clr<- rainbow(k)
# finalP<- list()
# 
# 
# for (i in 1:k){
#   
#   param<- params
#   param$alpha<- param$alpha+0.1*i
#   print(i)
#   
#   finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
# }
#   
# save(finalP, file='finalP_rho.rData')
# 
# curve(vol(x, finalPrices), from=-0.1, to=0.1, col='black',
#       xlab="Log strike", ylab="Implied Vol", ylim= c(0.1, 0.35) ,lwd=2)
# 
# for (i in 1:k){
#   curve(vol(x, finalP[[i]]),from=-0.1, to=0.1, add=TRUE, col=clr[i], lwd=2)
# }
# 
# 
# 
# legend( 'topright',
#         legend=c ('rho= -0.90', 
#                   'rho= -0.80',
#                   'rho= -0.70',
#                   'rho= -0.60',
#                   'rho= -0.50'), 
#         col = c('black', clr), lty=1)
# 
# 
# 
# 
# title( 'Sensitivity on Rho')





k<-5
clr<- rainbow(k)
finalP<- list()

#params$rho<-0
#finalPrices <- hybridScheme(params)(paths, n, 1, T)
#summary(finalPrices)

#for (i in 1:k){
#  
#  param<- params
#  param$alpha<- param$alpha+0.1*i
#  print(i)
  
#  finalP[[i]]<- hybridScheme(param)(paths, n,1,T)
#}

#save(finalP, file='finalP_rho.rData')

curve(vol(x, finalPrices), from=-0.35, to=0.35, col='black',
      xlab="Log strike", ylab="Implied Vol", ylim= c(0.15, 0.4) ,lwd=2)

#for (i in 1:(k-1)){
#  curve(vol(x, finalP[[i]]),from=-0.35, to=0.15, add=TRUE, col=clr[i], lwd=2)
#}



#legend( 'topleft',
#        legend=c ('rho= 0', 
#                  'rho= 0.10',
#                  'rho= 0.20',
 #                 'rho= 0.30',
#                  'rho= 0.40'), 
 #       col = c('black', clr), lty=1)




#title( 'Sensitivity on Rho')

