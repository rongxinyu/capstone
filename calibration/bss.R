
bss<- function (bssparams, paths, n, kappa,t){
  
  res<-function (k){
  finalPrices<- hybridScheme(bssparams)(paths,n,kappa, t)
  return(ImpliedVol(k, finalPrices)^2*t)
  }  
  return (res)
}
