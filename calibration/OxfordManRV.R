library(xts)


OxfordManRV.clean<- function (rv.data){
  # clean and regularize the realized vol data from Oxford-Man Realized Vol Institution
  
  # There are many different estimates of realized variance, all of them very similar.
  # We keep the realized kernel estimates denoted by ".rk", and delete others. 
  
  colnums <- which(sapply(rv.data, function(x) grep(".rk",x))>0)
  col.names <- names(colnums)
  
  rv1 <- rv.data[,col.names]
  index.names <- rv1[2,]
  
  datesRaw <- rv.data[-(1:2),1]
  dates <- strptime(datesRaw,"%Y%m%d")
  
  rv.list <- NULL
  index.names <- as.matrix(index.names)
  
  n <- length(index.names)
  for (i in 1:n){
    tmp.krv1 <- xts(as.numeric(rv1[-(1:2),i]),order.by=dates) 
    rv.list[[i]] <- tmp.krv1[(tmp.krv1!="")&(tmp.krv1!="0")]
  }
  
  names(rv.list)<- index.names 
  
  save(rv.list, file="oxfordRV.rData")
  
  return (rv.list)
}


OxfordManRV.param<- function ( rv.list,symbol, endDate){

  rv<- rv.list[[symbol]]
  rv<- rv[index(rv)< endDate]
  sig<- sqrt(as.numeric(rv))
  

  x <- 1:100
  dlsig2 <- function(lag){mean((diff(log(sig),lag=lag))^2)}
  dlsig2Vec <- function(x){sapply(x,dlsig2)}
  
  fit.lm <- lm(log(dlsig2Vec(x)) ~ log(x))
  
  nu <- sqrt(exp(coef(fit.lm)[1]))
  h <- coef(fit.lm)[2]/2
  
  eta<- 2*nu* sqrt(gamma(1.5-h)/(gamma(h+.5)*gamma(2-2*h)) )
  return (list(h= as.numeric(h), eta= as.numeric(eta), nu= as.numeric(nu), sig= rv))  
  
  
}