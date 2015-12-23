
library(quantmod)
assetPrice<- function (symbol, from=NULL, to=NULL){
  return (get(getSymbols(symbol,from= from, to=to)))
}

assetDailyReturn<- function (symbol, from, to){
  # return daily return 
  
  price<- Cl(assetPrice(symbol, from, to))
  
  return (diff(log(price)))
}

ret<-assetDailyReturn('^GSPC', from='2000-01-01', to='2010-01-01')
head(ret)