## JPM SCTO Strategy

## This Strategy can beat SPY only when there is market crush like 2008 in the life period. Or SPY is has better return.



##Clean the environment 
rm(list=ls())

## Start 
require(quantmod)
require(PerformanceAnalytics)

ini.date="2000-01-01"

## download data and calculate the return
symbols <- c("XLB", "XLE", "XLF", "XLI", "XLK", "XLP", "XLU", "XLV", "XLY", "RWR", "SHY")
getSymbols(symbols, from=ini.date)
prices <- list()
for(i in 1:length(symbols)) {
  prices[[i]] <- Ad(get(symbols[i]))  
}
prices <- do.call(cbind, prices)
colnames(prices) <- gsub("\\.[A-z]*", "", colnames(prices))
returns <- na.omit(Return.calculate(prices))  ## calculate the returns (daily)



sctoStrat <- function(returns, cashAsset = "SHY", lookback = 4, annVolLimit = .2,
                      topN = 5, scale = 252) {
  ep <- endpoints(returns, on = "months")
  weights <- list()
  cashCol <- grep(cashAsset, colnames(returns))
  
  #remove cash from asset returns
  cashRets <- returns[, cashCol]
  assetRets <- returns[, -cashCol]
  for(i in 2:(length(ep) - lookback)) {
    retSubset <- assetRets[ep[i]:ep[i+lookback]]
    
    #forecast is the cumulative return of the lookback period
    forecast <- Return.cumulative(retSubset)
    
    #annualized (realized) volatility uses a 22-day lookback period
    annVol <- StdDev.annualized(tail(retSubset, 22))
    
    #rank the forecasts (the cumulative returns of the lookback)
    rankForecast <- rank(forecast) - ncol(assetRets) + topN
    
    #weight is inversely proportional to annualized vol
    weight <- 1/annVol
    
    #zero out anything not in the top N assets
    weight[rankForecast <= 0] <- 0
    
    #normalize and zero out anything with a negative return
    weight <- weight/sum(weight)
    weight[forecast < 0] <- 0
    
    #compute forecasted vol of portfolio
    forecastVol <- sqrt(as.numeric(t(weight)) %*% 
                          cov(retSubset) %*% 
                          as.numeric(weight)) * sqrt(scale)
    
    #if forecasted vol greater than vol limit, cut it down
    if(as.numeric(forecastVol) > annVolLimit) {
      weight <- weight * annVolLimit/as.numeric(forecastVol)
    }
    weights[[i]] <- xts(weight, order.by=index(tail(retSubset, 1)))
  }

  #replace cash back into returns
  returns <- cbind(assetRets, cashRets)
  weights <- do.call(rbind, weights)
  
  #cash weights are anything not in securities
  weights$CASH <- 1-rowSums(weights)
  
  #compute and return strategy returns
  stratRets <- Return.portfolio(R = returns, weights = weights)
  return(stratRets)      
}

scto4_20 <- sctoStrat(returns)
getSymbols("SPY", from = ini.date)
spyRets <- Return.calculate(Ad(SPY))
comparison <- na.omit(cbind(scto4_20, spyRets))
colnames(comparison) <- c("strategy", "SPY")
charts.PerformanceSummary(comparison)
apply.yearly(comparison, Return.cumulative)
stats <- rbind(table.AnnualizedReturns(comparison),
               maxDrawdown(comparison),
               CalmarRatio(comparison),
               SortinoRatio(comparison)*sqrt(252))
round(stats, 3)
