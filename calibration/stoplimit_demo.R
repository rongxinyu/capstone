> # define the strategy
  > strategy(strat.st, store=TRUE)

> #one indicator
  > add.indicator(strat.st, name = "MACD", 
                  +   		  arguments = list(x=quote(Cl(mktdata))),
                  + 			  label='_' 
                  + )
[1] "macd"

> #two signals
  > add.signal(strat.st,name="sigThreshold",
               + 		   arguments = list(column="signal._",
                                       + 				   			relationship="gt",
                                       + 							threshold=0,
                                       + 							cross=TRUE),
               + 		   label="signal.gt.zero"
               + )
[1] "macd"

> add.signal(strat.st,name="sigThreshold",
             + 		   arguments = list(column="signal._",
                                     + 				            relationship="lt",
                                     + 							threshold=0,
                                     + 							cross=TRUE),
             + 	       label="signal.lt.zero"
             + )
[1] "macd"

> ####
  > # add rules
  > 
  > # entry
  > add.rule(strat.st,name='ruleSignal', 
             + 		 arguments = list(sigcol="signal.gt.zero",
                                   + 				         sigval=TRUE, 
                                   + 						 orderqty=100, 
                                   + 						 ordertype='market', 
                                   + 						 orderside='long', 
                                   + 						 threshold=NULL),
             + 	     type='enter',
             + 		 label='enter',
             + 		 storefun=FALSE
             + )
[1] "macd"

> #alternatives for risk stops:
  > # simple stoplimit order, with threshold multiplier
  > #add.rule(strat.st,name='ruleSignal', arguments = list(sigcol="signal.gt.zero",sigval=TRUE, orderqty='all', ordertype='stoplimit', orderside='long', threshold=-.05,tmult=TRUE, orderset='exit2'),type='chain', parent='enter', label='risk',storefun=FALSE)
  > # alternately, use a trailing order, also with a threshold multiplier
  > #add.rule(strat.st,name='ruleSignal', arguments = list(sigcol="signal.gt.zero",sigval=TRUE, orderqty='all', ordertype='stoptrailing', orderside='long', threshold=-1,tmult=FALSE, orderset='exit2'),	type='chain', parent='enter', label='trailingexit')
  > 
  > # exit
  > add.rule(strat.st,name='ruleSignal', 
             + 		 arguments = list(sigcol="signal.lt.zero",
                                   + 				          sigval=TRUE, 
                                   + 						  orderqty='all', 
                                   + 						  ordertype='market', 
                                   + 						  orderside='long', 
                                   + 						  threshold=NULL,
                                   + 						  orderset='exit2'),
             +          type='exit',
             + 		 label='exit'
             + )
[1] "macd"

> #end rules
  > ####
  > 
  > getSymbols(stock.str,from=initDate)
[1] "AAPL"

> start_t<-Sys.time()

> out<-applyStrategy(strat.st , portfolios=portfolio.st,parameters=list(nFast=fastMA, nSlow=slowMA, nSig=signalMA,maType=maType),verbose=TRUE)