rm(list=ls())

book<- function(LL, logging){
  
  Price<- -LL:LL
  buySize<- c(rep(5,LL-8), 5,4,4,3,3,2,2,1, rep(0, LL+1))
  sellSize<- c(rep(0,LL), 0,1,2,2,3,3,4,4,5,rep(5, LL-8))
  book<- data.frame(Price, buySize, sellSize )
  book<-list(Book=book, Log=NA, logging=F)
  
  if(logging) {
    eventLog<- as.data.frame(matrix(0,nr= 0, nc=2))
    colnames(eventLog)<-c("Types", "Prices")
    
    book$Log<-eventLog
    book$logging=T
  }
  return (book)
}

bid<- function(obj){
  max(obj$Price[obj$buySize>0])
}

ask<- function(obj){
  min(obj$Price[obj$sellSize>0])
}

spread<- function(obj){
  ask(obj)-bid(obj)
}

midPrice<-function(obj){
  0.5*(bid(obj)+ask(obj))
}

bidDepth<-function(obj){
  length(obj$Price[obj$Price<=bid(obj)])
}

askDepth<- function(obj){
  length(obj$Price[obj$Price>=ask(obj)])
}

swithLog<-function(obj, status){ 
  obj$logging<<-status
  return (obj)
}

### Although we may have a very long book, 
### we only consider changes occuring near the bid-ask

pick<- function (L) sample(1:L, 1)

limitBuy<- function(obj, L, price=NA){
  if(is.na(price)) {
    bid= bid(obj$Book)
    price=bid-pick(L)+1
  }
  obj$Book $buySize[obj$Book$Price==price]<-
    obj$Book $buySize[obj$Book$Price==price]+1
 ################ 
  if(obj$logging) 
  {
    obj$Log[nrow(obj$Log)+1, ]= c("LB", price)
  }
  return(obj)
}

limitSell<- function (obj, L, price=NA){
  if(is.na(price)){
    price= ask(obj$Book)+pick(L)-1
  }
  obj$Book$sellSize[obj$Book$Price==price]<-
    obj$Book$sellSize[obj$Book$Price==price]+1
  
  if(obj$logging) {
    obj$Log[nrow(obj$Log)+1, ]= c("LS", price)
  }
  
  obj
}

marketBuy<-function(obj){
  price=ask(obj$Book)
  obj$Book$sellSize[obj$Book$Price==price]<-
    obj$Book$sellSize[obj$Book$Price==price]-1
  
  if(obj$logging) 
    {
    obj$Log[nrow(obj$Log)+1, ]= c("MB", price)
  }
  obj
  
}

marketSell<-function(obj){
  price=bid(obj$Book)
  obj$Book$buySize[obj$Book$Price==price]<-
    obj$Book$buySize[obj$Book$Price==price]-1
  
  if(obj$logging) {
    obj$Log[nrow(obj$Log)+1, ]= c("MS", price)
  }
  obj
}

### cancel order is a little tricky
cancelableBuy<-function (obj, L){
  sum(obj$buySize[obj$Price>bid(obj)-L])
}

cancelableSell<-function (obj, L){
  sum(obj$sellSize[obj$Price<ask(obj)+L])
}

cancelBuy<-function(obj, L, price=NA){
  if(is.na(price)){
    nb= cancelableBuy (obj$Book, L)
    a=pick(nb)
    tmp= cumsum(rev(obj$Book$buySize))
    p=length(tmp[tmp>=a])
    price=obj$Book$Price[p]
  }
  
  obj$Book$buySize[p]<-obj$Book$buySize[p]-1
  
  if(obj$logging)
    {
    obj$Log[nrow(obj$Log)+1, ]= c("CB", price)
  }
  
  obj
}

cancelSell<-function (obj, L, price =NA){
  if (is.na(price)){
    ns=cancelableSell(obj$Book, L)
    a=pick(ns)
    tmp= cumsum(obj$Book$sellSize)
    p= length(tmp[tmp<a])+1
    price=obj$Book$Price[p]
  }
  
  obj$Book$sellSize[p]<-obj$Book$sellSize[p]-1
  
  if (obj$logging){
    obj$Log[nrow(obj$Log)+1, ]= c("CS", price)
  }
  obj
}


event<-function(obj, L, mu, alpha, delta){
  nb=cancelableBuy(obj$Book, L)
  ns=cancelableSell( obj$Book, L)
  prob=c(L*alpha, L*alpha, mu*.5, mu*.5, nb*delta, ns*delta)
  prob=prob/sum(prob)
 
  p= sample(1:6, 1, replace=T, prob=prob)
  switch(p,
         limitBuy(obj , L),
         limitSell( obj, L),
         marketBuy(obj),
         marketSell(obj),
         cancelBuy(obj,L),
         cancelSell(obj, L)
         )

}


bookShape<- function(obj, band){
  a=length(obj$Price[obj$Price<=bid(obj)])
  b=length(obj$Price[obj$Price<=ask(obj)])
  mid= floor((a+b)/2.0)
  
  obj[(mid-band):(mid+band),]
  
}

bookPlot<-function (obj, band){
  a=length(obj$Price[obj$Price<=bid(obj)])
  b=length(obj$Price[obj$Price<=ask(obj)])
  mid= floor((a+b)/2.0)
  
  q=c(obj$buySize[(mid-band):mid],obj$sellSize[(mid+1):(mid+band)] )
  price=obj$Price[mid+(-band:band)]
  plot(price, q, col="red",type="l",xlab="Price",ylab="Quantity")
  
}

b1<- book(1000,F)
alpha<-1
mu<-10
delta<-0.2
L<-100
for(i in 1:1000) {b1<-event(b1, L, mu, alpha, delta)}
numEvent<-100000
avgBook<- b1$Book
avgBook$buySize<-avgBook$buySize/numEvent
avgBook$sellSize<- avgBook$sellSize/numEvent
for (i in 2:numEvent){
  b1<-event(b1, L, mu, alpha, delta)
  avgBook$buySize<- avgBook$buySize+ b1$Book$buySize/numEvent
  avgBook$sellSize<- avgBook$sellSize+ b1$Book$sellSize/numEvent
}
#bookShape(avgBook, L)
bookPlot(avgBook, L)



