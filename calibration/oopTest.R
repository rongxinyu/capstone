rm(list=ls())
coords <-function(a, b){ 
  coords <- list(x = a, y = b) 
  class(coords) <- "coords"
  return (coords)
}  
p <- coords(103,6)
print<-function(obj){
  sprintf("(%.2f, %.2f)",obj$x, obj$y )
}
print(p)