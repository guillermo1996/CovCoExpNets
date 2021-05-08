standardizing<- function(x){
  m <- mean(x)
  d <- sd(x)
  x <- (x-m)/d
  return(x)
}
