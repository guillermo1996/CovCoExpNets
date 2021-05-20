#' Normalize a vector
#'
#' Normalizes a vector by subtracting its mean and dividing by its standard deviation
#'
#' @param x numeric vector
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' normalize(1:4)
normalize <- function(x){
  m <- mean(x)
  d <- stats::sd(x)
  x <- (x-m)/d
  return(c(mean = m, standard.deviation = d, x))
  #return(x)
}

#' Scaling, centering and normalizing a matrix
#'
#' Scales a non-negative numerical matrix by applying the logarithm function
#' and centers and normalizes its data. When applying the logarithm function
#' it adds one to avoid calculating the logarithm of zero.
#'
#' @param data non-negative numerical matrix with column names
#'
#' @return numerical matrix
#' @export
#'
#' @examples
#' m <- matrix(1:4, nr= 2, nc = 2)
#' colnames(m)<-c("1", "2")
#' scn(m)
scn <- function(data){
  data <- log2(data + 1)
  preProc <- caret:: preProcess(data, method = c("center", "scale"))
  data <- stats::predict(preProc, data)
  return(data)
}

#' Remove redundant predictors
#'
#' Reduces the number of predictors by eliminating those with near-zero variance
#' and from pairs that have high correlation (>= 90%) eliminate one of the two predictors.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#'
#' @return numerical matrix
#' @export
#'
#' @examples
#' nycfl = nycflights13::flights
#' # We won´t use those non numerical predictors
#' nycfl = na.omit(nycfl[,-c(1,10,12,13,14,19)])
#' ncol(nycfl)
#' nycfl2 <- t(rRedundantPredictors(t(nycfl)))
#' ncol(nycfl2)
rRedundantPredictors <- function(data){
  r.data <- apply(data, 1, round)
  ind <-caret::nearZeroVar(r.data,freqCut = 20, uniqueCut = 5)
  if(length(ind) > 0)
    data <- data[-ind,]
  correlacion <- stats::cor(t(data))
  eliminar <- caret::findCorrelation(correlacion, cutoff=0.9)
  if(length(eliminar > 0))
    data <- data[-eliminar,]
  return(data)
}



