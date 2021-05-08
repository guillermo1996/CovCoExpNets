#' Normalize a vector
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
  return(x)
}

#' Scaling, centering and normalizing a matrix
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
#' scn(matrix(1:4, nr= 2, nc = 2))
scn <- function(data){
  data <- log2(data + 1)
  preProc <- caret:: preProcess(data, method = c("center", "scale"))
  data <- predict(preProc, data)
  return(data)
}
