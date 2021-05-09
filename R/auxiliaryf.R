#' Calculate correlation
#'
#' Given a coexpression matrix x and a gene (which is represented as a numerical vector) y from that matrix,
#' this function returns a subset of x with the genes most correlated with y.
#' The size of this subset is tam
#'
#' @param x coexpression matrix
#' @param y numeric vector
#' @param tam number
#'
#' @return
#' @export
calculateCorrelation <- function(x,y, tam){

  corr <- stats::cor(x=t(x), y = t(x[which(rownames(x)==y),]))
  ind <- order(abs(corr), decreasing = T)[1:(tam+1)]
  df <- data.frame(gen.principal = y, genes = utils::head(rownames(corr)[ind], (tam+1)), cor = corr[ind])

  return(df)
}

#' Calculate the Jaccard index
#'
#' Given a dataframe df formed by sets of genes, this function calculates the Jaccard index
#' of each set with the rest and returns a summary with the information of all the calculated indexes.
#' Among the information returned by this function is the median of all the indexes.
#' The closer the median is to zero, the more independent the sets will be.
#'
#' @param df dataframe
#' @param tam size of the sets
#'
#' @return summary of Jaccard indexes
#' @export
indiceJaccard <- function(df, tam){
  indJac <- numeric(sum(1:((nrow(df)/(tam+1))-1)))
  cont <- 1

  for (i in 1:(nrow((nrow(df)/(tam+1)))-1)) {
    for (j in (i+1):(nrow((nrow(df)/(tam+1))))) {
      indJac[cont] <- length(intersect(df[(i+(i-1)*tam):(i+(i-1)*tam+tam),2], df[(j+(j-1)*tam):(j+(j-1)*tam+tam),2])) / length(union(df[(i+(i-1)*tam):(i+(i-1)*tam+tam),2], df[(j+(j-1)*tam):(j+(j-1)*tam+tam),2]))
      cont <- cont + 1
    }
  }

  return(summary(indJac))
}
