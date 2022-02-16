#' Calculate correlation
#'
#' Given a coexpression matrix x and a gene (which is represented as a numerical vector) y from that matrix,
#' this function returns a subset of x with the genes most correlated with y.
#' The size of this subset is tam
#'
#' @param x coexpression matrix
#' @param y factor
#' @param tam number
#'
#' @return data.frame
#' @export
calculateCorrelation <- function(x,y, tam){

  corr <- stats::cor(x=t(x), y = t(x[which(rownames(x)==y),]))
  ind <- order(abs(corr), decreasing = T)[1:(tam+1)]
  df <- data.frame(main.gen = y, genes = utils::head(rownames(corr)[ind], (tam+1)), cor = corr[ind])

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
jaccardIndex <- function(df, tam){
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


#' Selects the indexes of the rows in x that match the rows in y
#'
#' @param x char vector
#' @param y char vector
#'
#' @return numeric vector
#' @export
selectRows <- function(x,y){
  indexes = data.frame()
  indexes = which(x==y[1])

  for(n in 2:length(y)){
    indexes[n] = which(x==y[n])
  }
  return(indexes)
}



#' Calculates clusters by varying the size
#'
#' @param x coexpression matrix
#' @param y factor
#' @param covariate numeric vector
#'
#' @return dataframe with clusters
#' @export
calculateClusters <- function(x,y, covariate){

  corr <- stats::cor(x=t(x), y = t(x[which(rownames(x)==y),]))
  ind <- order(abs(corr), decreasing = T)
  df <- data.frame(gen.principal = y, genes = rownames(corr)[ind], cor = corr[ind])

  tam <- 2
  indx <- selectRows(rownames(x), df[1:tam,2])
  mydata <- data.frame(covariate = covariate, t(x[indx,]))
  mymodel = stats::lm(covariate~ .,data=mydata)
  adjusted_r2 <- as.numeric(summary(mymodel)[9])
  r2 <- adjusted_r2
  diff <- 1

  while (diff > 10^(-5)){
    tam <- tam+1
    indx <- selectRows(rownames(x), df[1:tam,2])
    mydata <- data.frame(covariate = covariate, t(x[indx,]))
    mymodel = stats::lm(covariate~ .,data=mydata)
    ad_r2 <- as.numeric(summary(mymodel)[9])

    diff <- ad_r2 -adjusted_r2

    adjusted_r2 <- ad_r2
    r2 <- c(r2, adjusted_r2)
  }

  r2=c(0,r2)
  if(diff > 0)
    return(cbind(df[1:tam,], Adjusted.R2 = r2))
  else
    return(cbind(df[1:(tam-1),], Adjusted.R2 = utils::head(r2,n=(tam-1))))
}


#' Running the function gprofiler
#'
#' @param selectedGenes dataframe with the genes selected as important by GLMNET algorithm
#' @param tam numeric vector containing the number of genes forming each cluster
#' @param data coexpression matrix
#' @param df dataframe with clusters
#'
#' @return The information obtained by running gprofiler
#' @export
#'
running.gprofiler <- function(selectedGenes, tam, data, df){
  all.genes = list()
  all.genes[[selectedGenes[1,1]]] = df[1:tam[1], 2]
  for (i in 2:nrow(selectedGenes)) {
    all.genes[[selectedGenes[i,1]]] = df[(1+sum(tam[1:(i-1)])):(sum(tam[1:i])), 2]
  }

  background <- rownames(data)

  output.gprofiler2 <- gprofiler2::gost(all.genes,
                                        correction_method="fdr",
                                        custom_bg = background,
                                        sources = c("GO","KEGG","REAC"),
                                        domain_scope = "custom",
                                        organism = "hsapiens",
                                        exclude_iea = F)
  return(output.gprofiler2)
}


#' Getting the r^2 value for a prediction
#'
#' @param real.values numeric vector of real values of the covariate
#' @param pred.values numeric vector of predicted values of the covariate
#'
#' @return The calculated r^2
#' @export
get.r2 <- function(real.values, pred.values){
  RSS = sum((pred.values - real.values)^2)
  TSS = sum((pred.values - mean(real.values))^2)
  return(1 - RSS/TSS)
}


#' Getting the adjusted r^2 value for a prediction
#'
#' @param real.values numeric vector of real values of the covariate
#' @param pred.values numeric vector of predicted values of the covariate
#' @param n number of samples
#' @param k number of predictors (length of selected genes)
#'
#' @return The calculated adjuted r^2
#' @export
get.r2_adj <- function(real.values, pred.values, k){
  RSS = sum((pred.values - real.values)^2)
  TSS = sum((pred.values - mean(real.values))^2)

  #coefic <- as.matrix(stats::coef(cvfit, s="lambda.min"))
  #k = length(which(coefic!=0)) - 1
  n = length(pred.values)

  r2_adj = 1 - (RSS/(n-k-1))/(TSS/(n-1))
  return(r2_adj)
}


#' Getting the RMSE for a model and a dataset
#'
#' @param real.values numeric vector of real values of the covariate
#' @param pred.values numeric vector of predicted values of the covariate
#'
#' @return The calculated adjuted RMSE
#' @export
get.RMSE <- function(real.values, pred.values){
  RSS = sum((pred.values - real.values)^2)
  n = length(pred.values)
  RMSE = sqrt(RSS/n)
  return(RMSE)
}


