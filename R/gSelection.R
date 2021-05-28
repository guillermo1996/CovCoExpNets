#' Calculates the best seed for GLMNET algorithm
#'
#' @param data expression matrix
#' @param covariate numeric vector
#'
#' @return a number representing the best seed
#' @export
bestSeed <- function(data,covariate){
  seed = 0
  counter = 1;
  seeds = numeric(10)
  genesSelect = numeric(10)

  for(i in 9:0){
    if(seed == 987654321)
      seed = seed +1
    else
      seed = seed*10+i
    set.seed(seed)

    ind.train.10 <- sample(1:ncol(data), 0.8*ncol(data))
    data.train.10 <-  t(data[,ind.train.10])
    covariate.train.10 <-  covariate[ind.train.10]

    cvfit.10<- glmnet::cv.glmnet(data.train.10, covariate.train.10, alpha=1, family = "gaussian")

    coefic.10 <- as.matrix(stats::coef(cvfit.10,s="lambda.min"))
    seeds[counter] <- seed
    genesSelect[counter] <- length(which(coefic.10!=0))-1
    counter = counter + 1

  }

  dat <- data.frame(seeds, genesSelect)

  seed <- dat$seeds[which.min(dat$genesSelect)]
  return(seed)
}


#' Run GLMNET algorithm with the best seed
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param seed number
#'
#' @return The generated glmnet object
#' @export
glmnetGenes <- function(data,covariate, seed){
  set.seed(seed)

  ind.train <- sample(1:ncol(data), 0.8*ncol(data))

  data.train <-  t(data[,ind.train])
  covariate.train <-  covariate[ind.train]

  data.test <- t(data[,-ind.train])
  covariate.test <- covariate[-ind.train]
  cvfit<- glmnet::cv.glmnet(data.train, covariate.train, alpha=1, family = "gaussian")
  return(cvfit)

}


#' Detect main genes
#'
#' Given an expression matrix and a covariate, this function calculates the variables
#' that best predict that covariate using glmnet algorithm
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param cvfit GMLNET object
#'
#' @return Dataframe containing the variables selected by the glmnet algorithm and the coefficients associated with these variables
#' @export
#'
detectGenes <- function(data, covariate, cvfit){

  coefic <- as.matrix(stats::coef(cvfit,s="lambda.min"))
  non.zero.rows <- rownames(coefic)[which(coefic!=0)]
  selected.genes <- data.frame(Genes = non.zero.rows[-1], Coeficientes = coefic[which(coefic!=0)][-1])

  return(selected.genes)
}

#' Co-expression network with fixed cluster sizes
#'
#' Calculate a coexpression network using the genes selected by the glmnet algorithm
#' as the best predictors of a covariate.
#'
#' @param data expression matrix
#' @param genes dataframe with the genes selected by the glmnet algorithm and their coefficients
#' @param tam number: indicates the number of genes in each cluster
#'
#' @return A dataframe where for each gene selected by glmnet appears the tam genes most correlated with it and that correlation
#' @export
coexpressionNetworkFixed <- function(data, genes, tam){
  df <- calculateCorrelation(data, genes[1,1], tam)
  for (i in 2:nrow(genes)) {
    df <- rbind(df, calculateCorrelation(data, genes[i,1], tam))
  }
  return(df)
}

#' Co-expression network with variable cluster sizes
#'
#' @param data expression matrix
#' @param selectedGenes dataframe with the genes selected by the glmnet algorithm and their coefficients
#' @param covariate numeric vector
#' @param seed number
#'
#' @return A dataframe where for each gene selected by glmnet appears the tam genes most correlated with it and a vector with cluster sizes
#' @export
#'
coexpressionNetworkVariable <- function(data, selectedGenes, covariate, seed){

  df <- calculateClusters(data, selectedGenes[1,1], selectedGenes[,1], covariate, seed)
  tam <- nrow(df)
  for (i in 2:nrow(selectedGenes)) {
    df_g <- calculateClusters(data, selectedGenes[i,1], selectedGenes[,1], covariate, seed)
    tam <- c(tam, nrow(df_g))
    df <- rbind(df, df_g)
  }
  return(list(df, tam))

}
