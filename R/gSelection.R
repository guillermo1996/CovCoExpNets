#' Calculates the best seed for GLMNET algorithm (DEPRECATED)
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param useRMSE boolean parameter
#'
#' @return a number representing the best seed
#' @export
bestSeed <- function(data, covariate, useRMSE = FALSE){
  .Deprecated("geneSelection", msg = "This algorthm is deprecated.\nUse 'geneSelection()' instead.\n Not yet supported to use list of expression matrices.")
  seed = 0
  counter = 1;
  seeds = numeric(10)
  errorsTest = numeric(10)
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

    if(useRMSE){
      predict.cvfit.test.10 <- stats::predict(cvfit.10, newx = data.test.10, type = "response", s = "lambda.min")
      errorsTest[counter] <- MLmetrics::RMSE(predict.cvfit.test.10, covariate.test.10)
    }else{
      coefic.10 <- as.matrix(stats::coef(cvfit.10,s="lambda.min"))
      seeds[counter] <- seed
      genesSelect[counter] <- length(which(coefic.10!=0))-1
    }
    counter = counter + 1

  }

  if(useRMSE){
    dat <- data.frame(seeds, errorsTest)
    seed <- dat$seeds[which.min(dat$errorsTest)]
  }else{
    dat <- data.frame(seeds, genesSelect)
    seed <- dat$seeds[which.min(dat$genesSelect)]
  }

  return(seed)
}


#' Reduce the predictors to a subset
#'
#' Estimates the most common subset of genes that the GLMNET algorithm considers
#' important to predict a covariate. The inputs can be either a list object
#' generated from previous SuCoNets functions or individual condition variables.
#'
#' @param data list of numeric expression matrices
#' @param covariate list of numeric vectors
#' @param t number of times to repeat the train/test split. Defaults to 10
#' @param k number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param n numeric value. Minimum ratio of gene appearance to be selected as a
#'   predictor
#' @param train.split numeric percentage of samples to take as train
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#'
#' @return list of vectors containing the genes selected for each condition. If
#'   only one data matrix is provided, it will return the reduced matrix, not a
#'   list.
#' @export
geneSelection <- function(data, covariate, t = 10, k = 10, n = 0.5, train.split = 0.8, sample.prob = c()){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    sample.prob = list(sample.prob)
  }

  genes.subset.combined = foreach(i = 1:length(data), .combine = "c") %:% when(i) %do%{
    data.condition = data[[i]]
    covariate.condition = covariate[[i]]
    sample.prob.condition = sample.prob[[i]]
    set.seed(0)
    seeds = sample(999999999, t)

    genes.subset = foreach(j = 1:t, .combine = "rbind") %:% when(i) %dopar%{
      set.seed(seeds[j])

      ind.train.t = sample(1:ncol(data.condition),
                           train.split*ncol(data.condition),
                           prob = sample.prob.condition)
      data.train.t = t(data.condition[, ind.train.t])
      covariate.train.t = covariate.condition[ind.train.t]

      cvfit.t = glmnet::cv.glmnet(data.train.t, covariate.train.t, nfolds = k, alpha = 1, family = "gaussian")
      detectGenes(cvfit.t)
    }

    genes.subset = table(genes.subset$Genes, dnn = "Genes") %>%
      as.data.frame() %>%
      mutate(Freq = Freq / t) %>%
      arrange(desc(Freq))

    genes.subset = genes.subset %>% filter(Freq >= n) %>% pull(Genes)
    genes.subset <- list(as.character(genes.subset))
    names(genes.subset) <- names(data)[i]
    genes.subset
  }

  if(length(genes.subset.combined) == 1){
    return(genes.subset.combined[[1]])
  }else{
    return(genes.subset.combined)
  }
}


#' Reduce the predictors to a subset using bootstrap
#'
#' Estimates the most common subset of genes that the GLMNET algorithm considers
#' important to predict a covariate using bootstrap. The inputs can be either a
#' list object generated from previous SuCoNets functions or individual
#' condition variables.
#'
#' @param data list of numeric expression matrices
#' @param covariate list of numeric vectors
#' @param B number of bootstrap resamples (default 20)
#' @param k number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param n numeric value. Minimum ratio of gene appearance to be selected as a
#'   predictor
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#'
#' @return list of vectors containing the genes selected for each condition. If
#'   only one data matrix is provided, it will return the reduced matrix, not a
#'   list.
#' @export
geneSelectionBootstrap <- function(data, covariate, B = 20, k = 10, n = 0.5, sample.prob = c()){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    sample.prob = list(sample.prob)
  }

  genes.subset.combined = foreach(i = 1:length(data), .combine = "c") %dopar% {
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    sample.prob.i = sample.prob[[i]]

    set.seed(0)
    seeds = sample(999999999, B)

    genes.subset = foreach(j = 1:B, .combine = "rbind") %dopar%{
      set.seed(seeds[j])
      ind.train = sample(1:ncol(data.i), ncol(data.i), replace = TRUE, prob = sample.prob.i)

      data.train.B = t(data.i[, ind.train])
      covariate.train.B = covariate.i[ind.train]

      cvfit.B = glmnet::cv.glmnet(data.train.B, covariate.train.B, nfolds = k, alpha=1, family = "gaussian")

      detectGenes(cvfit.B)
    }

    genes.subset = table(genes.subset$Genes, dnn = "Genes") %>%
      as.data.frame() %>%
      mutate(Freq = Freq / B) %>%
      arrange(desc(Freq))

    genes.subset = genes.subset %>% filter(Freq >= n) %>% pull(Genes)
    genes.subset <- list(as.character(genes.subset))
    names(genes.subset) <- names(data)[i]
    genes.subset
  }

  if(length(genes.subset.combined) == 1){
    return(genes.subset.combined[[1]])
  }else{
    return(genes.subset.combined)
  }
}


#' Run GLMNET algorithm with the selected genes
#'
#' The inputs can be either a list object generated from previous SuCoNets
#' functions or individual condition variables.
#'
#' @param data list of numeric expression matrices
#' @param covariate list of numeric vectors
#' @param genes.subset list of selected predictors for each condition
#' @param k number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param train.split numeric percentage of samples to take as train
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#'
#' @return list of the generated glmnet object. If only one data matrix is
#'   provided, it will return only one object, not a list.
#' @export
glmnetGenesSubset <- function(data, covariate, genes.subset, k = 10, train.split = 0.8, sample.prob = c()){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    genes.subset = list(genes.subset)
    sample.prob = list(sample.prob)
  }

  cvfit.combined = foreach(i = 1:length(data), .combine = "c") %:% when(i) %dopar%{
    data.condition = data[[i]]
    covariate.condition = covariate[[i]]
    genes.subset.condition = genes.subset[[i]]
    sample.prob.condition = sample.prob[[i]]

    set.seed(0)
    ind.train = sample(1:ncol(data.condition), train.split*ncol(data.condition), prob = sample.prob.condition)

    data.train = t(data.condition[genes.subset.condition, ind.train])
    covariate.train = covariate.condition[ind.train]

    cvfit = glmnet::cv.glmnet(data.train, covariate.train, kfold = k, alpha = 1, family = "gaussian")
    cvfit <- list(cvfit)
    names(cvfit) <- names(data)[i]
    cvfit
  }

  if(length(cvfit.combined) == 1){
    return(cvfit.combined[[1]])
  }else{
    return(cvfit.combined)
  }
}


#' Run GLMNET algorithm with the best seed (DEPRECATED)
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
#' This function calculates the variables that best predict a covariate using
#' glmnet algorithm. The inputs can be either a list object generated from
#' previous SuCoNets functions or individual condition variables.
#'
#' @param cvfit list of GMLNET object
#'
#' @return list of dataframe containing the variables selected by the glmnet
#'   algorithm and the coefficients associated with these variables. If only one
#'   object was supplied, the output will be one object and not a list.
#' @export
#'
detectGenes <- function(cvfit){
  if(!is(cvfit, "list")){
    cvfit = list(cvfit)
  }

  selected.genes.combined <- foreach(i = 1:length(cvfit), .combine = "c") %dopar%{
    cvfit.condition = cvfit[[i]]

    coefic = as.matrix(stats::coef(cvfit.condition, s="lambda.min"))
    non.zero.rows = rownames(coefic)[which(coefic!=0)]
    selected.genes = list(data.frame(Genes = non.zero.rows[-1], Coeficientes = coefic[which(coefic!=0)][-1]))
    names(selected.genes) <- names(cvfit)[i]
    selected.genes
  }

  if(length(selected.genes.combined) == 1){
    return(selected.genes.combined[[1]])
  }else{
    return(selected.genes.combined)
  }
}


#' Study of the stability of the selection of genes with which we are going to
#' work. (DEPRECATED)
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param selectedGenes dataframe with the genes selected as important by GLMNET
#'   algorithm
#'
#' @return Summary of gene frequency information
#' @export
#'
stabilitySelection <- function(data, covariate, selectedGenes){
  gs = as.character()
  seed = 0
  counter = 1;
  seeds = numeric(10)
  genesSelect = numeric(10)
  errorsTrain = numeric(10)
  errorsTest = numeric(10)

  for(i in 9:0){
    if(seed == 987654321)
      seed = seed +1
    else
      seed = seed*10+i
    set.seed(seed)

    ind.train.10 <- sample(1:ncol(data), 0.8*ncol(data))

    data.train.10 <-  t(data[,ind.train.10])
    data.test.10 <-  t(data[,-ind.train.10])
    covariate.train.10 <-  covariate[ind.train.10]
    covariate.test.10 <- covariate[-ind.train.10]

    cvfit.10<- glmnet::cv.glmnet(data.train.10, covariate.train.10, alpha=1, family = "gaussian")

    coefic.10 <- as.matrix(stats::coef(cvfit.10,s="lambda.min"))
    non.zero.rows.10 <- rownames(coefic.10)[which(coefic.10!=0)]
    selected.genes.10 <- data.frame(Genes = non.zero.rows.10[-1], Coeficientes = coefic.10[which(coefic.10!=0)][-1])

    predict.cvfit.test.10 <- stats::predict(cvfit.10, newx = data.test.10,type = "response",s = "lambda.min")
    predict.cvfit.train.10 <- stats::predict(cvfit.10, newx = data.train.10,type = "response",s = "lambda.min")

    gs = c(gs, as.character(selected.genes.10[,1]))
    seeds[counter] <- seed
    genesSelect[counter] <- nrow(selected.genes.10)
    errorsTrain[counter] <- MLmetrics::RMSE(predict.cvfit.train.10, covariate.train.10)
    errorsTest[counter] <- MLmetrics::RMSE(predict.cvfit.test.10, covariate.test.10)
    counter = counter + 1

  }

  gs = as.factor(gs)
  gs.our.genes <- data.frame(table(gs))
  id <- selectRows(as.character(gs.our.genes[,1]), as.character(selectedGenes[,1]))
  gs.our.genes <- gs.our.genes[id,]
  # as.data.frame(gs.our.genes[order(gs.our.genes[,2], decreasing = T),])
  #
  #
  # df.gene.occurrences <- data.frame(num.gene.occurrences = gs.our.genes[,2])
  return(gs.our.genes)

}


#' Co-expression network with fixed cluster sizes
#'
#' Calculate a coexpression network using the genes selected by the glmnet
#' algorithm as the best predictors of a covariate. The inputs can be either a
#' list object generated from previous SuCoNets functions or individual
#' condition variables.
#'
#' @param data expression matrix
#' @param selected.genes list of dataframes with the genes selected as important
#'   by GLMNET algorithm
#' @param tam number: indicates the number of genes in each cluster
#'
#' @return A dataframe where for each gene selected by glmnet appears the tam
#'   genes most correlated with it and that correlation. If only one object was
#'   supplied, the output will be one object and not a list.
#' @export
coexpressionNetworkFixed <- function(data, selected.genes, tam){
  if(!is(data, "list")){
    data = list(data)
    selected.genes = list(selected.genes)
  }

  r = foreach(i = 1:length(data), .combine="c") %dopar% {
    data.i = data[[i]]
    selected.genes.i = selected.genes[[i]]

    df = calculateCorrelation(data.i, selected.genes.i[1, 1], tam)

    for(j in 2:nrow(selected.genes.i)){
      df = rbind(df, calculateCorrelation(data.i, selected.genes.i[j, 1], tam))
    }

    df <- list(df)
    names(df) <- names(data)[i]
    df
  }

  if(length(data) == 1){
    return(r[[1]])
  }else{
    return(r)
  }
}

#' Co-expression network with variable cluster sizes
#'
#' Calculate a coexpression network using the genes selected by the glmnet
#' algorithm as the best predictors of a covariate. The inputs can be either a
#' list object generated from previous SuCoNets functions or individual
#' condition variables.
#'
#' @param data list of expression matrices
#' @param covariate list of numeric vectors
#' @param selected.genes list of dataframes with the genes selected as important
#'   by GLMNET algorithm
#'
#' @return A dataframe where for each gene selected by glmnet appears the tam
#'   genes most correlated with it and a vector with cluster sizes. If only one
#'   object was supplied, the output will be one object and not a list.
#' @export
#'
coexpressionNetworkVariable <- function(data, covariate, selected.genes){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    selected.genes = list(selected.genes)
  }

  r = foreach(i = 1:length(data), .combine="c") %dopar% {
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    selected.genes.i = selected.genes[[i]]

    df = calculateClusters(data.i, selected.genes.i[1, 1], covariate.i)
    tam = nrow(df)

    for(j in 2:nrow(selected.genes.i)){
      df_g = calculateClusters(data.i, selected.genes.i[j, 1], covariate.i)
      tam = c(tam, nrow(df_g))
      df = rbind(df, df_g)
    }

    df = list(list(df, tam))
    names(df) = names(data)[i]
    df
  }

  if(length(data) == 1){
    return(r[[1]])
  }else{
    return(r)
  }
}
