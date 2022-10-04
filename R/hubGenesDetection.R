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
extractModelGenes <- function(cvfit, genes.freq = NA){
  return.list = is(cvfit, "list")
  if(!return.list){
    cvfit = list(cvfit)
    genes.freq = list(genes.freq)
  }

  given.genes.freq = !missing(genes.freq)

  extracted.genes.combined = foreach(i = 1:length(cvfit), .combine = "c") %dopar%{
    cvfit.i = cvfit[[i]]

    coefic = as.matrix(stats::coef(cvfit.i, s="lambda.min"))
    non.zero.rows = rownames(coefic)[which(coefic!=0)]
    extracted.genes = data.frame(Genes = non.zero.rows[-1], Coeficientes = coefic[which(coefic!=0)][-1]) %>%
      arrange(desc(abs(Coeficientes)))

    if(given.genes.freq){
      genes.freq.i = genes.freq[[i]]
      genes.coefficients = genes.freq.i %>%
        group_by(Genes) %>%
        summarise(Avg_Coefficients = mean(Coeficientes), Appearances = n())

      extracted.genes = left_join(extracted.genes, genes.coefficients, by = "Genes") %>%
        arrange(desc(abs(Avg_Coefficients)))
    }

    extracted.genes = list(extracted.genes)
    names(extracted.genes) = if(return.list) names(data)[i] else "Condition"
    extracted.genes
  }

  return(returnList(return.list, extracted.genes.combined))
}


#' GLMNET repetitions
#'
#' Executes GLMNET a set number of times (t) and extracts the genes selected as
#' relevant and its coefficient. The inputs can be either a list object
#' generated from previous CovCoExpNets functions or individual condition variables.
#'
#' It can also measure the RMSE of each iteration if train.spli != 1 or if extra data
#' was provided in data.test.extra and covariate.test.extra
#'
#' @param data list of numeric expression matrices
#' @param covariate list of numeric vectors
#' @param t number of times to repeat the train/test split. Defaults to 10
#' @param k.folds number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param train.split numeric percentage of samples to take as train
#' @param iter.RMSE boolean, whether to measure the RMSE of all iterations
#' @param data.test.extra list of numeric expression matrices to calculate the RMSE
#' @param covariate.test.extra list of numeric vectors to calculate the RMSE
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#' @param seed fixed seed to generate split the dataset in train/test
#' @param glmnet.family see \link[glmnet]{cv.glmnet} family parameter.
#'
#' @return list of vectors containing the genes selected for each condition. If
#'   only one data matrix is provided, it will return the reduced matrix, not a
#'   list.
#' @export
#'
#geneFrequency(data, age, data.test.extra = sapply(data, function(x) x[, 1:100]), covariate.test.extra = sapply(age, function(x) x[1:100]), train.split = 0.5, iter.RMSE = F)
geneFrequency <- function(data, covariate, t = 10, k.folds = 10, train.split = 1,
                          iter.RMSE = F, data.test.extra, covariate.test.extra,
                          sample.prob = c(), seed = sample(1:999999, 1), glmnet.family = "gaussian"){
  return.list = is(data, "list")
  if(!return.list){
    data = list(data)
    covariate = list(covariate)
    sample.prob = list(sample.prob)
    data.test.extra = list(data.test.extra)
    covariate.test.extra = list(covariate.test.extra)
  }

  if(missing(data.test.extra)) data.test.extra = rep(list(NA), length(data))
  if(missing(covariate.test.extra)) covariate.test.extra = rep(list(NA), length(data))

  genes.freq.combined = foreach(i = 1:length(data), .combine = "c") %do%{
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    sample.prob.i = sample.prob[[i]]

    data.test.i = data.test.extra[[i]]
    covariate.test.i = covariate.test.extra[[i]]

    genes.freq = foreach(j = 1:t, .combine = "rbind") %dopar%{
      set.seed(seed + j + i)
      ind.train = sample(1:ncol(data.i), train.split*ncol(data.i), prob = sample.prob.i)

      data.train = t(data.i[, ind.train])
      covariate.train = covariate.i[ind.train]

      cvfit.t = glmnet::cv.glmnet(data.train, covariate.train, nfolds = k.folds, alpha = 1, family = glmnet.family)

      if(iter.RMSE){
        data.test = t(data.i[, -ind.train])
        covariate.test = covariate.i[-ind.train]
        if(!is.na(data.test.i[[1]]) & !is.na(covariate.test.i[[1]])){
          data.test = rbind(data.test, t(data.test.i))
          covariate.test = c(covariate.test, covariate.test.i)
        }

        predict.t = predict(cvfit.t, s = "lambda.min", newx = data.test)
        bind_cols(extractModelGenes(cvfit.t), iter = j, RMSE = MLmetrics::RMSE(predict.t, covariate.test))
      }else{
        bind_cols(extractModelGenes(cvfit.t), iter = j)
      }
    }

    genes.freq = list(genes.freq)
    names(genes.freq) = if(return.list) names(data)[i] else "Condition"
    genes.freq
  }

  return(returnList(return.list, genes.freq.combined))
}


#' Reduce the predictors to a subset
#'
#' Estimates the most common subset of genes that the GLMNET algorithm considers
#' important to predict a covariate. The input must be the output of
#' \link[CovCoExpNets]{geneFrequency}
#'
#' @param genes.freq list of vectors containing the genes selected for each condition.
#'   Shold be the output of \link[CovCoExpNets]{geneFrequency}
#' @param mrfa numeric, either the minimum frequency of appearance of a gene (0 <= mrda <= 1)
#'   or the minimum number of times that a gene must appear in the
#'   glmnet repetitions to be extracted (1 < mrfa < t)
#' @param force.mrfa boolean, whether to force the mrfa to be the given value. If set
#'   to TRUE, the algorithm will keep decreasing the mrfa if no genes are being extracted
#'
#' @return list of vectors containing the genes selected for each condition. If
#'   only one data matrix is provided, it will return the reduced matrix, not a
#'   list.
#' @export
reduceGenes <- function(genes.freq, mrfa = 0.9, force.mrfa = T, relative = T){
  return.list = is(genes.freq, "list")
  if(!return.list){
    genes.freq <- list(genes.freq)
  }

  if(mrfa == 1 & missing(relative)) logger::log_warn("Specified mrfa = 1 is ambiguous. Please specify with the parameter 'relative = TRUE' if you are in the range 0 <= mrfa <= 1, or 'relative = FALSE' if 1 <= mrfa <= t. Defaults to TRUE")

  min.mrfa = if(mrfa <= 1 & relative) 0 else 1
  red.interval = if(mrfa <= 1 & relative) -0.1 else -1

  genes.subset.combined <- foreach(i = 1:length(genes.freq), .combine = "c") %dopar%{
    genes.freq.i = genes.freq[[i]]
    max.iter = if(mrfa <= 1 & relative) max(genes.freq.i$iter) else 1

    genes.freq.i = genes.freq.i %>%
      group_by(Genes) %>%
      summarise(Freq = n() / max.iter) %>%
      ungroup() %>%
      arrange(desc(Freq))

    for(min.freq in seq(mrfa, min.mrfa, red.interval)){
      genes.freq.pruned = genes.freq.i %>%
        filter(Freq >= min.freq) %>%
        pull(Genes) %>% as.character() %>% list()

      if(force.mrfa | length(genes.freq.pruned[[1]]) > 1) break
    }

    names(genes.freq.pruned) = if(return.list) names(genes.freq)[i] else "Condition"
    genes.freq.pruned
  }

  return(returnList(return.list, genes.subset.combined))
}


#' Run GLMNET algorithm with the selected genes
#'
#' The inputs can be either a list object generated from previous CovCoExpNets
#' functions or individual condition variables. Aditionally, a list of
#' predictors must be provided to prune the training datsets.
#'
#' @param data list of numeric expression matrices
#' @param covariate list of numeric vectors
#' @param genes.subset list of selected predictors for each condition
#' @param k.folds number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param train.split numeric percentage of samples to take as train
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#' @param seed fixed seed to generate split the dataset in train/test
#'
#' @return list of the generated glmnet object. If only one data matrix is
#'   provided, it will return only one object, not a list.
#' @export
glmnetGenesSubset <- function(data, covariate, genes.subset, k.folds = 10, train.split = 1.0, sample.prob = c(), seed = sample(1:999999, 1), glmnet.family = "gaussian", return.genes.subset = F){
  return.list = is(data, "list")
  if(!return.list){
    data = list(data)
    covariate = list(covariate)
    genes.subset = list(genes.subset)
    sample.prob = list(sample.prob)
  }

  cvfit.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    genes.subset.i = genes.subset[[i]]
    sample.prob.i = sample.prob[[i]]

    set.seed(seed)
    ind.train = sample(1:ncol(data.i), train.split*ncol(data.i), prob = sample.prob.i)

    data.train = t(data.i[genes.subset.i, ind.train])
    covariate.train = covariate.i[ind.train]

    cvfit = glmnet::cv.glmnet(data.train, covariate.train, nfolds = k.folds, alpha = 1, family = glmnet.family)
    if(return.genes.subset){
      cvfit = list(list(cvfit = cvfit, genes.relevant = extractModelGenes(cvfit)))
    }else{
      cvfit = list(cvfit)
    }
    names(cvfit) <- if(return.list) names(data)[i] else "Condition"
    cvfit
  }

  return(returnList(return.list, cvfit.combined))
}
