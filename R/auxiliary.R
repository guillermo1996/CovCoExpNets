#' Return function for all list functions
#'
#' Interal function to return either the whole list or only the first element.
#' It is employed when a function can return both a single item or a list of items
#' (usually, when several conditions are processed at the same time)
#'
#' @param return.list boolean. Whether to return a list or not
#' @param data list
#'
#' @return the list or the first element of it
#' @export
returnList <- function(return.list, data){
  if(return.list){
    return(data)
  }else{
    return(data[[1]])
  }
}

#' Split data and covariate with bootstrapping
#'
#' Given a data matrix and a covariate vector (or a list of them), it
#' splits them into training set and testing set using a bootstrapping technique
#'
#' @param data numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param covariate list of numeric vectors
#' @param seed fixed seed to generate split the dataset in train/test
#'
#' @return named list with data.train, data.test, covariate.train and covariate.test
#'   for every condition given.
#' @export
splitDataBootstrap <- function(data, covariate, seed = sample(1:9999999, 1)){
  return.list = is(data, "list")
  if(!return.list){
    data = list(data)
    covariate = list(covariate)
  }

  data.split.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    covariate.i = covariate[[i]]

    set.seed(seed + i)
    ind.train = sample(1:ncol(data.i), ncol(data.i), replace = T)

    data.i.train = data.i[, ind.train]
    data.i.test = data.i[, -ind.train]
    covariate.i.train = covariate.i[ind.train]
    covariate.i.test = covariate.i[-ind.train]

    data.split = list(list(data.train = data.i.train,
                           data.test = data.i.test,
                           covariate.train = covariate.i.train,
                           covariate.test = covariate.i.test))
    names(data.split) = if(return.list) names(data)[i] else "Condition"
    data.split
  }

  data.train = list()
  data.test = list()
  covariate.train = list()
  covariate.test = list()

  for(i in 1:length(data)){
    data.i = data.split.combined[[i]]
    data.name = names(data.split.combined)[i]
    data.train[[data.name]] = data.i$data.train
    data.test[[data.name]] = data.i$data.test
    covariate.train[[data.name]] = data.i$covariate.train
    covariate.test[[data.name]] = data.i$covariate.test
  }

  data.split = list(data.train = data.train,
                    data.test = data.test,
                    covariate.train = covariate.train,
                    covariate.test = covariate.test)

  return(returnList(return.list, data.split))
}


#' Split data and covariate with cross-validation
#'
#' Given a data matrix and a covariate vector (or a list of them), it
#' splits them into training set and testing sets following a cross-validation
#' splitting process.
#'
#'
#' @param data numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param covariate list of numeric vectors
#' @param seed fixed seed to generate split the dataset in train/test
#'
#' @return list by folds of data.train, data.test, covariate.train and covariate.test for
#'   every condition given.
#' @export
splitDataCV <- function(data, covariate, k.folds, seed = sample(1:9999999, 1)){
  return.list = is(data, "list")
  if(!return.list){
    data = list(data)
    covariate = list(covariate)
  }

  data.folds = foreach(k = 1:k.folds, .combine = "c") %do%{
    data.combined = foreach(i = 1:length(data), .combine = "c") %do% {
      data.i = data[[i]]
      covariate.i = covariate[[i]]

      set.seed(seed)
      folds <- cut(seq(1, ncol(data.i)), breaks = k.folds, labels = FALSE)
      folds <- sample(folds)

      ind.train <- which(folds != k)
      data.i.train = data.i[, ind.train]
      data.i.test = data.i[, -ind.train]
      covariate.i.train = covariate.i[ind.train]
      covariate.i.test = covariate.i[-ind.train]

      data.split = list(list(data.train = data.i.train,
                             data.test = data.i.test,
                             covariate.train = covariate.i.train,
                             covariate.test = covariate.i.test))
      names(data.split) = if(return.list) names(data)[i] else "Condition"
      data.split
    }

    data.train = list()
    data.test = list()
    covariate.train = list()
    covariate.test = list()

    for(i in 1:length(data)){
      data.name = names(data.combined)[i]
      data.i = data.combined[[i]]
      if(return.list){
        data.train[[data.name]] = data.i$data.train
        data.test[[data.name]] = data.i$data.test
        covariate.train[[data.name]] = data.i$covariate.train
        covariate.test[[data.name]] = data.i$covariate.test
      }else{
        data.train = data.i$data.train
        data.test = data.i$data.test
        covariate.train = data.i$covariate.train
        covariate.test = data.i$covariate.test
      }
    }

    data.fold <- list(list(data.train = data.train,
                           data.test = data.test,
                           covariate.train = covariate.train,
                           covariate.test = covariate.test))
    names(data.fold) = k
    data.fold
  }

  return(data.folds)
}


#' Summarize a metric across all repetitions
#'
#'
#' @param data numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param covariate list of numeric vectors
#' @param seed fixed seed to generate split the dataset in train/test
#'
#' @return list by folds of data.train, data.test, covariate.train and covariate.test for
#'   every condition given.
#' @export
summaryVars <- function(df = NULL, variable, groups = NULL, conf.interval = 0.95, na.rm = T){
  customLength = function(column, na.rm = T){
    if(na.rm){
      return(sum(!is.na(column)))
    }else{
      return(length(column))
    }
  }

  func = function(df, var, na.rm){
    data.frame(counter = customLength(df[[var]], na.rm = na.rm),
               mean = mean(df[[var]], na.rm = na.rm),
               sd = sd(df[[var]], na.rm = na.rm)
    )
  }

  sum.data = df %>%
    group_by_at(groups) %>%
    do(func(., variable, na.rm))
  sum.data = rename(sum.data, !!variable := mean)
  sum.data$se <- sum.data$sd / sqrt(sum.data$counter)
  sum.data$ci <- sum.data$se*qt(conf.interval/2 + 0.5, sum.data$counter - 1)

  return(sum.data)
}


#' Extracts the relevant genes for a particular iteration
#'
#' @param gene.freq data.frame containing the GLMNET repetitions results.
#'   Should be the output of \link[CovCoExpNets]{geneFrequency}
#' @param iter numeric, particular to iteration to extract relevant
#'   genes from.
#'
#' @return the genes selected as relevant for a particular iteration.
#' @export
extractGenesIter <- function(gene.freq, iteration){
  gene.freq %>%
    filter(iter == iteration) %>%
    pull(Genes) %>%
    as.character()
}


#' Converts a list of predictors into a One Hot Encoding list
#'
#' @param all.names the total list of predictors.
#' @param selected.names the predictors that are going to be encoded.
#'
#' @return named vector of 1s and 0s
#' @export
toOneHotEncoding = function(all.names, selected.names){
  if(is(selected.names, "factor")) selected.names = as.character(selected.names)

  one.hot = rep(0, length(all.names))
  names(one.hot) = all.names
  one.hot[selected.names] = 1
  one.hot
}


#' Calculates the jaccard index of two sets
#'
#' @param a set of elements of the first vector.
#' @param b set of elements of the second vector.
#'
#' @return jaccard index between the two sets.
#' @export
jaccardIndex <- function(a, b){
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
