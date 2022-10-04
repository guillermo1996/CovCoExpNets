#' Getting the r^2 value for a prediction
#'
#' @param real.values numeric vector of real values of the covariate
#' @param pred.values numeric vector of predicted values of the covariate
#'
#' @return The calculated r^2
#' @export
getR2 <- function(real.values, pred.values){
  RSS = sum((pred.values - real.values)^2)
  TSS = sum((real.values - mean(real.values))^2)
  return(1 - RSS/TSS)
}

#' Getting the adjusted r^2 value for a prediction
#'
#' @param r2 numeric, r2 value of the prediction
#' @param n number of samples
#' @param k number of predictors
#'
#' @return The calculated adjuted r^2
#' @export
get.r2_adj <- function(r2, n, p){
  r2_adj = 1 - ((1-r2)*(n-1)/(n-p-1))
  return(r2_adj)
}


#' Evaluates the CovCoExpNets models
#'
#' Given the models, testing datasets and the subset of genes, this
#' function generates the main metrics RMSE, returned predictors, r2
#' and adjusted r2.
#'
#' If provided with the train dataset, it also calculates the metrics
#' for it.
#'
#' The optional parameters "m" and "d" represents the mean and standard
#' deviation of the input covariate as returned by \link[CovCoExpNets]{normalize}
#'
#' @param cvfit list of models generated with \link[CovCoExpNets]{glmnetGenesSubset}
#' @param data.test list of numeric expression matrices for testing
#' @param covariate.test list of numeric vectors for testing
#' @param genes.subset list of selected predictors for each condition
#' @param data.train list of numeric expression matrices for training
#' @param covariate.test list of numeric vectors for training
#' @param m list of means for the covariate
#' @param d list of standard deviations for the covariate
#'
#' @return The calculated adjuted r^2
#' @export
evaluateModel <- function(cvfit, data.test, covariate.test, genes.subset, data.train = NA, covariate.train = NA, m = NA, d = NA){
  return.list = is(cvfit, "list")
  if(!return.list){
    cvfit = list(cvfit)
    data.test = list(data.test)
    covariate.test = list(covariate.test)
    genes.subset = list(genes.subset)
    data.train = list(data.train)
    covariate.train = list(covariate.train)
    m = list(m)
    d = list(d)
  }

  if(missing(data.train)) data.train = rep(list(NA), length(cvfit))
  if(missing(covariate.train)) covariate.train = rep(list(NA), length(cvfit))
  if(missing(m)) m = rep(list(0), length(cvfit))
  if(missing(d)) d = rep(list(1), length(cvfit))

  foreach(i = 1:length(cvfit), .combine = "rbind") %dopar%{
    cvfit.i = cvfit[[i]]
    data.test.i = data.test[[i]]
    covariate.test.i = covariate.test[[i]]
    genes.subset.i = genes.subset[[i]]

    data.train.i = data.train[[i]]
    covariate.train.i = covariate.train[[i]]
    m.i = m[[i]]
    d.i = d[[i]]

    condition.name = if(return.list) names(cvfit)[i] else "Condition"
    relevant.genes = extractModelGenes(cvfit.i)$Genes

    print(is.na(data.train.i[[1]]))
    if(!is.na(data.train.i[[1]] & !is.na(covariate.train.i[[1]]))){
      predict.train = predict(cvfit.i, s = "lambda.min", newx = t(data.train.i[genes.subset.i, ]), type = "response")
      rmse.train = MLmetrics::RMSE(predict.train*d.i + m.i, covariate.train.i*d.i + m.i)
      r2.train = getR2(covariate.train.i, predict.train)
      r2_adj.train = get.r2_adj(r2.train, ncol(data.train.i), length(relevant.genes))
      Samples.train = ncol(data.train.i)
    }else{
      rmse.train = NA; r2.train = NA; r2_adj.train = NA; Samples.train = NA
    }

    predict.test = predict(cvfit.i, s = "lambda.min", newx = t(data.test.i[genes.subset.i, ]), type = "response")
    rmse.test = MLmetrics::RMSE(predict.test*d.i + m.i, covariate.test.i*d.i + m.i)
    r2.test = getR2(covariate.test.i, predict.test)
    r2_adj.test = get.r2_adj(r2.test, ncol(data.test.i), length(relevant.genes))

    avg.covariate = mean(covariate.train.i)*d.i + m.i
    rmse.null = MLmetrics::RMSE(avg.covariate, covariate.test.i*d.i + m.i)

    data.frame(Condition = condition.name,
               rmse.test = rmse.test,
               rmse.train = rmse.train,
               rmse.bootstrap = 0.632*rmse.test + 0.368*rmse.train,
               rmse.null = rmse.null,
               Samples.train = Samples.train,
               Samples.test = ncol(data.test.i),
               Initial.predictors = nrow(data.test.i),
               Returned.predictors = length(relevant.genes),
               r2.train = r2.train,
               r2.test = r2.test,
               r2_adj.train = r2_adj.train,
               r2_adj.test = r2_adj.test)
    # boot = i,
    # n = n.j)
  }
}
#evaluateModel(cvfit, data.test, covariate.test, genes.subset, data.train, covariate.train)


#' Generate and evaluates GLMNET models
#'
#' From a given train and test datasets, it generates the GLMNET models
#' and evaluates the most relevant metrics: RMSE, returned predictors, r2
#' and adjusted r2.
#'
#' The optional parameters "m" and "d" represents the mean and standard
#' deviation of the input covariate as returned by \link[CovCoExpNets]{normalize}
#'
#' @param data.train list of numeric expression matrices for training
#' @param covariate.test list of numeric vectors for training
#' @param data.test list of numeric expression matrices for testing
#' @param covariate.test list of numeric vectors for testing
#' @param m list of means for the covariate
#' @param d list of standard deviations for the covariate
#' @param glmnet.family see \link[glmnet]{cv.glmnet} family parameter.
#'
#' @return data.frame with the most relevant metrics for each condition given
#' @export
evaluateModelGLMNET <- function(data.train, covariate.train, data.test, covariate.test, glmnet.family = "gaussian", m = NA, d = NA){
  return.list = is(data.train, "list")
  if(!return.list){
    data.train = list(data.train)
    covariate.train = list(covariate.train)
    data.test = list(data.test)
    covariate.test = list(covariate.test)
    m = list(m)
    d = list(d)
  }

  if(missing(m)) m = rep(list(0), length(data.train))
  if(missing(d)) m = rep(list(1), length(data.train))

  foreach(i = 1:length(data.train), .combine = "rbind") %dopar%{
    data.train.i = data.train[[i]]
    covariate.train.i = covariate.train[[i]]
    data.test.i = data.test[[i]]
    covariate.test.i = covariate.test[[i]]
    m.i = m[[i]]
    d.i = d[[i]]

    condition.name = if(return.list) names(data.train)[i] else "Condition"
    cvfit.i = glmnet::cv.glmnet(t(data.train.i), covariate.train.i, alpha = 1, family = glmnet.family)
    relevant.genes = extractModelGenes(cvfit.i)$Genes

    predict.train = predict(cvfit.i, s = "lambda.min", newx = t(data.train.i), type = "response")
    rmse.train = MLmetrics::RMSE(predict.train*d.i + m.i, covariate.train.i*d.i + m.i)
    r2.train = getR2(covariate.train.i, predict.train)
    r2_adj.train = get.r2_adj(r2.train, ncol(data.train.i), length(relevant.genes))

    predict.test = predict(cvfit.i, s = "lambda.min", newx = t(data.test.i), type = "response")
    rmse.test = MLmetrics::RMSE(predict.test*d.i + m.i, covariate.test.i*d.i + m.i)
    r2.test = getR2(covariate.test.i, predict.test)
    r2_adj.test = get.r2_adj(r2.test, ncol(data.test.i), length(relevant.genes))

    avg.covariate = mean(covariate.train.i)*d.i + m.i
    rmse.null = MLmetrics::RMSE(avg.covariate, covariate.test.i*d.i + m.i)

    data.frame(Condition = condition.name,
               rmse.test = rmse.test,
               rmse.train = rmse.train,
               rmse.bootstrap = 0.632*rmse.test + 0.368*rmse.train,
               rmse.null = rmse.null,
               Samples.train = ncol(data.train.i),
               Samples.test = ncol(data.test.i),
               Initial.predictors = nrow(data.train.i),
               Returned.predictors = length(relevant.genes),
               r2.train = r2.train,
               r2.test = r2.test,
               r2_adj.train = r2_adj.train,
               r2_adj.test = r2_adj.test)
    # boot = i,
    # n = n.j)
  }
}
#evaluateModelGLMNET(data.train, covariate.train, data.test, covariate.test, m = m, d = d)


#' Summarize a metric across all repetitions
#'
#'
#' @param df numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param covariate list of numeric vectors
#' @param seed fixed seed to generate split the dataset in train/test
#'
#' @return list by folds of data.train, data.test, covariate.train and covariate.test for
#'   every condition given.
#' @export
groupMetrics <- function(df, metrics = c("rmse.test", "Returned.predictors", "r2_adj.train")){
  df.group = data.frame(Condition = unique(df$Condition))
  r = foreach(i = 1:length(metrics)) %do%{
    var.name = metrics[i]

    if(var.name %in% c("Returned.predictors")){
      digits = 0
    }else if(var.name %in% c("r2_adj.train")){
      digits = 3
    }else{
      digits = 1
    }

    sum.df = summaryVars(df, var.name, c("Condition"), na.rm = T) %>% ungroup() %>% select(!!sym(var.name), ci) %>% round(digits)
    names(sum.df) = c(paste0(var.name, ".mean"), paste0(var.name, ".ci"))

    df.group = cbind(df.group, sum.df)
  }

  return(df.group)
}


# extractCondition <- function(genes.freq, diag = "default"){
#   return.list = is(genes.freq, "list")
#   if(!return.list){
#     genes.freq = list(genes.freq)
#   }
#   heatmap.combined = foreach(i = 1:length(genes.freq), .combine = "c") %do%{
#     genes.freq.i = genes.freq[[i]]
#     condition.name = names(genes.freq)[i]
#     max.iter = max(genes.freq.i$iter)
#
#     heatmap = matrix(0, nrow = max.iter, ncol = max.iter)
#     r = foreach(j = 1:max.iter) %:% foreach(k = j:max.iter) %do%{
#       genes.j = extractGenesIter(genes.freq.i, j)
#       genes.k = extractGenesIter(genes.freq.i, k)
#
#
#       if(j == k){
#         if(is.na(diag)){
#           heatmap[j, k] = NA
#         }else if(diag == "default"){
#           heatmap[j, k] = length(genes.j)
#         }else{
#           heatmap[j, k] = diag
#         }
#       }else{
#         heatmap[j, k] = length(intersect(genes.j, genes.k))/length(genes.k)
#         heatmap[k, j] = length(intersect(genes.j, genes.k))/length(genes.j)
#       }
#     }
#
#     heatmap = list(heatmap)
#     names(heatmap) = condition.name
#     heatmap
#   }
#
#   return(returnList(return.list, heatmap.combined))
# }
