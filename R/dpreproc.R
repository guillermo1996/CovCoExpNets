#' Normalize a vector
#'
#' Normalizes a vector by subtracting its mean and dividing by its standard
#' deviation
#'
#' @param x numeric vector or list of numeric vectors.
#'
#' @return numeric vector or list of numeric vectors.
#' @export
#'
#' @examples
#' normalize(1:4)
normalize <- function(x){
  if(!is(x, "list")) x <- list(x)
  x.combined <- foreach(i = 1:length(x), .combine = "c") %do%{
    x.condition <- x[[i]]
    m <- mean(x.condition)
    d <- stats::sd(x.condition)
    x.condition <- (x.condition-m)/d
    x.condition <- list(c(mean = m, standard.deviation = d, x.condition))
    names(x.condition) <- names(x)[i]
    x.condition
  }

  #return(c(mean = m, standard.deviation = d, x))
  if(length(x.combined) == 1){
    return(x.combined[[1]])
  }else{
    return(x.combined)
  }
}

#' Scaling, centering and normalizing a matrix
#'
#' Scales a non-negative numerical matrix by applying the logarithm function and
#' centers and normalizes its data. When applying the logarithm function it adds
#' one to avoid calculating the logarithm of zero.
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

#' Remove redundant predictors (DEPRECATED)
#'
#' Reduces the number of predictors by eliminating those with near-zero variance
#' and from pairs that have high correlation (>= 90%) eliminate one of the two
#' predictors.
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
  .Deprecated("dataPreprocess()", msg = "This function is deprecated.\nUse 'dataPreprocess()' instead. Use at your own risk.")
  r.data <- apply(data, 1, round)
  ind <-caret::nearZeroVar(r.data,freqCut = 20, uniqueCut = 5)
  if(length(ind) > 0)
    data <- data[-ind,]
  correlation <- stats::cor(t(data))
  remove <- caret::findCorrelation(correlation, cutoff=0.9)
  if(length(remove > 0))
    data <- data[-remove,]
  return(data)
}


#' Remove near zero variance predictors
#'
#' Identifies and eliminates predictors that have very few unique values and the
#' ratio of the frecuency of the most common value to the frequency of the
#' second most common value is large. More info at \link[caret]{nearZeroVar}.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#'
#' @return reduced numerical matrix
#' @export
applyNearZeroVar <- function(data){
  r.data = apply(data, 1, round)
  ind = caret::nearZeroVar(r.data, freqCut = 20, uniqueCut = 5)
  if(length(ind) > 0){
    data <- data[-ind, ]
  }
  return(data)
}


#' Remove predictors following an activation threshold
#'
#' Identifies and eliminates predictors that do not fulfill the activation
#' threshold requirement. By default, a predictor should have a valur greater
#' that 0.1 in at least the 80% of the samples.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param threshold numeric: minimum TPM value (after applying
#'   \link[SuCoNets]{scn}) to consider a predictor 'active'
#' @param perc numeric: minimum percentage of activated predictors in all
#'   samples.
#'
#' @return reduced numerical matrix
#' @export
applyActiveThreshold <- function(data, threshold = 0.1, perc = 0.8){
  remove = apply(data, 1, function(x) length(x[x > threshold])/length(x) > perc)
  remove = which(unname(remove))
  if(length(remove > 0)){
    data <- data[-remove, ]
  }
  return(data)
}


#' Remove correlated predictors
#'
#' Identifies and eliminates one of two predictors that have high correlation
#' (>= 90%).
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param cutoff numeric: value for the pair-wise absolute correlation cutoff
#'   (\link[caret]{findCorrelation}).
#'
#' @return reduced numerical matrix
#' @export
applyCorrSkew <- function(data, cutoff = 0.9){
  correlation = stats::cor(t(data))
  remove = caret::findCorrelation(correlation, cutoff = cutoff)
  if(length(remove > 0)){
    data <- data[-remove, ]
  }
  return(data)
}


#' Splits data matrix
#'
#' Creates a lists of numeric matrices. Each matrix corresponds to a column in
#' the 'df.bool_data' matrix. The names of the conditions can be provoided or
#' extracted automatically from the 'df.bool_data' column names.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param df.bool_data boolean matrix with samples as rows and conditions as
#'   columns. The data matrix will be splitted according to this boolean matrix
#' @param names names assigned to each different numerical matrix in the
#'   returned list. If not provided, will be selected from the 'df.bool_data'
#'   matrix
#'
#' @return named list of numerical matrices
#' @export
dataSplit <- function(data, df.bool_data, names){
  if(missing(names)) names <- colnames(df.bool_data)
  len.names = length(names)
  data.combined = foreach(i = 1:len.names) %dopar%{
    condition = names[i]
    data[, rownames(df.bool_data[df.bool_data[, condition] == 1, ])]
  }
  names(data.combined) <- names
  data.combined
}


#' Preprocess the input data matrices
#'
#' Applies all the preprocessing required for the SuCoNets pipeline to work.
#' Each individual process can be deactivadted (but not recommended). It also
#' includes the possibility of applying a Gene name correction for GTExPortal
#' data.
#'
#' @param data numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param doSCN boolean to apply \link[SuCoNets]{scn}
#' @param doNZV boolean to apply \link[SuCoNets]{applyNearZeroVar}
#' @param doAT boolean to apply \link[SuCoNets]{applyActiveThreshold}
#' @param doCorr boolean to apply \link[SuCoNets]{applyCorrSkew}
#' @param threshold numeric: minimum TPM value (after applying
#'   \link[SuCoNets]{scn}) to consider a predictor 'active'. Defaults to 0.1
#' @param perc numeric: minimum percentage of activated predictors in all
#'   samples. Defaults to 0.8
#' @param cutoff numeric: value for the pair-wise absolute correlation cutoff
#'   (\link[caret]{findCorrelation}). Defaults to 0.9
#' @param correctGTExPortal: boolean to apply a correction to the GTExPortal
#'   gene naming convention. Defaults to True
#'
#' @return named list of reduced numerical matrices. If only one data matrix is
#'   provided, it will return the reduced matrix, not a list.
#' @export
dataPreprocess <- function(data, doSCN = TRUE, doNZV = TRUE, doAT = TRUE, doCorr = TRUE,
                           threshold = 0.1, perc = 0.8, cutoff = 0.9, correctGTExPortal = TRUE){
  if(!is(data, "list")) data <- list(data)

  data.combined <- foreach(i = 1:length(data), .combine = "c") %dopar% {
    data.i = data[[i]]

    if(doSCN) data.i <- scn(data.i)
    if(doNZV) data.i <- applyNearZeroVar(data.i)
    if(doAT) data.i <- applyActiveThreshold(data.i, threshold, perc)
    if(doCorr) data.i <- applyCorrSkew(data.i, cutoff)
    if(correctGTExPortal){
      names <- stringr::str_split(row.names(data.i), "[.]", simplify = TRUE)[, 1]
    }else names <- row.names(data.i)

    row.names(data.i) <- names

    data.i <- list(data.i)
    names(data.i) <- names(data)[i]
    data.i
  }

  if(length(data.combined) == 1){
    return(data.combined[[1]])
  }else{
    return(data.combined)
  }
}
