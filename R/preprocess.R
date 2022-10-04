#' Normalize a vector
#'
#' Normalizes a vector by subtracting its mean and dividing by its standard
#' deviation
#'
#' @param x numeric vector or list of numeric vectors
#'
#' @return numeric vector or list of numeric vectors
#' @export
#'
#' @examples
#' normalize(1:4)
normalize <- function(x){
  return.list = is(x, "list")
  if(!return.list) x = list(x)

  x.combined = foreach(i = 1:length(x), .combine = "c") %do%{
    x.i = x[[i]]
    m = mean(x.i)
    d = stats::sd(x.i)
    x.i = (x.i - m)/d

    x.i = list(list(mean = m, standard.deviation = d, covariate = x.i))
    names(x.i) = if(return.list) names(x)[i] else "Condition"
    x.i
  }

  return(returnList(return.list, x.combined))
}


#' Log2 + 1 to all matrix
#'
#' Applies a factor log2 + 1 to all values in the data matrix.
#'
#' @param data non-negative numerical matrix
#'
#' @return numerical matrix
#' @export
applyLog2 <- function(data){
  logger::log_debug("Applying log2(TPM + 1)")
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    data.i = log2(data.i + 1)

    data.i = list(data.i)
    names(data.i) = if(return.list) names(data)[i] else "Condition"
    data.i
  }

  return(returnList(return.list, data.combined))
}



#' Scaling, centering and normalizing a matrix
#'
#' Scales a non-negative numerical matrix by centering and normalizing its data.
#'
#' @param data non-negative numerical matrix with column names
#' @param by.columns boolean. Whether to normalize by columns of rows, by default FALSE
#'
#' @return numerical matrix
#' @export
applyZScore <- function(data, by.columns = FALSE){
  logger::log_debug(paste0("Applying ZScore by ", ifelse(by.columns, "columns", "rows")))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    if(!by.columns) data.i = t(data.i)
    preProc = caret:: preProcess(data.i, method = c("center", "scale"))
    data.i = stats::predict(preProc, data.i)

    if(!by.columns) data.i = t(data.i)
    data.i = list(data.i)
    names(data.i) = if(return.list) names(data)[i] else "Condition"
    data.i
  }

  return(returnList(return.list, data.combined))
}



#' Remove near zero variance predictors
#'
#' Identifies and eliminates predictors that have very few unique values and the
#' ratio of the frecuency of the most common value to the frequency of the
#' second most common value is large. More info at \link[caret]{nearZeroVar}.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param freqCut numerical, the cutoff for the ratio of the most common value to the second most common value
#' @param uniqueCut numerical, the cutoff for the percentage of distinct values out of the number of total samples
#' @param sigDigits numerical, number of significant digits to round
#'
#' @return reduced numerical matrix
#' @export
applyNearZeroVar <- function(data, freqCut = 20, uniqueCut = 5,  sigDigits = 1){
  logger::log_debug(paste0("Applying nearZeroVar. Rounding for ",
                           sigDigits,
                           " significant numbers. Parameters: freqCut = ",
                           freqCut,
                           " \\ uniqueCut = ",
                           uniqueCut))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    r.data.i = if(sigDigits == 0) data.i else apply(data.i, 1, signif, digits = 1)

    remove.i = caret::nearZeroVar(r.data.i, freqCut = 20, uniqueCut = 5)

    logger::log_debug(paste0("\tFor condition ", condition.name, ", removed ", length(remove.i), " predictors (", nrow(data.i), " -> ", nrow(data.i)-length(remove.i), ")."))
    if(length(remove.i > 0)){
      data.i <- data.i[-remove.i, ]
    }

    data.i = list(data.i)
    names(data.i) = if(return.list) names(data)[i] else "Condition"
    data.i
  }

  return(returnList(return.list, data.combined))
}


#' Remove predictors following an activation threshold
#'
#' Identifies and eliminates predictors that do not fulfill the activation
#' threshold requirement. By default, a predictor should have a value greater
#' that 0.1 in at least the 80% of the samples.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param threshold numeric: minimum TPM value to consider a predictor 'active'
#' @param perc numeric: minimum percentage of activated predictors in all
#'   samples.
#'
#' @return reduced numerical matrix
#' @export
applyActiveThreshold <- function(data, threshold = 0.1, perc = 0.8){
  logger::log_debug(paste0("Applying active threshold. Requirement: ", threshold, " activation in at least ", perc*100, "%  of samples."))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    remove.i = apply(data.i, 1, function(x) length(x[x >= threshold])/length(x) <= perc)
    remove.i = which(unname(remove.i))

    condition.name = if(return.list) names(data)[i] else "condition"
    logger::log_debug(paste0("\tFor condition ", condition.name, ", removed ", length(remove.i), " predictors (", nrow(data.i), " -> ", nrow(data.i)-length(remove.i), ")."))

    if(length(remove.i > 0)){
      data.i <- data.i[-remove.i, ]
    }

    data.i = list(data.i)
    if(return.list) names(data.i) = names(data)[i]
    data.i
  }
  return(returnList(return.list, data.combined))
}


#' Corrects GTEx naming convention
#'
#' Removes the version term of the ENSEMBL naming convention of GTEx data.
#' More information [here](https://www.ensembl.org/Help/Faq?id=488)
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
correctGTExNames <- function(data){
  logger::log_debug(paste0("Correcting GTEx gene names"))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]

    names <- stringr::str_split(row.names(data.i), "[.]", simplify = TRUE)[, 1]
    row.names(data.i) <- names

    data.i <- list(data.i)
    names(data.i) <- names(data)[i]
    data.i
  }

  return(returnList(return.list, data.combined))
}


#' Remove non protein-coding genes
#'
#' Identifies and eliminates genes which are not protein coding.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#'
#' @return reduced numerical matrix
#' @export
applyProteinCodingRequirement <- function(data){
  library(biomaRt)

  mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

  logger::log_debug(paste0('Applying "protein" coding requirement'))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %do%{
    data.i = data[[i]]
    condition.name = if(return.list) names(data)[i] else "condition"
    gene.list = rownames(data.i)

    genesNames <- biomaRt::getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"),
                                 filters="ensembl_gene_id",
                                 values=unique(gene.list),
                                 mart=mart,
                                 useCache=FALSE)
    genesNames <- genesNames[match(unique(gene.list), genesNames$ensembl_gene_id), ]
    genesNames <- genesNames[which(genesNames$gene_biotype %in% "protein_coding"), ]
    new.data = data.i[genesNames$ensembl_gene_id, ]
    rownames(new.data) = genesNames$external_gene_name
    if(!all(data.frame(ensemmbl_gene_id = rownames(data.i[genesNames$ensembl_gene_id, ]), external_gene_name = rownames(new.data)) == genesNames[, c("ensembl_gene_id", "external_gene_name")])){
      logger::log_error(paste0("Not a valid conversion was found for condition ", condition.name))
    }

    if(sum(rownames(new.data) == "") > 0){
      logger::log_warn(paste0("\tFound ", sum(rownames(new.data) == ""), " empty gene names for condition ", condition.name), ".")
      new.data = new.data[rownames(new.data) != "", ]
    }

    remove.i = nrow(data.i) - nrow(new.data)
    logger::log_debug(paste0("\tFor condition ", condition.name, ", removed ", remove.i, " predictors (", nrow(data.i), " -> ", nrow(new.data), ")."))

    new.data = list(new.data)
    if(return.list) names(new.data) = names(data)[i]
    new.data
  }

  return(returnList(return.list, data.combined))
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
  logger::log_debug(paste0("Applying correlation removal. Parameters: cutoff = ", cutoff, ". "))
  return.list = is(data, "list")
  if(!return.list) data = list(data)

  data.combined = foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    condition.name = if(return.list) names(data)[i] else "condition"
    correlation.i = stats::cor(t(data.i))
    remove.i = caret::findCorrelation(correlation.i, cutoff = cutoff)

    logger::log_debug(paste0("\tFor condition ", condition.name, ", removed ", length(remove.i), " predictors (", nrow(data.i), " -> ", nrow(data.i)-length(remove.i), ")."))
    if(length(remove.i > 0)){
      data.i <- data.i[-remove.i, ]
    }

    data.i = list(data.i)
    if(return.list) names(data.i) = names(data)[i]
    rm(correlation.i)
    gc()
    data.i
  }

  return(returnList(return.list, data.combined))
}


#' Splits data matrix by condition
#'
#' Creates a lists of numeric matrices. Each matrix corresponds to a column in
#' the 'bool.matrix'. The names of the conditions can be provided or
#' extracted automatically from the 'bool.matrix' column names.
#'
#' @param data numerical matrix with predictors as rows and samples as columns
#' @param bool.matrix boolean matrix with samples as rows and conditions as
#'   columns. The data matrix will be splitted according to this boolean matrix
#' @param names names assigned to each different numerical matrix in the
#'   returned list. If not provided, will be selected from the 'bool.matrix'
#'
#' @return named list of numerical matrices
#' @export
splitByCondition <- function(data, bool.matrix, names = NA){
  if(is.na(names)) names <- colnames(bool.matrix)
  len.names = length(names)
  data.combined = foreach(i = 1:len.names) %dopar%{
    condition = names[i]
    data[, rownames(bool.matrix[bool.matrix[, condition] == 1, ])]
  }
  names(data.combined) <- names
  data.combined
}


#' Preprocess the input data matrices
#'
#' Applies all the preprocessing required for the CovCoExpNets pipeline to work.
#' Each individual process can be deactivadted (but not recommended).
#'
#' If input data is from the GTEx foundation, it is requiered to execute the naming
#' correction to process the protein coding genes.
#'
#' @param data numerical matrix (or list of numerical matrices) with predictors
#'   as rows and samples as columns
#' @param doLog2 boolean to apply \link[CovCoExpNets]{applyLog2}
#' @param doAT boolean to apply \link[CovCoExpNets]{applyActiveThreshold}
#' @param doNZV boolean to apply \link[CovCoExpNets]{applyNearZeroVar}
#' @param doGTExCorr boolean to apply \link[CovCoExpNets]{correctGTExNames}
#' @param doProtCod boolean to apply \link[CovCoExpNets]{applyProteinCodingRequirement}
#' @param doZSco boolean to apply \link[CovCoExpNets]{applyZScore}
#' @param doCorr boolean to apply \link[CovCoExpNets]{applyCorrSkew}
#'
#' @param threshold numeric, required for \link[CovCoExpNets]{applyActiveThreshold}. minimum TPM value to consider a predictor 'active'.
#'   Defaults to 0.1
#' @param perc numeric, required for \link[CovCoExpNets]{applyActiveThreshold}. minimum percentage of activated predictors in all
#'   samples. Defaults to 0.8
#' @param freqCut numerical, required for \link[CovCoExpNets]{applyNearZeroVar}. The cutoff for the ratio of the most common value to the second most common value
#' @param uniqueCut numerical, required for \link[CovCoExpNets]{applyNearZeroVar}. The cutoff for the percentage of distinct values out of the number of total samples
#' @param sigDigits numerical, required for \link[CovCoExpNets]{applyNearZeroVar}. Number of significant digits to round
#' @param by.columns boolean. Whether to normalize by columns of rows, by default FALSE
#' @param cutoff numeric: value for the pair-wise absolute correlation cutoff (\link[caret]{findCorrelation}). Defaults to 0.9
#'
#' @return named list of reduced numerical matrices. If only one data matrix is
#'   provided, it will return the reduced matrix, not a list.
#' @export
dataPreprocess <- function(data, doLog2 = T, doAT = T, doNZV = T, doGTExCorr = T, doProtCod = T, doZSco = T, doCorr = F,
                           threshold = 0.1, perc = 0.8, freqCut = 20, uniqueCut = 5,  sigDigits = 1, by.columns = F, cutoff = 0.9){
  if(doLog2) data <- applyLog2(data)
  if(doAT) data <- applyActiveThreshold(data, threshold = threshold, perc = perc)
  if(doNZV) data <- applyNearZeroVar(data, freqCut = freqCut, uniqueCut = uniqueCut, sigDigits = sigDigits)
  if(doGTExCorr) data <- correctGTExNames(data)
  if(doProtCod) data <- applyProteinCodingRequirement(data)
  if(doCorr) data <- applyCorrSkew(data, cutoff = cutoff)
  if(doZSco) data <- applyZScore(data, by.columns)

  return(data)
}
