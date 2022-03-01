#' Comparison between real values of the covariate and predicted values
#'
#' The inputs can be either a list object generated from previous SuCoNets
#' functions or individual condition variables.
#'
#' @param data list of expression matrices
#' @param covariate list of numeric vectors
#' @param m list of means of the original covariates
#' @param d list of standard deviations of the original covariates
#' @param cvfit list of objects generated with GLMNET algorithm
#' @param genes.subset list of selected predictors for each condition
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#'
#' @return plot
#' @export
#'
#' @examples
#' comparisonActualPredictedCovariate(data, covariate, m, d, cvfit, genes.subset)
#' comparisonActualPredictedCovariate(data[[1]], covariate[[1]], m[[1]], d[[1]], cvfit[[1]], genes.subset[[1]])
comparisonActualPredictedCovariate <- function(data, covariate, m, d, cvfit, genes.subset, sample.prob = c()){
  if(!is(data, "list")){
    data <- list(data)
    covariate <- list(covariate)
    m <- list(m)
    d <- list(d)
    cvfit <- list(cvfit)
    genes.subset <- list(genes.subset)
    sample.prob <- list(sample.prob)
  }

  r = foreach(i = 1:length(data)) %:% when(i) %do%{
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    m.i = m[[i]]
    d.i = d[[i]]
    cvfit.i = cvfit[[i]]
    genes.subset.i = genes.subset[[i]]
    sample.prob.i = sample.prob[[i]]
    set.seed(0)

    ind.train = sample(1:ncol(data.i), 0.8*ncol(data.i), prob = sample.prob.i)

    data.train = t(data.i[genes.subset.i, ind.train])
    covariate.train = covariate.i[ind.train]
    data.test = t(data.i[genes.subset.i, -ind.train])
    covariate.test = covariate.i[-ind.train]

    predict.cvfit.train = stats::predict(cvfit.i, newx = data.train, type = "response", s = "lambda.min")
    predict.cvfit.test = stats::predict(cvfit.i, newx = data.test, type = "response", s = "lambda.min")

    covariate.real.train = covariate.train*d.i + m.i
    covariate.real.test = covariate.test*d.i + m.i
    covariate.predicha.train = predict.cvfit.train*d.i + m.i
    covariate.predicha.test = predict.cvfit.test*d.i + m.i

    dats = data.frame(x = c(covariate.real.test, covariate.real.train),
                      y = c(covariate.predicha.test, covariate.predicha.train),
                      type = as.factor(c(rep("Test", ncol(data.i)-length(ind.train)), rep("Train", length(ind.train)))))

    p = ggplot2::ggplot(dats, ggplot2::aes(x=dats[,1],y=dats[,2], color=type)) +
      #ggplot2::geom_abline(ggplot2::aes(fill = "black"), intercept = 0, slope = 1) +
      ggplot2::geom_abline(ggplot2::aes(slope=1, intercept=0), show.legend = F) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method='lm', se = FALSE) +
      ggplot2::xlab("Real covariate") +
      ggplot2::ylab("Predicted covariate") +
      ggplot2::labs(title="Comparison between real values of the covariate and predicted values")
    if(length(data) != 1) p = p +
      ggplot2::labs(title=paste0("Comparison real covariate vs predicted - ", names(data)[i]))

    print(p)
  };
}


#' Distribution of individuals according to the covariate
#'
#' The inputs can be either a list object generated from previous SuCoNets
#' functions or individual condition variables
#'
#' @param data list of expression matrices
#' @param covariate list of numeric vectors
#' @param selected.genes list of dataframes with the genes selected as important
#'   by GLMNET algorithm
#' @param m list of meas of the original covariates
#' @param d list of standard deviations of the original covariates
#'
#' @return plot
#' @export
#'
#' @examples
#' distributionIndividualsCovariate(data, covariate, selected.genes, m, d)
#' distributionIndividualsCovariate(data[[1]], covariate[[1]], selected.genes[[1]], m[[1]], d[[1]])
distributionIndividualsCovariate <- function(data, covariate, selected.genes, m, d){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    selected.genes = list(selected.genes)
    m = list(m)
    d = list(d)
  }

  r = foreach(i = 1:length(data)) %do% {
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    selected.genes.i = selected.genes[[i]]
    m.i = m[[i]]
    d.i = d[[i]]

    e = min(covariate.i*d.i + m.i):max(covariate.i*d.i + m.i)
    dfGraf = data.frame(e)
    dfGraf$col = grDevices::colorRampPalette(c("blue", "red"))(length(e))

    ind = selectRows(rownames(data.i), selected.genes.i[, 1])
    data.genes = data.i[ind, ]
    pca.result = stats::prcomp(t(data.genes))

    pcas = as.data.frame(pca.result$x, stringsAsFactors = F)
    pcas = cbind(covariate.genes = as.numeric(covariate.i*d.i + m.i), pcas)
    pcas$col = dfGraf$col[pcas$covariate.genes - 19]

    graphics::plot(x = pcas$PC1,
                   y = pcas$PC2,
                   col = pcas$col,
                   pch = 20,
                   xlab = "First Principal Component",
                   ylab = "Second Principal Component",
                   panel.first = graphics::grid())
  };
}

#' Histogram of gene frequencies (DEPRECATED)
#'
#' @param gs.our.genes data.frame
#'
#' @return plot
#' @export
histGeneFreq <- function(gs.our.genes){
  df.gene.occurrences <- data.frame(num.gene.occurrences = gs.our.genes[,2])
  ggplot2::ggplot(data=df.gene.occurrences, ggplot2::aes(num.gene.occurrences)) +
    ggplot2::geom_histogram(col="light blue",fill="light blue",alpha=0.75, breaks=seq(0, 10, by=1)) +
    ggplot2::labs(title="Histogram gene frequency", x = "Num of times that the genes are repeated", y="Num of genes")
}


#' Bar Plot of the RMSE values for each condition
#'
#' The inputs can be either a list object generated from previous SuCoNets
#' functions or individual condition variables.
#'
#' @param data list of expression matrices
#' @param covariate list of numeric vectors
#' @param genes.subset list of selected predictors for each condition
#' @param B number of bootstrap resamples. Defaults to 20
#' @param sample.prob list of numeric vectors as weights when applying
#'   \link[base]{sample}
#' @param k number of kfolds to execute in the GLMNET algorithm
#'   (\link[glmnet]{cv.glmnet}). Defaults to 10
#' @param do.boxplot boolean to do a boxplot or a barplot. Defaults to False
#'
#' @return ggplot2 object
#' @export
plotRMSE <- function(data, covariate, genes.subset, m = 0, d = 1, B = 20, sample.prob = c(), k = 10, do.boxplot = FALSE){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    genes.subset= list(genes.subset)
    m = list(m)
    d = list(d)
  }else{
    if(!is(m, "list") || !is(d, "list")){
      m = rep(list(0), length(data))
      d = rep(list(1), length(data))
    }
  }

  df.rmse = foreach(i = 1:length(data), .combine = "rbind") %dopar% {
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    genes.subset.i = genes.subset[[i]]
    sample.prob.i = sample.prob[[i]]
    m.i = m[[i]]
    d.i = d[[i]]

    set.seed(0)
    seeds = sample(999999999, B)
    df.rmse.i <- foreach(j = 1:B, .combine = "rbind") %do%{
      set.seed(seeds[j])
      ind.train.B = sample(1:ncol(data.i), ncol(data.i), replace = TRUE, prob = sample.prob.i)

      data.train.B = t(data.i[genes.subset.i, ind.train.B])
      covariate.train.B = covariate.i[ind.train.B]

      data.test.B = t(data.i[genes.subset.i, -ind.train.B])
      covariate.test.B = covariate.i[-ind.train.B]

      cvfit.B = glmnet::cv.glmnet(data.train.B, covariate.train.B, alpha=1, kfold=k, family = "gaussian")

      pred.train.B = predict(cvfit.B, s = "lambda.min", newx = data.train.B, type = "response")
      pred.test.B = predict(cvfit.B, s = "lambda.min", newx = data.test.B, type = "response")

      # if(!missing(m) & !missing(d)){
      #   pred.train.B = pred.train.B*d.i + m.i
      #   pred.test.B = pred.test.B*d.i + m.i
      #   covariate.train.B = covariate.train.B*d.i + m.i
      #   covariate.test.B = covariate.train.B*d.i + m.i
      rmse.train.B = MLmetrics::RMSE(pred.train.B*d.i + m.i, covariate.train.B*d.i + m.i)
      rmse.test.B = MLmetrics::RMSE(pred.test.B*d.i + m.i, covariate.test.B*d.i + m.i)

      if(is.null(names(data)[i])) cond.name = "Condition" else cond.name = names(data)[i]

      data.frame(Condition = cond.name, rmse.train = rmse.train.B, rmse.test = rmse.test.B)
    }

    df.rmse.i
  }

  if(!do.boxplot){
    df.rmse.group = df.rmse %>%
      dplyr::group_by(Condition) %>%
      dplyr::summarize(mean.rmse.train = mean(rmse.train),
                       sd.rmse.train = sd(rmse.train),
                       mean.rmse.test = mean(rmse.test),
                       sd.rmse.test = sd(rmse.test))

    rmse.plot = ggplot2::ggplot(df.rmse.group, ggplot2::aes(x = Condition, y = mean.rmse.test)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", color = "black") +
      #ggplot2::geom_bar(data = df.rmse.group, ggplot2::aes(x = Condition, y = mean.rmse.train)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean.rmse.test - sd.rmse.test,
                                          ymax = mean.rmse.test + sd.rmse.test),
                             width = .2,
                             position = ggplot2::position_dodge(.9)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  limits = c(0, 1.07*(max(df.rmse.group$mean.rmse.test) +
                                                        max(df.rmse.group$sd.rmse.test))))
  }else{
    rmse.plot = ggplot2::ggplot(df.rmse, ggplot2::aes(x = Condition, y = rmse.test)) +
      ggplot2::geom_boxplot(fill = "steelblue")
      #ggplot2::stat_summary(fun = mean, geom="point", size=2, color="red") +
  }

  rmse.plot = rmse.plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", angle = -75, vjust = 0.5, hjust = 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = .5)) +
    ggplot2::theme(aspect.ratio = 9/20) +
    ggplot2::labs(y = "RMSE") +
    ggplot2::ggtitle(bquote("RMSE by condition for test sample"))

  if(length(data) == 1) rmse.plot = rmse.plot + ggplot2::theme(aspect.ratio = 20/9)
  return(rmse.plot)
}


#' Bar Plot of the RMSE values for each tissue (DEPRECATED)
#'
#' @param tissues list of strings containing the tissue names
#' @param dataset_folder path to the dataset folder
#' @param models_folder path to the models folder
#' @param B number of bootstrap resamples (default 20)
#'
#' @return plot
#' @export
plotR2adj <- function(data, covariate, genes.subset, B = 20, sample.prob = c(), k = 10, plot.split = "test", do.boxplot = FALSE){
  if(!is(data, "list")){
    data = list(data)
    covariate = list(covariate)
    genes.subset= list(genes.subset)
  }

  df.r2_adj = foreach(i = 1:length(data), .combine = "rbind") %dopar% {
    data.i = data[[i]]
    covariate.i = covariate[[i]]
    genes.subset.i = genes.subset[[i]]
    sample.prob.i = sample.prob[[i]]

    set.seed(0)
    seeds = sample(999999999, B)
    df.r2_adj.i <- foreach(j = 1:B, .combine = "rbind") %do%{
      set.seed(seeds[j])
      ind.train.B = sample(1:ncol(data.i), ncol(data.i), replace = TRUE, prob = sample.prob.i)

      data.train.B = t(data.i[genes.subset.i, ind.train.B])
      covariate.train.B = covariate.i[ind.train.B]

      data.test.B = t(data.i[genes.subset.i, -ind.train.B])
      covariate.test.B = covariate.i[-ind.train.B]

      cvfit.B = glmnet::cv.glmnet(data.train.B, covariate.train.B, alpha=1, kfold=k, family = "gaussian")
      selected.genes.B <- detectGenes(cvfit.B)

      pred.train.B = predict(cvfit.B, s = "lambda.min", newx = data.train.B, type = "response")
      pred.test.B = predict(cvfit.B, s = "lambda.min", newx = data.test.B, type = "response")

      r2_adj.train.B = get.r2_adj(pred.train.B, covariate.train.B, length(selected.genes[, 1]))
      r2_adj.test.B = get.r2_adj(pred.test.B, covariate.test.B, length(selected.genes[, 1]))

      if(is.null(names(data)[i])) cond.name = "Condition" else cond.name = names(data)[i]

      data.frame(Condition = cond.name, r2_adj.train = r2_adj.train.B, r2_adj.test = r2_adj.test.B)
    }

    df.r2_adj.i
  }

  if(!do.boxplot){
    df.r2_adj.group = df.rmse %>%
      dplyr::group_by(Condition)

    if(plot.split == "test"){
      df.r2_adj.group = df.r2_adj.group %>%
        dplyr::summarize(mean.r2_adj = mean(r2_adj.test),
                         sd.r2_adj = sd(r2_adj.test))
    }else{
      df.r2_adj.group = df.r2_adj.group %>%
        dplyr::summarize(mean.r2_adj = mean(r2_adj.train),
                         sd.r2_adj = sd(r2_adj.train))
    }

    r2_adj.plot = ggplot2::ggplot(df.r2_adj.group, ggplot2::aes(x = Condition, y = mean.r2_adj.test)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean.r2_adj.test - sd.r2_adj.test,
                                          ymax = mean.r2_adj.test + sd.r2_adj.test),
                             width = .2,
                             position = ggplot2::position_dodge(.9)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  limits = c(0, 1.07*(max(df.rmse.group$mean.rmse.test) +
                                                        max(df.rmse.group$sd.rmse.test))))
  }else{
    rmse.plot = ggplot2::ggplot(df.rmse, ggplot2::aes(x = Condition, y = rmse.test)) +
      ggplot2::geom_boxplot(fill = "steelblue")
      #ggplot2::stat_summary(fun = mean, geom="point", size=2, color="red") +
  }

  rmse.plot = rmse.plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", angle = -75, vjust = 0.5, hjust = 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold")) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = .5)) +
    ggplot2::theme(aspect.ratio = 9/20) +
    ggplot2::labs(y = "RMSE") +
    ggplot2::ggtitle(bquote("RMSE by condition"))

  if(length(data) == 1) rmse.plot = rmse.plot + ggplot2::theme(aspect.ratio = 20/9)
  return(rmse.plot)
}
# boot.r2_adj <- function(tissues, dataset_folder, models_folder, B = 20){
#   tissues_short = gsub(" ", "", tissues, fixed = TRUE)
#   len_tissue = length(tissues)
#
#   df_r2_adj <- foreach(i = 1:len_tissue, .combine = "rbind") %:% when(i) %dopar%{
#     bt = tissues_short[i]
#     cat("(", i, "/", len_bt ,") Comenzamos con ", bt, "\n", sep="")
#
#     data <- readRDS(paste0(datasets_folder, "data_", bt, ".rds"))
#     age <- readRDS(paste0(datasets_folder, "age_", bt, ".rds"))
#     m <- age[1]
#     d <- age[2]
#     age <- age[-c(1,2)]
#     cvfit <- readRDS(paste0(models_folder, "cvfit_", bt, ".rds"))
#     genes.subset <- readRDS(paste0(models_folder, "genes.subset_", bt, ".rds"))
#
#     r2_adj.bt.train <- numeric(B)
#     r2_adj.bt.test <- numeric(B)
#     foreach(j = 1:B) %do%{
#       ind.train <- sample(1:ncol(data), ncol(data), replace = TRUE)
#
#       data.train <- t(data[genes.subset, ind.train])
#       covariate.train <- age[ind.train]
#       data.test <- t(data[genes.subset, -unique(ind.train)])
#       covariate.test <- age[-unique(ind.train)]
#
#       cvfit <- glmnet::cv.glmnet(data.train, covariate.train, alpha=1, family = "gaussian")
#       selected.genes <- detectGenes(cvfit)
#
#       pred.train <- glmnet::predict.glmnet(cvfit, s = "lambda.min", newx = data.train, type = "response")
#       r2_adj.bt.train[j] <- get.r2_adj(covariate.train, pred.train, length(selected.genes[, 1]))
#
#       pred.test <- glmnet::predict.glmnet(cvfit, s = "lambda.min", newx = data.test, type = "response")
#       r2_adj.bt.test[j] <- get.r2_adj(covariate.test, pred.test, length(selected.genes[, 1]))
#     }
#
#     data.frame(Tissue = str_split(tissues[i], " - ", simplify = TRUE)[, 2],
#                "Train r2_adj" = mean(r2_adj.bt.train),
#                "Train std" = sd(r2_adj.bt.train),
#                "Test r2_adj" = mean(r2_adj.bt.test),
#                "Test std" = sd(r2_adj.bt.test))
#   }
#
#   r2_adj_plot <- ggplot(data = df_r2_adj, aes(x = as.character(Tissue), y = Train.r2_adj)) +
#     geom_bar(stat="identity", position=position_dodge(), fill="steelblue", color="black") +
#     geom_errorbar(aes(ymin=Train.r2_adj-Train.std, ymax=Train.r2_adj+Train.std), width=.2,
#                   position=position_dodge(.9)) +
#     coord_cartesian(ylim=c(0, 1)) +
#     theme(axis.text.x = element_text(face = "bold", angle = -75, vjust = 0.5, hjust=0)) +
#     theme(axis.title.y = element_text(angle = 90, vjust = 1.5, hjust=1)) +
#     theme(axis.text.y = element_text(face = "bold")) +
#     theme(aspect.ratio = 9/20) +
#     theme(axis.title.x=element_blank()) +
#     labs(y = "R^2 adjusted") +
#     ggtitle(bquote("R^2 adjusted por tejido en Train set")) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
#
#   r2_adj_plot
# }


#' Heatmap of the coincidences in genes selected for each condition
#'
#' The inputs can be either a list object generated from previous SuCoNets
#' functions or individual condition variables.
#'
#' @param selected.genes list of dataframes with the genes selected as important
#'   by GLMNET algorithm
#'
#' @return ggplot2 object
#' @export
plotCoincidences <- function(selected.genes){
  if(!is(selected.genes, "list")){
    print("This function needs the selected genes for different conditions. Only genes for one was supplied")
    return()
  }

  heat_map = matrix(0, nrow=length(selected.genes), ncol=length(selected.genes))
  foreach(i = 1:length(selected.genes)) %:%
    foreach(j = i:length(selected.genes)) %do% {
      common.genes = length(intersect(selected.genes[[i]][, 1], selected.genes[[j]][, 1]))
      heat_map[i, j] = common.genes
      heat_map[j, i] = common.genes
  };

  plot.coin = expand.grid(X = names(selected.genes), Y = names(selected.genes))
  plot.coin$coin = c(heat_map)
  diag(heat_map) = NA
  plot.coin$color_coin = c(heat_map)

  heat.plot <- ggplot2::ggplot(plot.coin, ggplot2::aes(x = X, y = Y, fill = color_coin)) +
    ggplot2::geom_tile(color = "black", lwd = 0.2, linetype = 1)+
    ggplot2::geom_text(ggplot2::aes(label = round(coin, 1))) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", angle = 75, vjust = 0.5, hjust = 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold", angle = 15, vjust = 0.5, hjust = 1)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_fill_continuous(name = "Coincidencias") +
    ggplot2::theme(aspect.ratio = 1)

  return(heat.plot)
}
