#' Comparison between real values of the covariate and predicted values
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param m mean of the original covariate
#' @param d standard deviation of the original covariate
#' @param cvfit object generated with GLMNET algorithm
#' @param seed number
#'
#' @return plot
#' @export
comparisonActualPredictedCovariate <- function(data, covariate, m, d, cvfit, seed, genes.subset){
  set.seed(seed)

  ind.train <- sample(1:ncol(data), 0.8*ncol(data))

  data.train <-  t(data[genes.subset, ind.train])
  covariate.train <-  covariate[ind.train]

  data.test <- t(data[as.character(genes.subset),-ind.train])
  covariate.test <- covariate[-ind.train]

  predict.cvfit.test <- stats::predict(cvfit, newx = data.test,type = "response",s = "lambda.min")
  predict.cvfit.train <- stats::predict(cvfit, newx = data.train,type = "response",s = "lambda.min")


  covariate.real.test <- covariate.test*d+m
  covariate.predicha.test <- predict.cvfit.test*d+m

  covariate.real.train <- covariate.train*d+m
  covariate.predicha.train <- predict.cvfit.train*d+m

  dats <- data.frame(x = c(covariate.real.test, covariate.real.train), y = c(covariate.predicha.test, covariate.predicha.train), type = as.factor(c(rep("Test", ncol(data)-length(ind.train)), rep("Train", length(ind.train)))))

  ggplot2::ggplot(dats, ggplot2::aes(x=dats[,1],y=dats[,2], color=type)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method='lm',se = FALSE) +
    ggplot2::xlab("Real covariate")+
    ggplot2::ylab("Predicted covariate")+
    ggplot2::labs(title="Comparison between real values of the covariate and predicted values")


}

#' Distribution of individuals according to the covariate
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param selectedGenes dataframe with the genes selected as important by GLMNET algorithm
#' @param m mean of the original covariate
#' @param d standard deviation of the original covariate
#'
#' @return plot
#' @export
distributionIndividualsCovariate <- function(data, covariate, selectedGenes, m, d){
  e <- min(covariate*d+m):max(covariate*d+m)
  dfGraf<- data.frame(e)
  dfGraf$col <- grDevices::colorRampPalette(c("blue", "red"))(length(e))

  ind <- selectRows(rownames(data),selectedGenes[,1])
  data.genes <-data[ind,]
  pca.result <- stats::prcomp(t(data.genes))

  pcas <- as.data.frame(pca.result$x,stringsAsFactors=F)
  pcas <- cbind(covariate.genes = as.numeric(covariate*d+m), pcas)
  pcas$col <- dfGraf$col[pcas$covariate.genes-19]

  graphics::plot(x = pcas$PC1, y = pcas$PC2, col=pcas$col, pch = 20, xlab = "First Principal Component", ylab = "Second Principal Component", panel.first = graphics::grid())

}

#' Histogram of gene frequencies
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


#' Bar Plot of the RMSE values for each tissue
#'
#' @param tissues list of strings containing the tissue names
#' @param dataset_folder path to the dataset folder
#' @param models_folder path to the models folder
#' @param B number of bootstrap resamples (default 20)
#'
#' @return plot
#' @export
boot.RMSE <- function(tissues, dataset_folder, models_folder, B = 20){
  tissues_short = gsub(" ", "", tissues, fixed = TRUE)
  len_tissue = length(tissues)

  df_rmse <- foreach(i = 1:len_tissue, .combine = "rbind") %:% when(i) %dopar%{
    bt = tissues_short[i]
    cat("(", i, "/", len_bt ,") Comenzamos con ", bt, "\n", sep="")

    data <- readRDS(paste0(datasets_folder, "data_", bt, ".rds"))
    age <- readRDS(paste0(datasets_folder, "age_", bt, ".rds"))
    m <- age[1]
    d <- age[2]
    age <- age[-c(1,2)]
    #cvfit <- readRDS(paste0(models_folder, "cvfit_", bt, ".rds"))
    genes.subset <- readRDS(paste0(models_folder, "genes.subset_", bt, ".rds"))

    rmse.bt.train <- numeric(B)
    rmse.bt.test <- numeric(B)
    foreach(j = 1:B) %do%{
      ind.train <- sample(1:ncol(data), ncol(data), replace = TRUE)

      data.train <- t(data[genes.subset ,ind.train])
      covariate.train <- age[ind.train]
      data.test <- t(data[genes.subset, -unique(ind.train)])
      covariate.test <- age[-unique(ind.train)]

      cvfit <- glmnet::cv.glmnet(data.train, covariate.train, alpha=1, family = "gaussian")
      selected.genes <- detectGenes("", "", cvfit)

      pred.train <- predict(cvfit, s = "lambda.min", newx = data.train, type = "response")
      rmse.bt.train[j] <- get.RMSE(covariate.train, pred.train)

      pred.test <- predict(cvfit, s = "lambda.min", newx = data.test, type = "response")
      rmse.bt.test[j] <- get.RMSE(covariate.test, pred.test)
    }

    data.frame(Tissue = str_split(tissues[i], " - ", simplify = TRUE)[, 2],
               "Train RMSE" = mean(rmse.bt.train),
               "Train std" = sd(rmse.bt.train),
               "Test RMSE" = mean(rmse.bt.test),
               "Test std" = sd(rmse.bt.test))
  }

  df_rmse = df_rmse %>% arrange(desc(Tissue))
  rmse_plot <- ggplot(data = df_rmse, aes(x = as.character(Tissue), y = Test.RMSE)) +
    geom_bar(stat="identity", position=position_dodge(), fill="steelblue", color="black") +
    geom_errorbar(aes(ymin=Test.RMSE-Test.std, ymax=Test.RMSE+Test.std), width=.2,
                  position=position_dodge(.9)) +
    theme(axis.text.x = element_text(face = "bold", angle = -75, vjust = 0.5, hjust=0)) +
    theme(axis.title.y = element_text(angle = 90, vjust = 1.5, hjust=1)) +
    theme(axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 9/20) +
    theme(axis.title.x=element_blank()) +
    labs(y = "RMSE") +
    ggtitle(bquote("RMSE por tejido")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))

  rmse_plot
}


#' Bar Plot of the RMSE values for each tissue
#'
#' @param tissues list of strings containing the tissue names
#' @param dataset_folder path to the dataset folder
#' @param models_folder path to the models folder
#' @param B number of bootstrap resamples (default 20)
#'
#' @return plot
#' @export
boot.r2_adj <- function(tissues, dataset_folder, models_folder, B = 20){
  tissues_short = gsub(" ", "", tissues, fixed = TRUE)
  len_tissue = length(tissues)

  df_r2_adj <- foreach(i = 1:len_tissue, .combine = "rbind") %:% when(i) %dopar%{
    bt = tissues_short[i]
    cat("(", i, "/", len_bt ,") Comenzamos con ", bt, "\n", sep="")

    data <- readRDS(paste0(datasets_folder, "data_", bt, ".rds"))
    age <- readRDS(paste0(datasets_folder, "age_", bt, ".rds"))
    m <- age[1]
    d <- age[2]
    age <- age[-c(1,2)]
    cvfit <- readRDS(paste0(models_folder, "cvfit_", bt, ".rds"))
    genes.subset <- readRDS(paste0(models_folder, "genes.subset_", bt, ".rds"))

    r2_adj.bt.train <- numeric(B)
    r2_adj.bt.test <- numeric(B)
    foreach(j = 1:B) %do%{
      ind.train <- sample(1:ncol(data), ncol(data), replace = TRUE)

      data.train <- t(data[genes.subset, ind.train])
      covariate.train <- age[ind.train]
      data.test <- t(data[genes.subset, -unique(ind.train)])
      covariate.test <- age[-unique(ind.train)]

      cvfit <- glmnet::cv.glmnet(data.train, covariate.train, alpha=1, family = "gaussian")
      selected.genes <- detectGenes("", "", cvfit)

      pred.train <- predict(cvfit, s = "lambda.min", newx = data.train, type = "response")
      r2_adj.bt.train[j] <- get.r2_adj(covariate.train, pred.train, length(selected.genes[, 1]))

      pred.test <- predict(cvfit, s = "lambda.min", newx = data.test, type = "response")
      r2_adj.bt.test[j] <- get.r2_adj(covariate.test, pred.test, length(selected.genes[, 1]))
    }

    data.frame(Tissue = str_split(tissues[i], " - ", simplify = TRUE)[, 2],
               "Train r2_adj" = mean(r2_adj.bt.train),
               "Train std" = sd(r2_adj.bt.train),
               "Test r2_adj" = mean(r2_adj.bt.test),
               "Test std" = sd(r2_adj.bt.test))
  }

  r2_adj_plot <- ggplot(data = df_r2_adj, aes(x = as.character(Tissue), y = Train.r2_adj)) +
    geom_bar(stat="identity", position=position_dodge(), fill="steelblue", color="black") +
    geom_errorbar(aes(ymin=Train.r2_adj-Train.std, ymax=Train.r2_adj+Train.std), width=.2,
                  position=position_dodge(.9)) +
    coord_cartesian(ylim=c(0, 1)) +
    theme(axis.text.x = element_text(face = "bold", angle = -75, vjust = 0.5, hjust=0)) +
    theme(axis.title.y = element_text(angle = 90, vjust = 1.5, hjust=1)) +
    theme(axis.text.y = element_text(face = "bold")) +
    theme(aspect.ratio = 9/20) +
    theme(axis.title.x=element_blank()) +
    labs(y = "R^2 adjusted") +
    ggtitle(bquote("R^2 adjusted por tejido en Train set")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))

  r2_adj_plot
}


#' Heatmap of the coincidences in genes selected for each tissue
#'
#' @param tissues list of strings containing the tissue names
#' @param dataset_folder path to the dataset folder
#' @param models_folder path to the models folder
#'
#' @return plot
#' @export
plotCoincidences <- function(tissues, dataset_folder, models_folder){
  tissues_short = gsub(" ", "", tissues, fixed = TRUE)
  heat_map = matrix(0, nrow=length(brain_tissues), ncol=length((brain_tissues)))

  r =
    foreach(i = 1:len_bt) %:%
      foreach(j = i:len_bt) %do% {
        tissue_i = brain_tissues_short[i]
        tissue_j = brain_tissues_short[j]

        selected.genes.i = readRDS(paste0(models_folder, "selected.genes_", tissue_i, ".rds"))
        selected.genes.j = readRDS(paste0(models_folder, "selected.genes_", tissue_j, ".rds"))

        common.genes = intersect(selected.genes.i[, 1], selected.genes.j[, 1])
        heat_map[i, j] = length(common.genes)
        heat_map[j, i] = length(common.genes)
        cat("Common between", tissue_i, "and", tissue_j, ":", length(common.genes), "\n")
      }

  plot_data <- expand.grid(X=str_replace(brain_tissues, "Brain - ", ""), Y=str_replace(brain_tissues, "Brain - ", ""))
  plot_data$Coincidencias <- c(heat_map)
  diag(heat_map) <- NA
  plot_data$Color_Coincidencias <- c(heat_map)

  heat_plot <- ggplot(plot_data, aes(x = X, y = Y, fill=Color_Coincidencias)) +
    geom_tile(color = "black", lwd = 0.2, linetype = 1) +
    geom_text(aes(label = round(Coincidencias, 1))) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits=rev) +
    theme(axis.text.x = element_text(face = "bold", angle = 75, vjust = 0.5, hjust=0)) +
    theme(axis.text.y = element_text(face = "bold", angle = 15, vjust = 0.5, hjust=1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    scale_fill_continuous(name = "Coincidencias") +
    theme(aspect.ratio = 1)

  heat_plot
}

## Bar plot of RMSE for all tissues
##
## @param tissue list of characters containing the names of the tissue to look for files
## @param model_folder path to the
##
## @return plot
## @export
# barPlotRMSE <- function(tissues, dataset_folder, models_folder){
#   tissues_short = gsub(" ", "", tissues, fixed = TRUE)
#   len_tissue = length(tissues)
#   rmse.train <- numeric(len_tissue)
#   rmse.test <- numeric(len_tissue)
#
#   r <- foreach(i = 1:len_tissue) %:% when(i) %do%{
#     bt = tissues_short[i]
#     cat("(", i, "/", len_bt ,") Comenzamos con ", bt, "\n", sep="")
#
#     data <- readRDS(paste0(datasets_folder, "data_", bt, ".rds"))
#     age <- readRDS(paste0(datasets_folder, "age_", bt, ".rds"))
#     m <- age[1]
#     d <- age[2]
#     age <- age[-c(1,2)]
#     cvfit <- readRDS(paste0(models_folder, "cvfit_", bt, ".rds"))
#     selected.genes <- readRDS(paste0(models_folder, "selected.genes_", bt, ".rds"))
#     all.genes <- readRDS(paste0(models_folder, "all.genes_", bt, ".rds"))
#
#     genes.subset = all.genes %>% filter(Freq >= 5) %>% select(Genes)
#
#     set.seed(0)
#     data.subset <- data[row.names(data) %in% genes.subset$Genes, ]
#     ind.train <- sample(1:ncol(data.subset), 0.8*ncol(data.subset))
#
#     data.train <- t(data.subset[,ind.train])
#     age.train <- age[ind.train]
#     data.test <- t(data.subset[, -ind.train])
#     age.test <- age[-ind.train]
#
#     pred.train <- predict(cvfit, s = "lambda.min", newx = data.train, type = "response")
#     pred.test <- predict(cvfit, s = "lambda.min", newx = data.test, type = "response")
#
#     rmse.train[i] <- get.RMSE(age.train, pred.train)
#     rmse.test[i] <- get.RMSE(age.test, pred.test)
#   }
#
#   rmse <- data.frame(cbind(str_split(brain_tissues, " - ", simplify = TRUE)[, 2], rmse.train, rmse.test))
#   colnames(rmse) <- c("Tissue", "Train", "Test")
#   rmse_long <- reshape2::melt(rmse, id.vars = c("Tissue"), variable.name = "Dataset", value.name = "rmse")
#   rmse_long$rmse <- as.numeric(pmax(rmse_long$rmse, 0))
#
#   rmse_plot <- ggplot(data = rmse_long, aes(x = Tissue, y = rmse, fill = Dataset)) +
#     geom_bar(stat="identity", position=position_dodge()) +
#     theme(axis.text.x = element_text(face = "bold", angle = -75, vjust = 0.5, hjust=0)) +
#     theme(axis.title.y = element_text(angle = 90, vjust = 1.5, hjust=1)) +
#     theme(axis.text.y = element_text(face = "bold")) +
#     theme(aspect.ratio = 9/20) +
#     theme(axis.title.x=element_blank()) +
#     ggtitle(bquote("RMSE por tejido y dataset")) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
#
#   rmse_plot
# }

