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
comparisonActualPredictedCovariate <- function(data, covariate, m, d, cvfit, seed){
  set.seed(seed)

  ind.train <- sample(1:ncol(data), 0.8*ncol(data))

  data.train <-  t(data[,ind.train])
  covariate.train <-  covariate[ind.train]

  data.test <- t(data[,-ind.train])
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
#' @param data expression matrix
#' @param covariate numeric vector
#' @param selectedGenes dataframe with the genes selected as important by GLMNET algorithm
#'
#' @return Summary of gene frequency information
#' @export
histGeneFreq <- function(data, covariate, selectedGenes){
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
  as.data.frame(gs.our.genes[order(gs.our.genes[,2], decreasing = T),])


  t <- data.frame(num.gene.occurrences = gs.our.genes[,2])

  ggplot2::ggplot(data=t, ggplot2::aes(num.gene.occurrences)) +
          ggplot2::geom_histogram(col="light blue",fill="light blue",alpha=0.75, breaks=seq(0, 10, by=1)) +
          ggplot2::labs(title="Histogram gene frequency", x = "Num of times that the genes are repeated", y="Num of genes")

  return(summary(gs.our.genes[,2]))

}


