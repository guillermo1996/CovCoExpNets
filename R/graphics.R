#' Comparison between real values of the covariate and predicted values
#'
#' @param data expression matrix
#' @param covariate numeric vector
#' @param m mean of the original covariate
#' @param d standard deviation of the original covariate
#' @param cvfit object generated with GLMNET algorithm
#' @param semilla number
#'
#' @return plot
#' @export
comparisonActualPredictedCovariate <- function(data, covariate, m, d, cvfit, semilla){
  set.seed(semilla)

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

  datos <- data.frame(x = c(covariate.real.test, covariate.real.train), y = c(covariate.predicha.test, covariate.predicha.train), type = as.factor(c(rep("Test", ncol(data)-length(ind.train)), rep("Train", length(ind.train)))))

  ggplot2::ggplot(datos, ggplot2::aes(x=datos[,1],y=datos[,2], color=type)) +
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

  ind <- seleccionarFilas(rownames(data),selectedGenes[,1])
  datos.genes <-data[ind,]
  pca.result <- stats::prcomp(t(datos.genes))

  pcas <- as.data.frame(pca.result$x,stringsAsFactors=F)
  pcas <- cbind(covariate.genes = as.numeric(covariate*d+m), pcas)
  pcas$col <- dfGraf$col[pcas$covariate.genes-19]

  graphics::plot(x = pcas$PC1, y = pcas$PC2, col=pcas$col, pch = 20, xlab = "First Principal Component", ylab = "Second Principal Component", panel.first = graphics::grid())

}
