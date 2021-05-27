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
  semilla = 0
  contador = 1;
  semillas = numeric(10)
  genesSelect = numeric(10)
  erroresTrain = numeric(10)
  erroresTest = numeric(10)

  for(i in 9:0){
    if(semilla == 987654321)
      semilla = semilla +1
    else
      semilla = semilla*10+i
    set.seed(semilla)

    ind.train.10 <- sample(1:ncol(data), 0.8*ncol(data))

    data.train.10 <-  t(data[,ind.train.10])
    data.test.10 <-  t(data[,-ind.train.10])
    covariate.train.10 <-  covariate[ind.train.10]
    covariate.test.10 <- covariate[-ind.train.10]

    cvfit.10<- glmnet::cv.glmnet(data.train.10, covariate.train.10, alpha=1, family = "gaussian")

    coefic.10 <- as.matrix(stats::coef(cvfit.10,s="lambda.min"))
    filas.dist.cero.10 <- rownames(coefic.10)[which(coefic.10!=0)]
    genes.seleccionados.10 <- data.frame(Genes = filas.dist.cero.10[-1], Coeficientes = coefic.10[which(coefic.10!=0)][-1])

    predict.cvfit.test.10 <- stats::predict(cvfit.10, newx = data.test.10,type = "response",s = "lambda.min")
    predict.cvfit.train.10 <- stats::predict(cvfit.10, newx = data.train.10,type = "response",s = "lambda.min")

    gs = c(gs, as.character(genes.seleccionados.10[,1]))
    semillas[contador] <- semilla
    genesSelect[contador] <- nrow(genes.seleccionados.10)
    erroresTrain[contador] <- MLmetrics::RMSE(predict.cvfit.train.10, covariate.train.10)
    erroresTest[contador] <- MLmetrics::RMSE(predict.cvfit.test.10, covariate.test.10)
    contador = contador + 1

  }

  gs = as.factor(gs)
  gs.nuestros.genes <- data.frame(table(gs))
  id <- seleccionarFilas(as.character(gs.nuestros.genes[,1]), as.character(selectedGenes[,1]))
  gs.nuestros.genes <- gs.nuestros.genes[id,]
  as.data.frame(gs.nuestros.genes[order(gs.nuestros.genes[,2], decreasing = T),])


  t <- data.frame(num.apariciones.genes = gs.nuestros.genes[,2])

  ggplot2::ggplot(data=t, ggplot2::aes(num.apariciones.genes)) +
          ggplot2::geom_histogram(col="light blue",fill="light blue",alpha=0.75, breaks=seq(0, 10, by=1)) +
          ggplot2::labs(title="Histograma frecuencia genes", x = "Num de veces que se repiten los genes", y="Num de genes")



  return(summary(gs.nuestros.genes[,2]))

}


