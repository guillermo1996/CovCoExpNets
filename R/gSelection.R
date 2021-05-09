#' Detect main genes
#'
#' Given an expression matrix and a covariate, this function calculates the variables
#' that best predict that covariate using glmnet algorithm
#'
#' @param data expression matrix
#' @param covariate numeric vector
#'
#' @return Dataframe containing the variables selected by the glmnet algorithm and the coefficients associated with these variables
#' @export
#'
detectGenes <- function(data, covariate){
  semilla = 0
  contador = 1;
  semillas = numeric(10)
  genesSelect = numeric(10)

  for(i in 9:0){
    if(semilla == 987654321)
      semilla = semilla +1
    else
      semilla = semilla*10+i
    set.seed(semilla)

    ind.train.10 <- sample(1:ncol(data), 0.8*ncol(data))
    data.train.10 <-  t(data[,ind.train.10])
    covariate.train.10 <-  covariate[ind.train.10]

    cvfit.10<- glmnet::cv.glmnet(data.train.10, covariate.train.10, alpha=1, family = "gaussian")

    coefic.10 <- as.matrix(stats::coef(cvfit.10,s="lambda.min"))
    semillas[contador] <- semilla
    genesSelect[contador] <- length(which(coefic.10!=0))-1
    contador = contador + 1

  }

  dat <- data.frame(semillas, genesSelect)

  semilla <- dat$semillas[which.min(dat$genesSelect)]

  set.seed(semilla)

  ind.train <- sample(1:ncol(data), 0.8*ncol(data))

  data.train <-  t(data[,ind.train])
  covariate.train <-  covariate[ind.train]

  data.test <- t(data[,-ind.train])
  covariate.test <- covariate[-ind.train]
  cvfit<- glmnet::cv.glmnet(data.train, covariate.train, alpha=1, family = "gaussian")
  BMTAR::print(cvfit)

  # Guardamos en un fichero el conjunto de genes seleccionados por glmnet para realizar la predicción:
  coefic <- as.matrix(stats::coef(cvfit,s="lambda.min"))
  filas.dist.cero <- rownames(coefic)[which(coefic!=0)]
  genes.seleccionados <- data.frame(Genes = filas.dist.cero[-1], Coeficientes = coefic[which(coefic!=0)][-1])

  # Predecimos la edad del conjunto de test
  predict.cvfit.test <- stats::predict(cvfit, newx = data.test,type = "response",s = "lambda.min")
  predict.cvfit.train <- stats::predict(cvfit, newx = data.train,type = "response",s = "lambda.min")
  # Predecimos la edad del conjunto de train


  # Calculamos el RMSE del modelo sobre los datos de train y sobre los datos de test
  tr <- MLmetrics::RMSE(predict.cvfit.train, covariate.train)
  te <- MLmetrics::RMSE(predict.cvfit.test, covariate.test)
  cat('\n',
      '-> RMSE del modelo sobre los datos de train: ', tr, '\n',
      '-> RMSE del modelo sobre los datos de test: ', te, '\n',
      sep = ''
  )
  return(genes.seleccionados)

}

#' Co-expression network
#'
#' Calculate a coexpression network using the genes selected by the glmnet algorithm
#' as the best predictors of a covariate.
#'
#' @param data expression matrix
#' @param genes dataframe with the genes selected by the glmnet algorithm and their coefficients
#'
#' @return A list where for each gene selected by glmnet appears the 50 genes most correlated with it and that correlation
#' @export
coexpressionNetwork <- function(data, genes){
  df <- calculateCorrelation(data, genes[1,1], 50)
  for (i in 2:nrow(genes)) {
    df <- rbind(df, calculateCorrelation(data, genes[i,1], 50))
  }

  all.genes = list()
  for (i in 1:nrow(genes)) {
    gen <- df[(i+(i-1)*50):(i+(i-1)*50+50), 2]
    corr <- df[(i+(i-1)*50):(i+(i-1)*50+50), 3]
    all.genes[[genes[i,1]]] = data.frame(gen,corr)
  }
  return(all.genes)
}
