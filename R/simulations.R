P = 10; Q1 = 0; Q2 = 2; Qc = 0; distribution = "uniform"; params = list(min = -2, max = 2, mean = 0, sd = 1); seed = sample(1:99999999, 1); share.coeffs = F
#' Generate simulations coefficients
#'
#' Generate the coefficients for the different predictors generated in a simulation.
#' The coefficients can be generated following a uniform or gaussian distribution.
#'
#' @param P number of predictors to simulate.
#' @param Q1 number of relevant predictors for the first covariate.
#' @param Q2 number of relevant predictors for the second covariate.
#' @param Qc number of common relevant predictors between the two predictors. It must
#'   be smaller than both Q1 and Q2.
#' @param distribution string indicating the distribution from which to take the
#'   coefficients. Only uniform and gaussian distributions are possible. By default uniform
#' @param params list with the parameters from which to generate the coefficients. If
#'   the distribution is uniform, it requires the minimum and maximum values. If the
#'   distribution is gaussian, it requires the mean and standard deviation.
#' @param seed numeric seed to ensure reproducibility
#' @param share.coeffs boolean to indicate if the common predictors between the two
#'   covariates should have the same coefficients. By default FALSE
#'
#' @return coefficients for each covariate and which predictor is relevant to each one
#' @export
generateWeights <- function(P, Q1, Q2, Qc, distribution = "uniform", params = list(min = -2, max = 2, mean = 0, sd = 1), seed = sample(1:99999999, 1), share.coeffs = F){
  set.seed(seed)
  splits.raw = split(sample(seq(1, P), Q1+Q2-Qc), rep(1:3, c(Q1-Qc, Q2-Qc, Qc)))

  splits = list()
  splits["Q1"] = splits.raw['1']
  splits["Q2"] = splits.raw['2']
  splits["Qc"] = splits.raw['3']

  relevant.Q1 = c(splits$Q1, splits$Qc)
  relevant.Q2 = c(splits$Q2, splits$Qc)

  if(distribution == "uniform"){
    if(all(c("min", "max") %in% names(params))){
      coefficients.Q1 = runif(P, min = params$min, max = params$max)
      coefficients.Q2 = runif(P, min = params$min, max = params$max)
    }else{
      stop('For uniform distributions, a parameter list must be provided with "min" and "max" value. Example: params = list(min = -2, max = 2)')
    }
  }else if(distribution == "gaussian"){
    if(all(c("mean", "sd") %in% names(params))){
      coefficients.Q1 = rnorm(P, mean = params$mean, sd = params$sd)
      coefficients.Q2 = rnorm(P, mean = params$mean, sd = params$sd)
    }else{
      stop('For gaussian distributions, a parameter list must be provided with "mean" and "sd" value. Example: params = list(mean = 0, sd = 1)')
    }
  }
  if(share.coeffs) coefficients.Q2 = coefficients.Q1

  weights.Q1 = rep(0, P)
  weights.Q2 = rep(0, P)

  weights.Q1[relevant.Q1] = coefficients.Q1[relevant.Q1]
  weights.Q2[relevant.Q2] = coefficients.Q2[relevant.Q2]

  return(list(weights.Q1 = weights.Q1, weights.Q2 = weights.Q2, splits = splits))
}

#' Generate simulation formula
#'
#' Generate the formula that will follow the simulations. It will produce a linear
#' relationship between the response variable y and the predictors x:
#' y ~ B1*x1 + B2*x2
#'
#' @param P number of predictors to simulate.
#' @param start numberic, where to begin the formula. By default, 1
#'
#' @return linear formula that the simulation will follow
#' @export
generateFormula <- function(P, start = 1){
  formula = "y ~ 0"
  for(i in 1:P){
    formula <- paste0(formula, " +x", i + start - 1)
  }
  return(as.formula(formula))
}


#' Generate simulation predictors
#'
#' Indicates how to create each predictor to the simulation algorithm.
#' It will produce a list of predictors with the instructions to be generated
#' as a normal-like distribution with mean 0 and standard deviation 1.
#'
#' @param P number of predictors to simulate.
#' @param start numberic, where to begin the formula. By default, 1
#'
#' @return list of instructions to generate the predictors.
#' @export
generateFixed <- function(P, start = 1){
  fixed = vector(mode = "list", length = P)
  names <- vector(mode = "list", length = P)
  for(i in 1:P){
    fixed[[i]] <- list(var_type = "continuous", mean = 0, sd = 1)
    names[[i]] <- paste0("x", i + start - 1)
  }
  names(fixed) <- names
  return(fixed)
}


#' Generate a bivariate simulation
#'
#' Generates the simulation object. It includes the data matrix, the numeric
#' response vectors and additional information, like the coefficients of the
#' relevant predictors for each covariate. It is based on the package \link[simglm]{simglm}.
#'
#' @param N number of samples to simulate.
#' @param P number of predictors to simulate.
#' @param Q1 number of relevant predictors for the first covariate.
#' @param Q2 number of relevant predictors for the second covariate.
#' @param Qc number of common relevant predictors between the two predictors. It must
#'   be smaller than both Q1 and Q2.
#' @param distribution string indicating the distribution from which to take the
#'   coefficients. Only uniform and gaussian distributions are possible. By default uniform
#' @param params list with the parameters from which to generate the coefficients. If
#'   the distribution is uniform, it requires the minimum and maximum values. If the
#'   distribution is gaussian, it requires the mean and standard deviation.
#' @param seed numeric seed to ensure reproducibility
#' @param share.coeffs boolean to indicate if the common predictors between the two
#'   covariates should have the same coefficients. By default FALSE
#'
#' @return simulation object from where we can extract the data matrix and the two
#'   response vectors.
#' @export
simulateBivariate <- function(N, P, Q1, Q2, Qc, distribution = "uniform", params = list(min = -2, max = 2), seed = sample(1:99999999, 1), share.coeffs = F){
  weights = generateWeights(P, Q1, Q2, Qc, distribution, params, seed, share.coeffs)
  weights.Q1 = weights$weights.Q1
  weights.Q2 = weights$weights.Q2
  splits = weights$splits

  sim_arguments.1 <- list(formula = as.formula(generateFormula(P)),
                          fixed = generateFixed(P),
                          sample_size = N,
                          error = list(variance = 1),
                          reg_weights = weights.Q1)


  sim_arguments.2 <- list(formula = as.formula(generateFormula(P)),
                          fixed = generateFixed(P),
                          sample_size = N,
                          error = list(variance = 1),
                          reg_weights = weights.Q2)

  set.seed(seed)
  sim.fixed = simglm::simulate_fixed(data = NULL, sim_arguments.1)

  sim.obj.1 = sim.fixed %>%
    simglm::simulate_error(sim_arguments.1) %>%
    simglm::generate_response(sim_arguments.1)

  sim.obj.2 = sim.fixed %>%
    simglm::simulate_error(sim_arguments.2) %>%
    simglm::generate_response(sim_arguments.2)

  data = sim.obj.1 %>% select(starts_with(c("x")))

  response.1 = sim.obj.1$y
  response.2 = sim.obj.2$y

  df.coefficients.1 <- data.frame(predictor = names(sim_arguments.1$fixed)[which(sim_arguments.1$reg_weights != 0)],
                                  coefficient = sim_arguments.1$reg_weights[which(sim_arguments.1$reg_weights != 0)])

  df.coefficients.2 <- data.frame(predictor = names(sim_arguments.2$fixed)[which(sim_arguments.2$reg_weights != 0)],
                                  coefficient = sim_arguments.2$reg_weights[which(sim_arguments.2$reg_weights != 0)])

  predictor.splits = list(Q1 = names(sim_arguments.1$fixed)[c(splits$Q1, splits$Qc)],
                          Q1.unique = names(sim_arguments.1$fixed)[splits$Q1],
                          Q2 = names(sim_arguments.1$fixed)[c(splits$Q2, splits$Qc)],
                          Q2.unique = names(sim_arguments.1$fixed)[splits$Q2],
                          Qc = names(sim_arguments.1$fixed)[splits$Qc])

  return(list(data = data,
              y1 = response.1,
              y2 = response.2,
              df.coefficients.1 = df.coefficients.1,
              df.coefficients.2 = df.coefficients.2,
              predictor.splits = predictor.splits))
}

#' Generate a bivariate simulation with large number of predictors
#'
#' Generates the simulation object when there are too many predictors for R to
#' handle properly. It splits the generation process in chunks of "max.predictors"
#' (by default, 15000). The main limitation is the formula size, which is limited
#' by default to a maximum value of arround 17000. Because of these limitations,
#' the relevant predictors will only be present in the first "max.predictors". That is,
#' if 1000 predictors are marked as relevant, they will be by default in the first 15000
#' predictors generated. Thus, the maximum number of relevant predictors is the argument
#' max.predictors. For more information about the simulation process, please refer to
#' \link[CovCoExpNets]{simulateBivariate}
#'
#' @param N number of samples to simulate.
#' @param P number of predictors to simulate.
#' @param Q1 number of relevant predictors for the first covariate.
#' @param Q2 number of relevant predictors for the second covariate.
#' @param Qc number of common relevant predictors between the two predictors. It must
#'   be smaller than both Q1 and Q2.
#' @param distribution string indicating the distribution from which to take the
#'   coefficients. Only uniform and gaussian distributions are possible. By default uniform.
#' @param params list with the parameters from which to generate the coefficients. If
#'   the distribution is uniform, it requires the minimum and maximum values. If the
#'   distribution is gaussian, it requires the mean and standard deviation.
#' @param seed numeric seed to ensure reproducibility.
#' @param share.coeffs boolean to indicate if the common predictors between the two
#'   covariates should have the same coefficients. By default FALSE.
#' @param max.predictors numeric maximum number of predictors per chunk of simulation.
#'   By default, 15000.
#'
#' @return simulation object from where we can extract the data matrix and the two
#'   response vectors.
#' @export
simulateBivariateLarge <- function(N, P, Q1, Q2, Qc, distribution = "uniform", params = list(min = -2, max = 2), seed = sample(1:99999999, 1), share.coeffs = F, max.predictors = 15000){
  if(P <= max.predictors){
    stop(paste0("The number of predictors given (", P, ") is not larger than the maximun number of predictors (", max.predictors, "). Please, use the \"simulateBivariate()\" function"))
  }
  P.splits = split(seq(1, P), ceiling(seq_along(seq(1: P))/max.predictors))

  sim.obj = simulateBivariate(N, length(P.splits[[1]]), Q1, Q2, Qc, distribution, params, seed, share.coeffs)

  extra.data = foreach(i = 2:length(P.splits), .combine = "cbind") %do%{
    P.length = length(P.splits[[i]])
    P.start = P.splits[[i]][1]

    sim_arguments <- list(formula = as.formula(generateFormula(P.length, P.start)),
                            fixed = generateFixed(P.length, P.start),
                            sample_size = N)
    sim.fixed = simglm::simulate_fixed(data = NULL, sim_arguments) %>% select(starts_with(c("x")))
    sim.fixed
  }

  sim.obj$data = cbind(sim.obj$data, extra.data)
  sim.obj
}

#' Generate an univariate simulation
#'
#' Generates the simulation object. It includes the data matrix, the numeric
#' response vector and additional information, like the coefficients of the
#' relevant predictors. It is based on the package \link[simglm]{simglm}. It is
#' a wrapper of the \link[CovCoExpNets]{simulateBivariate} function.
#'
#' @param N number of samples to simulate.
#' @param P number of predictors to simulate.
#' @param Q number of relevant predictors.
#' @param distribution string indicating the distribution from which to take the
#'   coefficients. Only uniform and gaussian distributions are possible. By default uniform.
#' @param params list with the parameters from which to generate the coefficients. If
#'   the distribution is uniform, it requires the minimum and maximum values. If the
#'   distribution is gaussian, it requires the mean and standard deviation.
#' @param seed numeric seed to ensure reproducibility.
#'
#' @return simulation object from where we can extract the data matrix and the
#'   response vector.
#' @export
simulateUnivariate <- function(N, P, Q, distribution = "uniform", params = list(min = -2, max = 2), seed = sample(1:99999999, 1)){
  sim.obj = simulateBivariate(N, P, Q, 0, 0, distribution, params, seed)
  sim.obj = sim.obj[c("data", "y1", "df.coefficients.1", "predictor.splits")]
  names(sim.obj) <- c("data", "y", "df.coefficients", "predictor.splits")

  return(sim.obj)
}


#' Generate an univariate simulation with large number of predictors
#'
#' Generates the simulation object when there are too many predictors for R to
#' handle properly. It splits the generation process in chunks of "max.predictors"
#' (by default, 15000). The main limitation is the formula size, which is limited
#' by default to a maximum value of arround 17000. Because of these limitations,
#' the relevant predictors will only be present in the first "max.predictors". That is,
#' if 1000 predictors are marked as relevant, they will be by default in the first 15000
#' predictors generated. Thus, the maximum number of relevant predictors is the argument
#' max.predictors. For more information about the simulation process, please refer to
#' \link[CovCoExpNets]{simulateUnivariate}. It is a wrapper of the
#' \link[CovCoExpNets]{simulateBivariateLarge} function.
#'
#' @param N number of samples to simulate.
#' @param P number of predictors to simulate.
#' @param Q number of relevant predictors.
#' @param distribution string indicating the distribution from which to take the
#'   coefficients. Only uniform and gaussian distributions are possible. By default uniform.
#' @param params list with the parameters from which to generate the coefficients. If
#'   the distribution is uniform, it requires the minimum and maximum values. If the
#'   distribution is gaussian, it requires the mean and standard deviation.
#' @param seed numeric seed to ensure reproducibility.
#' @param max.predictors numeric maximum number of predictors per chunk of simulation.
#'   By default, 15000.
#'
#' @return simulation object from where we can extract the data matrix and the
#'   response vector.
#' @export
simulateUnivariateLarge <- function(N, P, Q, distribution = "uniform", params = list(min = -2, max = 2), seed = sample(1:99999999, 1), max.predictors = 15000){
  if(P <= max.predictors){
    stop(paste0("The number of predictors given (", P, ") is not larger than the maximun number of predictors (", max.predictors, "). Please, use the \"simulateUnivariate()\" function"))
  }
  sim.obj = simulateBivariateLarge(N, P, Q, 0, 0, distribution, params, seed, max.predictors = max.predictors)
  sim.obj = sim.obj[c("data", "y1", "df.coefficients.1", "predictor.splits")]
  names(sim.obj) <- c("data", "y", "df.coefficients", "predictor.splits")

  return(sim.obj)
}


#' Generates the confussion matrix for a simulation
#'
#' Given a set of predictors selected by the model (set.calculated)
#' and the simulation object to extract the total number of predictors,
#' this function generates the confussion matrix when discriminating
#' relevant vs non-relevant predictors.
#'
#' @param set.calculated vector of predictors selected by the model.
#' @param sim.obj simulation object obtained from CovCoExpNets.
#' @param covariate which covariate to study. Either "Q1" or "Q2". Defaults
#'   to Q2.
#'
#' @return caret confussion matrix.
#' @export
calculateSimulationConfussionMatrix <- function(set.calculated, sim.obj, covariate = "Q1"){
  one.hot.calculated = toOneHotEncoding(colnames(sim.obj$data), set.calculated$Genes)
  one.hot.original = toOneHotEncoding(colnames(sim.obj$data), sim.obj$predictor.splits[[covariate]])

  cm = caret::confusionMatrix(data = one.hot.calculated %>% factor(levels = 0:1),
                              reference = one.hot.original %>% factor(levels = 0:1),
                              positive = "1")
  return(cm)
}


#' Calculates different simulation model's metrics
#'
#' Given a set of predictors selected by the model (set.calculated)
#' and the simulation object to extract the total number of predictors,
#' this function calculates the precision, sensitivity, specificity, jaccard
#' index and the returned predictors.
#'
#' @param set.calculated vector of predictors selected by the model.
#' @param sim.obj simulation object obtained from CovCoExpNets.
#' @param covariate which covariate to study. Either "Q1" or "Q2". Defaults
#'   to Q2.
#'
#' @return simulation model's metrics.
#' @export
calculateSimulationMetrics <- function(set.calculated, sim.obj, covariate = "Q1"){
  cm = calculateSimulationConfussionMatrix(set.calculated, sim.obj, covariate)
  data.frame(Precision = cm$byClass["Precision"],
             Sensitivity = cm$byClass["Sensitivity"],
             Specificity = cm$byClass["Specificity"],
             Jaccard.index = jaccardIndex(set.calculated$Genes, colnames(sim.obj$data)),
             Returned.predictors = nrow(set.calculated),
             row.names = "Simulation")
}
