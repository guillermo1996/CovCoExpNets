
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CovCoExpNets

<!-- badges: start -->

<!-- badges: end -->

Nowadays, the way to create a co-expression network is to use
Hierarchical Clustering. This package allows you to create a
co-expression network using the glmnet algorithm. In addition, the
network created will be relative to a covariate of the sample to be
studied, so we are creating a supervised coexpression network.

## Installation

CovCoExpNets can be installed with the following command:

``` r
devtools::install_github("carmen-maria-hernandez/SuCoNets")
```

# Example of use

## Hub gene detection

Suppose we have an expression matrix, `data`, with blood samples as
columns and genes as rows. Each sample is identified by the numerical
values taken by the genes. The particular covariate we are going to
study is the age of the samples’ donors, given in a numerical vector
denotated as `age`.

The hub gene detection algorithm is based on Lasso repetitions, where
the random effects of highly correlated variables are reduced by
applying a certain number of repetitions `t`. The genes selected as
relevant will constitute the final dataset to execute a Lasso
regularization again. From the final model generated, we extract the
relevant genes and their coefficients as their relevance to the given
covariate.

### Data preprocessing

The process starts by preprocessing the data. The specific details of
this process are shown in its respective tutorial ([Data
preparation](docs/Data_preparation)). As a summary, the dataset is first
transformed to a logarithmic scale. Then, we require a minimum
activation of 0.1 in at least 80% of the samples, followed by the
removall of all non protein coding genes and those with a low variation
across the samples. Lastly, we centralize and normalize the expression
matrix. All of these steps are encompassed in the `dataPreprocess()`
function, where every step can be individually modified or executed as
requested.

``` r
library(CovCoExpNets)
library(dplyr)

raw_data <- readRDS("raw_data.rds")
data <- dataPreprocess(raw_data, includeSexChromosomes = T)
```

As seen in the example, we can also filter by genes loacted in autosomal
chromosomes. The specific covariate to study, the age, also needs to be
normalized with the `normalize()` function:

``` r
raw_age <- readRDS("raw_age.rds")
age <- normalize(raw_age)

m <- age$mean
d <- age$standard.deviation
age <- age$covariate
```

The `normalize()` function also returns the mean and the standard
deviation of the vector. This data will be later used to restore the
predictions to an age scale.

### `Glmnet` repetitions

The next step in the hub gene detection is to execute the Lasso
repetitions. The chosen Lasso implementation is the `glmnet` algorithm.
With the function `geneFrequency()`, we set the number of repetitions
(use `t=10` as a baseline). We can also set the train/test split and the
initial seed to ensure reproducible results. We will then use the
function `reduceGenes()` to specify the minimum relative frequency of
appearance that a gene must have to be selected for the final model
execution. In this example, we will set this hyperparameter to `mrfa
= 0.9`.

``` r
genes_freq <- geneFrequency(data, age, t = 10, train.split = 0.8, seed = 1796)
genes_subset <- reduceGenes(genes_freq, mrfa = 0.9, relative = T)
```

### Final model generation

Once we have a first selection of relevant genes, we execute the
function `glmnetGenesSubset()` to calcualte the final model and extract
the relevant genes from
it.

``` r
glmnet_model <- glmnetGenesSubset(data, age, genes_subset, train.split = 0.8, seed = 1796, evaluate.model = T, m = m, d = d)

cvfit <- glmnet_model$cvfit
evaluation <- glmnet_model$evaluation
```

We set the same `train.split` and `seed` as for the `glmnet` repetitions
to ensure that the same training dataset was used to generate the model.
Since we set the argument `evaluate.model` to `TRUE`, the returned
predictor will contain both the generated model and the evaluation
(i.e. RMSE). With information about the mean and standard deviation
provided, the results will be in
years.

``` r
evaluation %>% mutate(across(where(is.numeric), round, 3)) %>% knitr::kable()
```

| Condition | rmse.test | rmse.train | rmse.bootstrap | rmse.null | Samples.train | Samples.test | Initial.predictors | Returned.predictors | r2.train | r2.test | r2\_adj.train | r2\_adj.test |
| :-------- | --------: | ---------: | -------------: | --------: | ------------: | -----------: | -----------------: | ------------------: | -------: | ------: | ------------: | -----------: |
| Condition |     5.676 |      0.322 |          3.706 |     9.004 |           204 |           51 |              15516 |                 173 |    0.999 |   0.598 |         0.993 |        1.163 |

### Genes extraction

To extract the relevant genes from the model, we use the function
`extractModelGenes`:

``` r
relevant_genes <- extractModelGenes(cvfit)
knitr::kable(relevant_genes %>% head(5))
```

| Genes    | Coefficients |
| :------- | -----------: |
| EDA2R    |    0.2877469 |
| BAIAP2L2 |    0.1573385 |
| ADRA2B   |  \-0.1251515 |
| GPR26    |  \-0.1232131 |
| GFAP     |    0.1187436 |

## Coexpression network generation

# Other tutorials

Other tutorials available to better understand how this package works
are found in the *docs* directory:

  - [Data preparation](docs/Data_preparation)
  - [Hub genes detection](docs/Hub_genes_detection)
  - [Covariate simulations](docs/Simulation_framework)
  - [`Glmnet` stability studies](docs/Stability_glmnet)

# Credits

This package is based on the package *glmnet*, available at the
following URL:
<https://cran.r-project.org/web/packages/glmnet/glmnet.pdf>. It is also
a fork of the package *SuCoNets*, developed by Carmen María Hernandez
Both *CovCoExpNets* and *SuCoNets* have been supervised by both Juan A.
Botía (Universidad de Murcia) (<https://github.com/juanbot>) and Alicia
Gómez Pascual (Universidad de Murcia), who also contributed to its
design.
