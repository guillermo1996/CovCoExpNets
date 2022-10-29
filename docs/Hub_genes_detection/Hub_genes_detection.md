-   [Hub Gene detection algorithm](#hub-gene-detection-algorithm)
    -   [Requisites](#requisites)
-   [Summary](#summary)
-   [Steps](#steps)
    -   [Step 1: Executing the `glmnet`
        repetitions](#step-1-executing-the-glmnet-repetitions)
    -   [Step 2: Minimum relative ratio of appearance
        threshold](#step-2-minimum-relative-ratio-of-appearance-threshold)
    -   [Step 3: Model generation](#step-3-model-generation)
    -   [Step 4: Genes and coefficients
        extraction](#step-4-genes-and-coefficients-extraction)
-   [Functions employed:](#functions-employed)
-   [Train/test split](#traintest-split)
    -   [Traditional split](#traditional-split)
    -   [Bootstrap](#bootstrap)
    -   [Hidden split](#hidden-split)
-   [Multiple conditions](#multiple-conditions)
-   [Predicting the sex](#predicting-the-sex)

Hub Gene detection algorithm
============================

In this tutorial, we will show the pipeline to extract the hub genes (or
predictors) for a given covariate. These genes will be used to generate
the co-expression networks.

Requisites
----------

The libraries required for this tutorial are `CovCoExpNets`, `magrittr`
for the `%>%` pipe-like operator, `dplyr` and `logger`. Additionally,
`CovCoExpNets` uses `foreach` and `doParallel` to execute the functions
in parallel.

``` r
library(CovCoExpNets)
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 4.1-2
#> Loading required package: foreach
#> Loading required package: doParallel
#> Loading required package: iterators
#> Loading required package: parallel
library(magrittr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(logger)

#doParallel::registerDoParallel(13)
cl <- makeCluster(13)
doParallel::registerDoParallel(cl)
```

As for the input data, we will use the preprocessed data obtained from
the [Data extraction tutorial](../Data_preparation). We need a data
matrix with genes as rows and samples as columns, as well as a numeric
vector for the covariate. In our case, we will use the data and age for
the Cortex tissue, but this tutorial will work if a list of data
matrices and a list of covariate vectors were used.

``` r
brain_data.path <- "~/GTEx_data/"

data = readRDS(paste0(brain_data.path, "data.combined.rds"))
age = readRDS(paste0(brain_data.path, "age.combined.rds"))

m <- lapply(age, function(x) x$mean)
d <- lapply(age, function(x) x$standard.deviation)
age <- lapply(age, function(x) x$covariate)

data = data[["Cortex"]]
age = age[["Cortex"]]
m = m[["Cortex"]]
d = d[["Cortex"]]
```

Summary
=======

The whole pipeline to obtain the hub genes is as follows:

``` r
genes.freq = CovCoExpNets::geneFrequency(data, age, t = 10, seed = 1796)
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9, relative = TRUE)
cvfit = CovCoExpNets::glmnetGenesSubset(data, age, genes.subset, seed = 1796, evaluate.model = F)
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)
```

The minimum requirements to execute the pipeline are the dataset and the
covariate vector, the *t* hyperparameter (by default 10) and the *mrfa*
hyperparameter (by default 0.9). The output is a list of the hub genes:

``` r
# First 25 relevant genes for Cortex tissue:
genes.relevant$Genes[1:25]
#>  [1] EDA2R      ADRA2B     BAIAP2L2   ZNF229     C12orf60   GPR26     
#>  [7] POC1A      HAP1       IQGAP3     FATE1      MCF2L2     IGF1      
#> [13] TLX2       ST6GALNAC2 SYT6       SCUBE1     ETV4       ABCB4     
#> [19] SOX7       LBHD1      RPSA       PARPBP     TSPOAP1    ZP3       
#> [25] ZKSCAN8P1 
#> 214 Levels: AARD ABCB4 ADAD2 ADAMTS7 ADRA2B AHRR ALDH3A1 ALOX15 ALOX15B ... ZP3
```

Steps
=====

Step 1: Executing the `glmnet` repetitions
------------------------------------------

The first step is to execute the `glmnet` repetitions a set number of
times (`t`). To do so, we will use the `geneFrequency` function from
`CovCoExpNets`.

``` r
genes.freq = CovCoExpNets::geneFrequency(data, age, t = 10, seed = 1796)
```

As seen in the example, only the `data` and `covariate` variables are
required by default. The result will be a data.frame with all genes
selected by each iteration and their coefficients.

``` r
head(genes.freq, 5)
#>      Genes Coefficients iter
#> 1    EDA2R    0.2500961    1
#> 2   ADRA2B   -0.1690084    1
#> 3 BAIAP2L2    0.1617018    1
#> 4   ZNF229   -0.1608639    1
#> 5    GPR26   -0.1472022    1
tail(genes.freq, 5)
#>       Genes  Coefficients iter
#> 2438 PCDHB8  4.974914e-05   10
#> 2439    LPA -9.606374e-06   10
#> 2440   LSM6  7.456103e-06   10
#> 2441  KYAT3 -1.796169e-06   10
#> 2442  TNNT2  1.633166e-06   10
```

Additionally, we can use the argument `iter.RMSE` to calculate the RMSE
of each `glmnet` repetition. To do so, we need to provide additional
testing data with `data.test.extra` and `covariate.test.extra` or to
input a `train.split < 1`.

Step 2: Minimum relative ratio of appearance threshold
------------------------------------------------------

In this step, we reduce the list of predictors returned in each
iteration. To do so, we use the hyperparameter *mrfa* to specify the
minimum relative frequency in which a predictor must appear. To do so,
we use the function `reduceGenes` from `CovCoExpNets`, where we count
how many times each gene is selected as relevant across the `glmnet`
repetitions and divide it by the total number of repetitions so that we
have the percentage of repetitions in which it was selected as relevant.
We require a minimum appearance threshold.

By default, we recommend an *mrfa* of 0.9:

``` r
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9, relative = T)

# First 25 genes:
genes.subset[1:25]
#>  [1] "AARD"      "ABCB4"     "ADAD2"     "ADAMTS7"   "ADRA2B"    "AHRR"     
#>  [7] "ALDH3A1"   "ALOX15"    "ALOX15B"   "ALPL"      "AMZ1"      "ANKRD18B" 
#> [13] "ANKRD33"   "APC"       "APOBEC3A"  "APOBEC3C"  "APOC1"     "ARHGAP11A"
#> [19] "ARL17A"    "ASCL2"     "ASF1B"     "ATAD3C"    "BAIAP2L2"  "C12orf54" 
#> [25] "C12orf60"
```

The output will be a set of genes that pass the threshold. In total,
there were 259 different genes across all repetitions, while 230 of
those passed the threshold.

Step 3: Model generation
------------------------

Next, we create the final `glmnet` model. We use the `glmnetGenesSubset`
function from `CovCoExpNets`.

``` r
cvfit = CovCoExpNets::glmnetGenesSubset(data, age, genes.subset, seed = 1796)

cvfit
#> 
#> Call:  glmnet::cv.glmnet(x = data.train, y = covariate.train, nfolds = k.folds,      alpha = 1, family = glmnet.family) 
#> 
#> Measure: Mean-Squared Error 
#> 
#>       Lambda Index Measure       SE Nonzero
#> min 0.001867    60 0.06773 0.008122     214
#> 1se 0.002972    55 0.07581 0.006457     218
```

The resulted final model can be evaluated in two ways:

1.  If the argument `evalute.model` and a `train.split` argument smaller
    than 1 are provided to the `glmnetGenesSubset()` function, the
    returned object will be a list of the model (`cvfit`) and the
    evaluation (`evaluation`).
2.  If the train/test split is done outside of `glmnetGenesSubset()`,
    then use the function `evaluateModel()`, with the additional test
    dataset as arguments (example in later
    [section](#Train/test-split)).

Step 4: Genes and coefficients extraction
-----------------------------------------

The final step is to extract the hub genes from the generated model. We
use the `extractModelGenes` function from `glmnet`.

``` r
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)

head(genes.relevant, 10)
#>       Genes Coefficients
#> 1     EDA2R    0.2532954
#> 2    ADRA2B   -0.1844454
#> 3  BAIAP2L2    0.1692918
#> 4    ZNF229   -0.1641643
#> 5  C12orf60   -0.1329406
#> 6     GPR26   -0.1187161
#> 7     POC1A   -0.1113262
#> 8      HAP1    0.1096197
#> 9    IQGAP3    0.1046102
#> 10    FATE1    0.1016612
```

The list of hub genes is found in the “Genes” column of `genes.relevant`
variables. The coefficients are also reported.

Functions employed:
===================

The required functions to execute the pipeline are:

-   **geneFrequency:** executes the `glmnet` repetitions. It has the
    following arguments:

    -   **data**: numeric data matrix with genes as rows and samples as
        columns. It also accepts a list of data matrices.
    -   **covariate:** numeric vector with the age of the sample’s
        donors. It also accepts a list of vectors.
    -   **t:** number of `glmnet` repetitions to execute. Its value
        depends on the stability of the `glmnet` repetitions. Deaults
        to 10.
    -   **k.folds:** (optional) number of k-folds in the `glmnet`
        cross-validation function `cv.glmnet`. Defaults to 10.
    -   **train.split:** (optional) percentage to split the training
        dataset for testing purposes. A value of 0.8 will set the 80% of
        the samples to train and 20% to estimate the test RMSE. Defaults
        to 1.
    -   **iter.RMSE:** (optional) whether to measure the RMSE of each
        repetition. Different solutions will not only have different
        returned predictors, but also different performance. Use this
        setting to prioritize better RMSE iterations when selecting the
        pruned predictors.
    -   **data.test.extra:** (optional) numeric data matrix with genes
        as rows and samples as columns. If the train/test split was
        executed before this step and you want to estimate the RMSE over
        the test dataset, provide the test dataset in here. By default,
        none.
    -   **covariate.test.extra:** (optional) numeric vector with the age
        of the sample’s donors. If the train/test split was executed
        before this step and you want to estimate the RMSE over the test
        dataset, provide the test dataset in here. By default, none.
    -   **sample.prob:** (optional) when executing the train/test split,
        it dictates the probability of each sample being executed. See
        [sample](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sample)
        for more information.
    -   **seed:** (optional) numeric seed to ensure reproducibility.
        Because of parallelization, even with an input seed results
        might differ.
    -   **glmnet.family:** (optional) family parameter used in the
        `cv.glmnet` function. See
        [`glmnet`](https://www.rdocumentation.org/packages/glmnet/versions/4.1-4/topics/cv.glmnet)
        for more details about the `family` parameter.

-   **reduceGenes:** executes the minimum relative frequency of
    appearance (*mrfa*) threshold. It has the following arguments:

    -   **genes.freq:** data.frame containing the genes that were
        selected in the `glmnet` repetitions. It must be obtained with
        the `geneFrequency` function from `CovCoExpNets`.
    -   **mrfa:** numeric minimum relative frequency of appearance. If
        between
    -   **force.mrfa:** (optional) whether to use only the *mrfa*
        parameter given. If set to FALSE, the parameter can be lowered
        until at least two predictors are selected. By default, TRUE.
    -   **relative:** (optional) whether to use a relative *mrfa* or
        not. If TRUE, the parameter will be considered as a percentage
        from 0 to 1. If FALSE, the parameter will be considered the
        plain number of times the predictor must appear from 1 to t. By
        default, TRUE.

-   **glmnetGenesSubset:** executes the final model training. It has the
    following arguments:

    -   **data**: numeric data matrix with genes as rows and samples as
        columns. It also accepts a list of data matrices.
    -   **covariate:** numeric vector with the age of the sample’s
        donors. It also accepts a list of vectors.
    -   **genes.subset:** character vector with the genes to be used to
        train the model. It must be the output of `reduceGenes` from
        `CovCoExpNets`. It also accepts a list of vectors.
    -   **k.folds:** (optional) number of k-folds in the `glmnet`
        cross-validation function `cv.glmnet`. Defaults to 10.
    -   **train.split:** (optional) percentage to split the training
        dataset for testing purposes. A value of 0.8 will set the 80% of
        the samples to train and 20% to estimate the test RMSE. Defaults
        to 1.
    -   **sample.prob:** (optional) when executing the train/test split,
        it dictates the probability of each sample being executed. See
        [sample](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sample)
        for more information.
    -   **seed:** (optional) numeric seed to ensure reproducibility.
        Because of parallelization, even with an input seed results
        might differ.
    -   **glmnet.family:** (optional) family parameter used in the
        `cv.glmnet` function. See
        [`glmnet`](https://www.rdocumentation.org/packages/glmnet/versions/4.1-4/topics/cv.glmnet)
        for more details about the `family` parameter.
    -   **evaluate.model:** (optional) whether to evaluate the generated
        models using train/test splits with the train.split argument.
    -   **m:** (optional) list of means for the covariate.
    -   **d:** (optional) list of standard deviations for the covariate.

-   **extractModelGenes:** extract the relevant predictors from a
    `glmnet` model. It has the following arguments.

    -   **cvfit:** `glmnet` model to extract the relevant predictors and
        their coefficients from.
    -   **genes.freq:** data.frame obtained from `geneFrequency`. If
        provided, it will also return the average coefficient of each
        gene across all glmnet repetitions. Defaults to none.

Train/test split
================

The train/test split can be done in several ways:

Traditional split
-----------------

Most of the functions are compatible with the use of separate train and
test dataset:

``` r
set.seed(1796)
train.idx = sample(1:ncol(data), 0.75*ncol(data))

data.train = data[, train.idx]
data.test = data[, -train.idx]

age.train = age[train.idx]
age.test = age[-train.idx]

genes.freq = geneFrequency(data.train, age.train, t = 10, seed = 1796, data.test.extra = data.test, covariate.test.extra = age.test, iter.RMSE = T)
genes.subset = reduceGenes(genes.freq, mrfa = 0.9)
cvfit = glmnetGenesSubset(data.train, age.train, genes.subset, seed = 1796, evaluate.model = F)
evaluateModel(cvfit, data.test, age.test, genes.subset, data.train = data.train, covariate.train = age.train, m = m, d = d)
#>   Condition rmse.test rmse.train rmse.bootstrap rmse.null Samples.train
#> 1 Condition  6.666868   1.033669       4.593851   9.30509           191
#>   Samples.test Initial.predictors Returned.predictors  r2.train  r2.test
#> 1           64              15516                 108 0.9900352 0.484671
#>   r2_adj.train r2_adj.test
#> 1    0.9769107    1.721461
```

Bootstrap
---------

The `CovCoExpNets` package allows for the generation of datasets via
bootstrapping technique:

``` r
data.bootstrap = splitDataBootstrap(data, age, seed = 1796)

data.train = data.bootstrap$data.train
data.test = data.bootstrap$data.test

age.train = data.bootstrap$covariate.train
age.test = data.bootstrap$covariate.test
```

Hidden split
------------

It is also possible to use the `train.split` argument inside most of the
functions to ignore the train/test split directly. It is important then
to keep a consistent `seed` across the different functions, or the
train/test splits will be different:

``` r
genes.freq = geneFrequency(data, age, t = 10, train.split = 0.75, seed = 1796, iter.RMSE = T)
genes.subset = reduceGenes(genes.freq, mrfa = 0.9)
glmnet.models = glmnetGenesSubset(data, age, genes.subset, train.split = 0.75, seed = 1796, evaluate.model = T, m = m, d = d)

cvfit = glmnet.models$cvfit

glmnet.models$evaluation
#>   Condition rmse.test rmse.train rmse.bootstrap rmse.null Samples.train
#> 1 Condition  7.234244   1.701644       5.198247   9.30509           191
#>   Samples.test Initial.predictors Returned.predictors r2.train   r2.test
#> 1           64              15516                  82 0.972995 0.3932257
#>   r2_adj.train r2_adj.test
#> 1    0.9524911    3.011936
```

Multiple conditions
===================

As mentioned in the first section, we can input several covariates at
the same time in the form of a list. The following pipeline will
generate the hub genes for both Cerebellar Hemisphere and Cerebellum
tissues.

``` r
data = readRDS(paste0(brain_data.path, "data.combined.rds"))
age = readRDS(paste0(brain_data.path, "age.combined.rds"))

m <- lapply(age, function(x) x$mean)
d <- lapply(age, function(x) x$standard.deviation)
age <- lapply(age, function(x) x$covariate)

data = data[c("Cerebellar Hemisphere", "Cerebellum")]
age = age[c("Cerebellar Hemisphere", "Cerebellum")]
m = m[c("Cerebellar Hemisphere", "Cerebellum")]
d = d[c("Cerebellar Hemisphere", "Cerebellum")]

genes.freq = CovCoExpNets::geneFrequency(data, age, t = 10, seed = 1796)
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9)
cvfit = CovCoExpNets::glmnetGenesSubset(data, age, genes.subset, seed = 1796)
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)
```

``` r
lapply(genes.relevant, function(x) head(x, 5))
#> [[1]]
#>     Genes Coefficients
#> 1 FAM126A   -0.1687862
#> 2   SYT16    0.1524898
#> 3   RGS19    0.1491559
#> 4   MTURN   -0.1444758
#> 5  SLC9A2   -0.1289828
#> 
#> [[2]]
#>     Genes Coefficients
#> 1    WWC2   -0.1477740
#> 2   PTPN3   -0.1415790
#> 3  PLAGL1   -0.1243174
#> 4   WNT5B   -0.1230792
#> 5 SLC12A1    0.1230452
```

Predicting the sex
==================

The `CovCoExpNets` package is also functional for binary covariates,
like the sex of the samples’ donors. Using the `glmnet.family` across
the different functions, we can generate the most relevant genes to
predict the sex. It is important to notice that we cannot use the same
dataset as before, since they included sexual genes. To predict the age,
we are interested only in autosomal genes.

``` r
brain_data.path <- "~/GTEx_data/"

data.autosomes = readRDS(paste0(brain_data.path, "data.autosomes.combined.rds"))
sex = readRDS(paste0(brain_data.path, "sex.combined.rds"))

data.autosomes = data.autosomes[["Cortex"]]
sex = sex[["Cortex"]]
```

``` r
genes.freq = CovCoExpNets::geneFrequency(data.autosomes, sex, t = 10, seed = 1796, glmnet.family = "binomial")
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9, relative = TRUE)
cvfit = CovCoExpNets::glmnetGenesSubset(data.autosomes, sex, genes.subset, seed = 1796, glmnet.family = "binomial")
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)

head(genes.relevant, 10)
#>         Genes Coefficients
#> 1        TLX2     4.256895
#> 2       ACOT4    -2.383916
#> 3      SPESP1    -2.369456
#> 4        HAAO     2.244320
#> 5        WEE1     2.091729
#> 6       ADAD2     1.973587
#> 7       NLRP2     1.768200
#> 8  CSGALNACT1     1.692188
#> 9        MYL7    -1.421576
#> 10   PCDHGA12     1.405528
```

``` r
stopCluster(cl)
```
