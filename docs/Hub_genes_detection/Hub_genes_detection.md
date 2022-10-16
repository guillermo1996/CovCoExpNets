-   [Hub Gene detection algorithm](#hub-gene-detection-algorithm)
    -   [Requisites](#requisites)
    -   [Summary](#summary)
    -   [Steps](#steps)
    -   [Functions employed:](#functions-employed)
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
-------

The whole pipeline to obtain the hub genes is as follows:

``` r
genes.freq = CovCoExpNets::geneFrequency(data, age, t = 10, seed = 1796)
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9, relative = TRUE)
cvfit = CovCoExpNets::glmnetGenesSubset(data, age, genes.subset, seed = 1796)
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)
```

The minimum requirements to execute the pipeline are the dataset and the
covariate vector, the *t* hyperparameter (by default 10) and the *mrfa*
hyperparameter (by default 0.9). The output is a list of the hub genes:

``` r
# First 25 relevant genes for Cortex tissue:
genes.relevant$Genes[1:25]
#>  [1] EDA2R     ADRA2B    ZNF229    BAIAP2L2  C12orf60  GPR26     POC1A    
#>  [8] HAP1      FATE1     IQGAP3    MCF2L2    IGF1      PARPBP    KISS1R   
#> [15] LBHD1     SYT6      ABCB4     ETV4      QPRT      SOX7      ADAD2    
#> [22] TLX2      TFCP2L1   ZKSCAN8P1 SCUBE1   
#> 182 Levels: AARD ABCB4 ADAD2 ADRA2B AHRR ALDH3A1 ALOX15 ALOX15B ALPL ... ZP3
```

Steps
-----

### Step 1: Executing the `glmnet` repetitions

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
#> 1    EDA2R    0.2404529    1
#> 2    GPR26   -0.1660976    1
#> 3   ZNF229   -0.1582700    1
#> 4   ADRA2B   -0.1576799    1
#> 5 BAIAP2L2    0.1458428    1
tail(genes.freq, 5)
#>       Genes  Coefficients iter
#> 2354 PCDHB8  4.974914e-05   10
#> 2355    LPA -9.606374e-06   10
#> 2356   LSM6  7.456103e-06   10
#> 2357  KYAT3 -1.796169e-06   10
#> 2358  TNNT2  1.633166e-06   10
```

### Step 2: Minimum relative ratio of appearance threshold

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
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9)

# First 25 genes:
genes.subset[1:25]
#>  [1] "AARD"      "ABCB4"     "ADAD2"     "ADRA2B"    "AHRR"      "ALDH3A1"  
#>  [7] "ALOX15"    "ALOX15B"   "ALPL"      "AMZ1"      "ANKRD18B"  "ANKRD33"  
#> [13] "APC"       "APOC1"     "ARHGAP11A" "ASCL2"     "ATAD3C"    "BAIAP2L2" 
#> [19] "C12orf60"  "C1orf53"   "C1QL4"     "C1QTNF12"  "CACFD1"    "CAPN14"   
#> [25] "CARD16"
```

The output will be a set of genes that pass the threshold. In total,
there were 269 different genes across all repetitions, while 191 of
those passed the threshold.

### Step 3: Model generation

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
#>        Lambda Index Measure       SE Nonzero
#> min 0.0005074    74 0.03315 0.003555     182
#> 1se 0.0012865    64 0.03635 0.003583     180
```

### Step 4: Genes and coefficients extraction

The final step is to extract the hub genes from the generated model. We
use the `extractModelGenes` function from `glmnet`.

``` r
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)

head(genes.relevant, 10)
#>       Genes Coefficients
#> 1     EDA2R   0.26777221
#> 2    ADRA2B  -0.17828527
#> 3    ZNF229  -0.17157346
#> 4  BAIAP2L2   0.16618988
#> 5  C12orf60  -0.14479918
#> 6     GPR26  -0.12854039
#> 7     POC1A  -0.11831739
#> 8      HAP1   0.10563272
#> 9     FATE1   0.10443456
#> 10   IQGAP3   0.09883007
```

The list of hub genes is found in the “Genes” column of `genes.relevant`
variables. The coefficients are also reported.

Functions employed:
-------------------

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
    -   **return.genes.subset:** (optional) whether to return the hub
        genes along with the `glmnet` model. Recommended to FALSE and
        use `extractModelGenes` from `CovCoExpNets` over the model later
        on. Defaults to FALSE.

-   **extractModelGenes:** extract the relevant predictors from a
    `glmnet` model. It has the following arguments.

    -   **cvfit:** `glmnet` model to extract the relevant predictors and
        their coefficients from.
    -   **genes.freq:** data.frame obtained from `geneFrequency`. If
        provided, it will also return the average coefficient of each
        gene across all glmnet repetitions. Defaults to none.

Multiple conditions
-------------------

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
#> 1 FAM126A   -0.2092029
#> 2  PLAGL1   -0.1728505
#> 3  SLC9A2   -0.1713773
#> 4   RGS19    0.1572788
#> 5   SYT16    0.1494378
#> 
#> [[2]]
#>     Genes Coefficients
#> 1    WWC2   -0.1478869
#> 2   PTPN3   -0.1418295
#> 3  PLAGL1   -0.1243939
#> 4 SLC12A1    0.1237764
#> 5   WNT5B   -0.1225396
```

Predicting the sex
------------------

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
#> 1        TLX2     4.449681
#> 2      SPESP1    -2.476102
#> 3        HAAO     2.348379
#> 4       ACOT4    -2.265550
#> 5      SMIM24     1.909714
#> 6       NLRP2     1.752245
#> 7  CSGALNACT1     1.643682
#> 8       ADAD2     1.606470
#> 9        HPGD     1.566646
#> 10      CHST9     1.331642
```

``` r
stopCluster(cl)
```
