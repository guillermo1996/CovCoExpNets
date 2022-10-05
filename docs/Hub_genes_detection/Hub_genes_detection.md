Hub Gene detection algorithm
============================

- [Hub Gene detection algorithm](#hub-gene-detection-algorithm)
  - [Requisites](#requisites)
  - [Summary](#summary)
  - [Steps](#steps)
    - [Step 1: Executing the `glmnet` repetitions](#step-1-executing-the-glmnet-repetitions)
    - [Step 2: Minimum relative ratio of appearance treshold](#step-2-minimum-relative-ratio-of-appearance-treshold)
    - [Step 3: Model generation](#step-3-model-generation)
    - [Step 4: Genes and coefficients extraction](#step-4-genes-and-coefficients-extraction)
  - [Functions employed:](#functions-employed)
  - [Multiple conditions](#multiple-conditions)


In this tutorial, we will show the pipeline to extract the hub genes (or
predictors) for a given covariate. These genes will be used to generate
the co-expression networks.

Requisites
----------

The libraries required for this tutorial are `CovCoExpNets`, `magrittr`
for the `%>%` pipe-like operator and `dplyr`. Additionally,
`CovCoExpNets` uses `foreach` and `doParallel` to execute the functions
in parallel.

``` r
library(CovCoExpNets)
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
doParallel::registerDoParallel(13)
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
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9)
cvfit = CovCoExpNets::glmnetGenesSubset(data, age, genes.subset, seed = 1796)
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)
```

The minimum requirements to execute the pipeline are the dataset and the
covariate vector, the *t* hyperparameter (by default 10) and the *mrfa*
hyperparameter (by defualt 0.9). The output is a list of the hub genes:

``` r
# First 25 relevant genes for Cortex tissue:
genes.relevant$Genes[1:25]
#>  [1] EDA2R      ADRA2B     BAIAP2L2   ZNF229     C12orf60   GPR26     
#>  [7] POC1A      IQGAP3     HAP1       FATE1      MCF2L2     IGF1      
#> [13] TLX2       ST6GALNAC2 ABCB4      LBHD1      ETV4       SCUBE1    
#> [19] SOX7       SYT6       QPRT       ZP3        ZKSCAN8P1  SFTPA2    
#> [25] RPSA      
#> 204 Levels: AARD ABCB4 ADAD2 ADAMTS7 ADRA2B AHRR ALDH3A1 ALOX15 ALOX15B ... ZP3
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
#> 1    EDA2R    0.2160115    1
#> 2    GPR26   -0.2121507    1
#> 3   ZNF229   -0.1527472    1
#> 4   ADRA2B   -0.1514803    1
#> 5 BAIAP2L2    0.1183130    1
tail(genes.freq, 5)
#>       Genes  Coefficients iter
#> 2288 PCDHB8  4.974914e-05   10
#> 2289    LPA -9.606374e-06   10
#> 2290   LSM6  7.456103e-06   10
#> 2291  KYAT3 -1.796169e-06   10
#> 2292  TNNT2  1.633166e-06   10
```

### Step 2: Minimum relative ratio of appearance treshold

In this step, we reduce the list of predictors returned in each
iteration. To do so, we use the hyperparameter *mrfa* to specify the
minimum relative frequency in which a predictor must appear. To do so,
we use the function `reduceGenes` from `CovCoExpNets`, where we count
how many times each gene is selected as relevant across the `glmnet`
repetitions and divide it by the total number of repetitions so that we
have the percentage of repetitions in which it was selected as relevant.
We require a minimum appearance threhsold.

By default, we recommend an *mrfa* of 0.9:

``` r
genes.subset = CovCoExpNets::reduceGenes(genes.freq, mrfa = 0.9)

# First 25 genes:
genes.subset[1:25]
#>  [1] "AARD"      "ABCB4"     "ADAD2"     "ADRA2B"    "AHRR"      "ALDH3A1"  
#>  [7] "ALOX15B"   "ALPL"      "AMZ1"      "ANKRD18B"  "APC"       "ARHGAP11A"
#> [13] "ASCL2"     "ATAD3C"    "BAIAP2L2"  "C12orf60"  "C1QL4"     "CARD16"   
#> [19] "CBY3"      "CCDC163"   "CCDC170"   "CCDC81"    "CD163L1"   "CD3E"     
#> [25] "CDK11A"
```

The output will be a set of genes that pass the threshold. In total,
there were 277 different genes across all repetitions, while 201 of
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
#> min 0.0005074    74 0.03150 0.002368     195
#> 1se 0.0008079    69 0.03368 0.002410     193
```

### Step 4: Genes and coefficients extraction

The final step is to extract the hub genes from the generated model. We
use the `extractModelGenes` function from `glmnet`.

``` r
genes.relevant = CovCoExpNets::extractModelGenes(cvfit)

head(genes.relevant, 10)
#>       Genes Coefficients
#> 1     EDA2R   0.27058261
#> 2    ADRA2B  -0.17935231
#> 3  BAIAP2L2   0.16943714
#> 4    ZNF229  -0.16852248
#> 5  C12orf60  -0.15005711
#> 6     GPR26  -0.12395506
#> 7     POC1A  -0.11090052
#> 8     FATE1   0.10987335
#> 9      HAP1   0.10321477
#> 10   IQGAP3   0.09903114
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
        reptition. Different solutions will not only have different
        returned predictors, but also different performance. Use this
        setting if you wish the prioritize better RMSE iterations when
        selecting the pruned predictors.
    -   **data.test.extra:** (optional) numeric data matrix with genes
        as rows and samples as columns. If the train/test split was
        executed before this step and you want to estimate the RMSE over
        the test dataset, provide the test dataset in here. By default,
        none.
    -   **covariate.test.extra:** (optional) numeric vector with the age
        of the sample’s donors. If the train/test split was executed
        before this step and you want to estimate the RMSE over the test
        dataset, provide the test dataset in here. By defualt, none.
    -   **sample.prob:** (optional) when executing the train/test split,
        it dictates the probability of each sample being executed. See
        [sample](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sample)
        for more information.
    -   **seed:** (optional) numeric seed to ensure reproducibility.
        Because of parallelization, even with an input seed results migh
        differ.
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
        Because of parallelization, even with an input seed results migh
        differ.
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
#> $`Cerebellar Hemisphere`
#>     Genes Coefficients
#> 1 FAM126A   -0.2047052
#> 2   SYT16    0.1718157
#> 3  SLC9A2   -0.1713823
#> 4   RGS19    0.1630940
#> 5  PLAGL1   -0.1605985
#> 
#> $Cerebellum
#>     Genes Coefficients
#> 1    WWC2   -0.1579401
#> 2   PTPN3   -0.1411976
#> 3  PLAGL1   -0.1276331
#> 4   WNT5B   -0.1257876
#> 5 SLC12A1    0.1145380
```
