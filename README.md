
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SuCoNets

<!-- badges: start -->

<!-- badges: end -->

Nowadays, the way to create a co-expression network is to use
Hierarchical Clustering. This package allows you to create a
co-expression network using the glmnet algorithm. In addition, the
network created will be relative to a covariate of the sample to be
studied, so we are creating a supervised coexpression network.

## Installation

You can install foofactors like so:

``` r
install.packages("SuCoNets")
```

## Example

Suppose we have an expression matrix, data, where the columns are blood
samples and the rows are genes, so that each sample is identified by the
numerical values taken by the genes. Suppose we also have a numeric
vector, age, with the age of each individual to whom each blood sample
corresponds.

An example of a typical execution of the functions contained in this
package would be as follows:

``` r
library(SuCoNets)
age <- normalize(age)
data <- scn(data)

data <- rRedundantPredictors(data)
genes.seleccionados <- detectGenes(data,age)

all.genes <- coexpressionNetwork(data, genes.seleccionados)
```
