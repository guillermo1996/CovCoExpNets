-   [Data extraction from the GTEx
    portal](#data-extraction-from-the-gtex-portal)
    -   [Requirements](#requirements)
-   [Steps](#steps)
    -   [Step 1: Data download](#step-1-data-download)
    -   [Step 2: Extract brain samples by
        tissues](#step-2-extract-brain-samples-by-tissues)
    -   [Step 3: Age & Sex extraction](#step-3-age-sex-extraction)
    -   [Step 4: Preprocessing](#step-4-preprocessing)
-   [Results](#results)
-   [Additional Notes](#additional-notes)
    -   [List of conditions](#list-of-conditions)

Data extraction from the GTEx portal
====================================

Requirements
------------

The requirements for this tutorial are the following libraries:

``` r
library(CovCoExpNets)
library(data.table)
library(magrittr)
library(logger)

#doParallel::registerDoParallel(13)
cl <- makeCluster(13)
doParallel::registerDoParallel(cl)
```

The package `data.table` is required to read the initial data files from
the GTEx, as well as for some data.frame manipulation later on. The
package `magrittr` allows to use of the pipe-like operator `%>%`, which
is employed across all tutorials and this package for its convenience.
In addition, the `CovCoExpNets` package uses the package `doParallel`
and `foreach` to execute as many steps in parallel as possible. Since we
work with different brain tissues and they do not interact with each
other, it is possible to parallelize most parts of this tutorial.
However, we need to specify the number of cores available.

Later on, we will also use the package `caret` to split the dataset by
tissues, but it is only required if you are interested in more than one
tissue or condition.

Steps
=====

Step 1: Data download
---------------------

The first step is to download the required data files. All the data for
this tutorial can be obtained from the [GTEx
Portal](https://www.gtexportal.org/home/datasets). Specifically, we need
three different files:

-   The subject information
    (GTEx\_Analysis\_v8\_Annotations\_SubjectPhenotypesDS.txt): it
    contains the information about the donors, in particular, the
    approximate age and the sex.

-   The sample information
    (GTEx\_Analysis\_v8\_Annotations\_SampleAttributesDS.txt): it
    contains the information of the samples. There are many fields in
    this document, but we are only interested in `SAMPID` with the ID of
    the sample, and both `SMTS` and `SMTSD` to indicate the tissue from
    which the sample was extracted.

-   The sample data
    (GTEx\_Analysis\_2017-06-05\_v8\_RNASeQCv1.1.9\_gene\_tpm.gct.gz):
    the file containing the TPM values. It consists in genes as rows and
    samples’s ID as columns.

``` r
# Download links
subject_info = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
sample_info = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

data.path = "~/GTEx_data/"
dir.create(data.path)
#> Warning in dir.create(data.path): '/home/grocamora/GTEx_data' ya existe

download.file(subject_info, paste0(data.path, "subject_info.txt"))
download.file(sample_info, paste0(data.path, "sample_info.txt"))

# The sample data can be downloaded from here: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
# It must be extracted in the same folder as the other files.
```

Step 2: Extract brain samples by tissues
----------------------------------------

The step to extract the brain data by tissues is a bit complex and can
be split in different steps:

1.  Extract the data samples’ ID for the brain tissues.

2.  Extract the TPM values for the samples extracted from the brain.

3.  Split the previous datasets by tissue using a dummy matrix.

``` r
# Extract the data samples from the brain
tissue.selection <- "Brain"
df.samples_data <- fread(paste0(data.path, "sample_info.txt"))
df.samples_data <- df.samples_data[df.samples_data$SMTS == tissue.selection] 

# Get the TPM data only for brain tissues
df.samples <- fread(paste0(data.path, "samples.gct"))
match_samples = intersect(colnames(df.samples), df.samples_data$SAMPID)
data <- df.samples[, ..match_samples]
data <- as.matrix(data, rownames = df.samples$Name)
  
rm(df.samples)
gc()
#>              used   (Mb) gc trigger    (Mb)   max used    (Mb)
#> Ncells    1806926   96.6    3041455   162.5    3041455   162.5
#> Vcells 1091354866 8326.4 2662336689 20312.1 2216893795 16913.6

# Construct a dummy matrix with samples as rows and tissues as columns.
# If a sample is from a given tissue, it will have a 1 in the corresponding
# cell.
library(caret)
#> Loading required package: ggplot2
#> Loading required package: lattice
dmy.matrix <- caret::dummyVars("SAMPID ~ SMTSD", data = df.samples_data)
df.bool_data <- data.frame(predict(dmy.matrix, newdata = df.samples_data))
condition.names <- sort(stringr::str_split(unique(df.samples_data$SMTSD), " - ", simplify = TRUE)[, 2])
colnames(df.bool_data) <- condition.names
rownames(df.bool_data) <- df.samples_data$SAMPID
df.bool_data <- df.bool_data[match_samples, ] 

# Split the data by tissues
data <- CovCoExpNets::splitByCondition(data, df.bool_data)
rm(df.samples_data)
gc() 
#>              used   (Mb) gc trigger    (Mb)   max used    (Mb)
#> Ncells    2324646  124.2    4675080   249.7    3041455   162.5
#> Vcells 1092760004 8337.1 2662336689 20312.1 2216893795 16913.6
```

In the end, we should have a `data` variable with a list containing the
data.frames for each tissue. Each element of the list corresponds to a
specific tissue, and they all have as rows the genes (56200) and the
samples’ ID as columns (from 139 to 255, depending on the tissue).

Step 3: Age & Sex extraction
----------------------------

The next step is to extract the age and sex of the samples’ donors for
each tissue. To do so, we use the `generateCovariate` function defined
next. It looks for the column names of the data matrix (containing the
samples’ ID) and matches them to the donors’ age and sex in the subject
information file:

``` r
generateCovariate <- function(data, df.subjects, covariate = "AGE"){
  return.list = is(data, "list")
  if(!return.list) data <- list(data)
  
  covariate.combined <- foreach(i = 1:length(data), .combine = "c") %dopar%{
    data.i = data[[i]]
    sample_donnors = c()

    for(name in colnames(data.i)){
      split_id = stringr::str_split(name, "-", simplify = TRUE)
      donnor_id = paste0(split_id[1], "-", split_id[2])
      sample_donnors = c(sample_donnors, donnor_id)
    }

    covariate.tissue = merge(sample_donnors, df.subjects, by.x="x", by.y="SUBJID")[, covariate]
    if(is(covariate.tissue, "character") & covariate == "AGE"){
      covariate.tissue = unname(sapply(covariate.tissue, function(x) as.numeric(substr(x, 0, 2)) + 5))
    }
    CovCoExpNets::convertToNamedList(covariate.tissue, names(data), i, return.list)
  }

  CovCoExpNets::returnList(return.list, covariate.combined)
}

df.subjects <- fread(paste0(data.path, "subject_info.txt"))
age <- generateCovariate(data, df.subjects, "AGE")
sex <- generateCovariate(data, df.subjects, "SEX")
  
rm(df.subjects)
gc()
#>              used   (Mb) gc trigger    (Mb)   max used    (Mb)
#> Ncells    2325247  124.2    4675080   249.7    3041455   162.5
#> Vcells 1092771695 8337.2 2662336689 20312.1 2216893795 16913.6
```

At the end of this step, we should have three variables: `data`, `age`
and `sex`. Each one is a list with an element for each tissue. The
`data` variables consists in the TPM values of the samples for that
given tissue, while the `age` and `sex` variables contain the biological
age and sex of the subjects for that given tissue.

> :warning: **Note**: Public available data from the GTEx project only
> provides the age of the subjects in 10 years ranges. For this
> tutorial, we set the age of the subject as the middle of the range.
> For example, if the range is 60-69, the age estimated is 65.

Step 4: Preprocessing
---------------------

The last step is to preprocess the data. The whole preprocessing
pipeline is the following:

1.  Execute Log2(TPM + 1) to reduce the variability of the data and to
    better approximate its distribution to a normal-like.
2.  Execute the activation threshold, in which we require for at least
    80% of the samples to be higher than 0.1.
3.  Eliminate the genes with a low variation.
4.  Require the genes to be Protein Coding. Previous to this step, we
    modified the GTEx naming convention from the GTEx project to match
    the ENSEMBL notation and we translate them into HGCN notation. In
    this step, we can also force the genes to be from autosomal genes.
5.  Z-Score normalization so that every gene follows a gaussian with
    mean 0 and standard deviation 1, so that every gene can be compared
    independent from the scale.

All of these steps are executed in the `dataPreprocess` function
included in `CovCoExpNets`, where every step can be individually
activated or deactivated. We also have to generate a dataset contanining
only non-sexual chromosomes’ genes. This is becasue when trying to
predict the sex of the subject, a single gene from a sexual chromosome
can almost perfectly identify it.

``` r
if(!file.exists(paste0(data.path, "data.autosomes.combined.rds"))){
  data.autosomes = CovCoExpNets::dataPreprocess(data, includeSexChromosomes = F)
  saveRDS(data.autosomes, paste0(data.path, "data.autosomes.combined.rds"))
}

if(!file.exists(paste0(data.path, "data.combined.rds"))){
  data = CovCoExpNets::dataPreprocess(data)
  saveRDS(data, paste0(data.path, "data.combined.rds"))
}
```

To preprocess the age vector, we execute the `normalize` function in
`CovCoExpNets`. It will return a list with three items: mean,
standard.deviation and covariate. The mean and the standard deviation
are stored so that we can recover the real age value to better interpret
the results

``` r
age = CovCoExpNets::normalize(age)
saveRDS(age, paste0(data.path, "age.combined_approx.rds"))
```

The sex vector does not need any normalization since it is only a binary
variable. However, it is recommended to set it to 0 and 1, instead of 1
and 2:

``` r
sex = sapply(sex, function(x) factor(x-1, levels = c(0, 1), labels=c("Male", "Female")))
saveRDS(sex, paste0(data.path, "sex.combined.rds"))
```

Results
=======

By the end of this tutorial, we would have generated the preprocessed
data matrices and ages of the subjects for all brain tissues.

``` r
data = readRDS(paste0(data.path, "data.combined.rds"))
data.autosomes = readRDS(paste0(data.path, "data.autosomes.combined.rds"))
age = readRDS(paste0(data.path, "age.combined.rds"))
sex = readRDS(paste0(data.path, "sex.combined.rds"))

# Number of tissues:
length(data)
#> [1] 13
length(age)
#> [1] 13
length(sex)
#> [1] 13

# Size of the datasets
tissue.samples = lapply(data, function(x) dim(x)[2])
tissue.genes = lapply(data, function(x) dim(x)[1])
tissue.genes.autosomes = lapply(data.autosomes, function(x) dim(x)[1])

df = data.table::rbindlist(list(Samples = tissue.samples, Genes = tissue.genes, Autosomal_genes = tissue.genes.autosomes)) %>% t %>% as.data.frame()
colnames(df) = c("Samples", "Genes", "Autosomal genes")
df
#>                                   Samples Genes Autosomal genes
#> Amygdala                              152 15121           14565
#> Anterior cingulate cortex (BA24)      176 15282           14723
#> Caudate (basal ganglia)               246 15276           14714
#> Cerebellar Hemisphere                 215 15301           14750
#> Cerebellum                            241 15484           14924
#> Cortex                                255 15516           14951
#> Frontal Cortex (BA9)                  209 15442           14875
#> Hippocampus                           197 15161           14607
#> Hypothalamus                          202 15500           14927
#> Nucleus accumbens (basal ganglia)     246 15264           14698
#> Putamen (basal ganglia)               205 15015           14462
#> Spinal cord (cervical c-1)            159 15210           14655
#> Substantia nigra                      139 15086           14528
```

Additional Notes
================

List of conditions
------------------

All the CovCoExpNets functions are prepared to deal with either a unique
condition or multiple.

In this tutorial, we worked with 13 different brain tissues at the same
time, executing their preprocessing in parallel. The functions always
expects named lists as input functions, however, they are also designed
to work with a single covariate. For example, our `data` variable
consists in a named list of 13 data matrices. When we execute the
preprocessing step, the function recognize it as a list, and executes
the preprocessing for every different element of the list. If we were to
provide one single data matrix instead, it will also identify that the
variable is not a list and will execute the pipeline only for that
specific nameless data matrix.

``` r
stopCluster(cl)
```
