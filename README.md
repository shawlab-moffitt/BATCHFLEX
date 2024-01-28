# BatchFLEX

An R package for correct batch effects in numeric data (typically RNA sequencing data). Companion to the BatchFlex Shiny Application at:

[shiny link]()

To install package:

install from CRAN:

install.packages("devtools")
install.packages("BiocManager")

```
install.packages("devtools")
install.packages("BiocManager")
#stringi
install.package("https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.3/stringi_1.8.3.tgz",
repos = NULL, type = 'source')

#https://bioconductor.org/packages/release/bioc/html/sva.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/sva_3.50.0.tgz",
repos = NULL, type="source")
#https://bioconductor.org/packages/release/bioc/html/limma.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/limma_3.58.1.tgz",
repos = NULL, type = 'source')
#https://bioconductor.org/packages/release/bioc/html/RUVSeq.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/RUVSeq_1.36.0.tgz",
repos = NULL, type = 'source')
#https://github.com/omnideconv/immunedeconv
devtools::install_github("omnideconv/immunedeconv")
#https://bioconductor.org/packages/release/bioc/html/Harman.html
install.packages("https://bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/Harman_1.30.0.tgz",
repos = NULL, type = 'source')
#https://www.bioconductor.org/packages/release/bioc/html/affyPLM.html
install.packages("https://www.bioconductor.org/packages/release/bioc/bin/macosx/big-sur-x86_64/contrib/4.3/affyPLM_1.78.0.tgz",
repos = NULL, type = 'source')

devtools::install_github('shawlab-moffitt/BATCHFLEX')
```

---
title: "Using BatchFlex to Perform Batch Correction"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using BatchFlex to Perform Batch Correction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Background 

Integrating multiple datasets for large-scale genomic or transcriptomic analysis is a difficult and tedious process that typically introduces significant batch effects. Batch effects can arise from many different types of technical variation from differing institutions, personnel, or instruments and often arise even in single studies when sample prep and analysis occurs over multiple days or weeks. Batch effects can significantly confound the data, making detection of true differences difficult. Although many tools are available to correct for these effects, not every method is appropriate in all instances. The BatchFLEX package is a comprehensive tool for simulating, merging, preprocessing, correcting, and evaluating data. BatchFLEX helps researchers determine if batch correction is necessary and helps researches choose the most optimal method of batch correction for their data. Finally, after batch correction, BatchFLEX allows for easy export of a corrected matrix and meta file for downstream analysis. 

## Methods

- Limma (`limma::removeBatchEffect`)
- ComBat (`sva::ComBat`)
- Mean Centering (`bapred::meancenter`)
- ComBatseq (`sva::ComBat_seq`)
- Harman (`Harman::harman`)
- RUVg (`RUVSeq::RUVg`)
- SVA (`sva::sva`)

## Using BatchFlex

```{r load-packages, include = FALSE}
library(BatchFLEX)
library(ggpubr)
```

### BatchFLEX with Simulated data
#### Simulated data used individually
For quick evaluation of the BatchFLEX tool, a `simulate_data` function is included within BatchFLEX. This tool can generate simulated data with multiple treatment effects and with multiple batch effects based on user specified parameters. Then, a meta file is automatically generated based on the batch information. If `simulate_data` is used alone, a list of the matrx and meta files is exported, which can then be used as inputs in `Batch_FLEX`. Defaults for `simulate_data` can be found within the help page by using `?simulate_data`. If `simulate_data` is called within the `Batch_FLEX` function, then the sim_mat and sim_meta files are automatically assigned to the mat and meta file for `Batch_FLEX`. Additionally, batch.1 and variable_of_interest is automatically assigned. In the example below, a single correction method and the rle evaluation method are being used for simplicity. 

The simulated data is based on the linear model framework introduced by Gagnon-Bartsch and Speed.

https://academic.oup.com/biostatistics/article/17/1/16/1744198


```{r}
simulated_data = simulate_data()
head(as.data.frame(simulated_data$sim_matrix))
head(simulated_data$sim_meta)
test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "pca", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", color_by = "batch")

test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["pca"]]
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["pca"]]
```

If desired, users can adjust the default simulation to include additional batches or additional treatments as well as many other parameters to test `BatchFLEX`. 

#### User adjusted simulation
```{r}
simulated_data <- simulate_data(num_samples = 1000,
                          num_genes = 2000,
                          num_control = 1800,
                          num_treatments = 2,
                          treatment_effect_multiplier = c(0, 1),
                          treatment_effect = 5,
                          treatment_sd = 1,
                          additional_batch = 4,
                          batch_proportion = 0.4,
                          batch_effect_multiplier = c(0.5, 1, 2, 4),
                          batch_effect = 5,
                          batch_effect_sd = 3,
                          batch_sample_sd_mean = 0,
                          batch_sample_sd_sd = 2,
                          epsilon_mean = 0,
                          epsilon_sd = 1)

test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "umap", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", variable_of_interest = "treatment")

test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["umap"]]
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["umap"]]

```


#### Simulate data within BatchFLEX
```{r}
test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("simulate_data", "batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "ev", batch.1 = "batch", variable_of_interest = "treatment", color_by = "batch")

test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["ev"]]
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["ev"]]
```


### BatchFLEX with merged data
`BathFLEX` comes preloaded with multiple datasets downloaded from GEO with meta files that have been manually adjusted for optimal merging. These files can be merged using `merge_data` to generate a matrix and meta file similar to the example matrix and meta file included in `BatchFLEX`. `merge_data` can be utilized as a stand alone function as shown below. If a named list is provided, `merge_data` will add a column to the merged_meta file specifying the named study that each sample came from.

#### Merge data used individually
```{r}
matrices_to_be_merged <- list("GSE112876" = BatchFLEX::GSE112876_matrix,
                              "GSE15907" = BatchFLEX::GSE15907_matrix,
                              "GSE37448" = BatchFLEX::GSE37448_matrix,
                              "GSE60336" = BatchFLEX::GSE60336_matrix,
                              "GSE75195" = BatchFLEX::GSE75195_matrix,
                              "GSE75202" = BatchFLEX::GSE75202_matrix,
                              "GSE75203" = BatchFLEX::GSE75203_matrix)
meta_files_to_be_merged <- list(BatchFLEX::GSE112876_meta,
                                BatchFLEX::GSE15907_meta,
                                BatchFLEX::GSE37448_meta,
                                BatchFLEX::GSE60336_meta,
                                BatchFLEX::GSE75195_meta,
                                BatchFLEX::GSE75202_meta,
                                BatchFLEX::GSE75203_meta)
test_merge <- merge_data(merge_matrix_files = matrices_to_be_merged, merge_meta_files = meta_files_to_be_merged)
head(as.data.frame(test_merge$merged_matrix))
head(test_merge$merged_meta)
```

If `merge_data` is supplied with an unnamed list, a column with a generic study name will be added. 

```{r}
matrices_to_be_merged <- list(BatchFLEX::GSE112876_matrix,
                              BatchFLEX::GSE15907_matrix,
                              BatchFLEX::GSE37448_matrix,
                              BatchFLEX::GSE60336_matrix,
                              BatchFLEX::GSE75195_matrix,
                              BatchFLEX::GSE75202_matrix,
                              BatchFLEX::GSE75203_matrix)
meta_files_to_be_merged <- list(BatchFLEX::GSE112876_meta,
                                BatchFLEX::GSE15907_meta,
                                BatchFLEX::GSE37448_meta,
                                BatchFLEX::GSE60336_meta,
                                BatchFLEX::GSE75195_meta,
                                BatchFLEX::GSE75202_meta,
                                BatchFLEX::GSE75203_meta)
test_merge_1 <- merge_data(merge_matrix_files = matrices_to_be_merged, merge_meta_files = meta_files_to_be_merged)
head(as.data.frame(test_merge$merged_matrix))
head(test_merge$merged_meta)
```

`merge_data` can also be used without adding any meta files. When this occurs, a basic meta file will be created specifying the sample names and the corresponding study they were retrieved from. 

```{r}
matrices_to_be_merged <- list("GSE112876" = BatchFLEX::GSE112876_matrix,
                              "GSE15907" = BatchFLEX::GSE15907_matrix,
                              "GSE37448" = BatchFLEX::GSE37448_matrix,
                              "GSE60336" = BatchFLEX::GSE60336_matrix,
                              "GSE75195" = BatchFLEX::GSE75195_matrix,
                              "GSE75202" = BatchFLEX::GSE75202_matrix,
                              "GSE75203" = BatchFLEX::GSE75203_matrix)
test_merge_2 <- merge_data(merge_matrix_files = matrices_to_be_merged)
head(as.data.frame(test_merge$merged_matrix))
head(test_merge$merged_meta)
```
#### Merge data within BatchFLEX

Similar to `simulate_data`, `merge_data` can be used within the `Batch_FLEX` function and automatically pipe the output to the `batch_correct` and `batch_evaluate` functions. 

```{r}
matrices_to_be_merged <- list(BatchFLEX::GSE112876_matrix,
                              BatchFLEX::GSE15907_matrix,
                              BatchFLEX::GSE37448_matrix,
                              BatchFLEX::GSE60336_matrix,
                              BatchFLEX::GSE75195_matrix,
                              BatchFLEX::GSE75202_matrix,
                              BatchFLEX::GSE75203_matrix)
meta_files_to_be_merged <- list(BatchFLEX::GSE112876_meta,
                                BatchFLEX::GSE15907_meta,
                                BatchFLEX::GSE37448_meta,
                                BatchFLEX::GSE60336_meta,
                                BatchFLEX::GSE75195_meta,
                                BatchFLEX::GSE75202_meta,
                                BatchFLEX::GSE75203_meta)
test_merge_BatchFLEX <- Batch_FLEX(
  Batch_FLEX_function = c("merge_data", "preprocess_matrix", "batch_correct", "batch_evaluate"), 
  merge_matrix_files = matrices_to_be_merged, 
  merge_meta_files = meta_files_to_be_merged, 
  batch.1 = "batchflex_study", 
  variable_of_interest = "MajorCellType",
  correction_method = "ComBat",
  evaluation_method = "ev",
  quantnorm = FALSE
  )
head(as.data.frame(test_merge_BatchFLEX[["merge_data"]][["merged_matrix"]]))
head(test_merge_BatchFLEX[["merge_data"]][["merged_meta"]])
head(as.data.frame(test_merge_BatchFLEX[["data_matrices"]][["Unadjusted_Log2_Norm"]]))
head(as.data.frame(test_merge_BatchFLEX[["data_matrices"]][["ComBat"]]))
test_merge_BatchFLEX[["batch_evaluation"]][["Unadjusted_Log2"]][["batch1"]][["Plots"]][["ev"]]
test_merge_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["ev"]]
```


### BatchFLEX with example input data
#### Example Data
`BatchFLEX` contains example data from several GEO studies that can be explored in the `?example_mat` help page. The number of samples was subset to include 48 to keep processing time and data size small. The data object `example_mat` includes the expression data for these 48 samples where the samples are in columns and genes in the rows. Sample level information can be found in `example_meta` such as the cell type the sample is and the study from which it came. 

```{r}
data(example_meta)
data(example_mat)
test_example <- Batch_FLEX(mat = preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE), 
           meta = BatchFLEX::example_meta, 
           batch.1 = "batchflex_study", 
           variable_of_interest = "MajorCellType",
           correction_method = "Limma",
           evaluation_method = "rle")
head(as.data.frame(BatchFLEX::example_mat))
head(as.data.frame(BatchFLEX::example_meta))
```

### Batch correct

#### Batch correct used individually
`BatchFLEX` can use multiple methods of correction either individually or simultaneously. These include "Limma", "ComBat", "ComBatseq", "Mean Centering", "Harman", "RUVg", and "SVA".

##### Limma
Limma can be used with a one batch or two batches and can include the variable of interest in the model matrix. 

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "Limma",
                                   batch.1 = "batchflex_study",
                                   batch.2 = NULL,
                                   variable_of_interest = "MajorCellType")
head(as.data.frame(adjusted_data))
```

##### ComBat
ComBat uses a single batch, but can also be used for non parametric data by changing the par.prior parameter.

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "ComBat",
                                   batch.1 = "batchflex_study",
                                   par.prior = TRUE)
head(as.data.frame(adjusted_data))
```

##### ComBatseq
ComBatseq is unique among the correction methods as it requies raw count data instead of log transformed or normalized data.
```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "ComBatseq",
                                   batch.1 = "batchflex_study",
                                   variable_of_interest = "MajorCellType")
head(as.data.frame(adjusted_data))
```

##### RUVg
RUVg requires a set of negative control genes known apriori not to be differentially expressed between the samples. Alternatively, a set of known positive or negative controls with known expression fold changes can be used and RUVg will be applied to the control-centered log counts. BatchFlex has a few sets of house keeping genes to select from for both mouse and human. Alternatively, another list could be provided by the user. RUVg also includes other technical factors that can be adjusted including k, drop, center, round, epsilon, and tolerance

Mouse:
- `data(eisenberg_mouse)`
- `data(lin500_mouse)`
- `data(hsiao_mouse)`

Human: 
- `data(eisenberg)`
- `data(lin500)`
- `data(hsiao)`


```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "RUVg",
                                   batch.1 = "batchflex_study",
                                   housekeeping = BatchFLEX::hsiao_mouse,
                                   k = 2,
                                   drop = 0,
                                   center = FALSE,
                                   round = FALSE,
                                   tolerance = 1e-8)
head(as.data.frame(adjusted_data))
```

##### SVA
SVA is unique as it does not require batch information and generates surrogate variables that are used for adjustment. The sva_nsv_method parameter is used to choose the method of surrogate variable esimation. Default is set to "leek".

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "SVA",
                                   variable_of_interest = "MajorCellType",
                                   sva_nsv_method = "leek")
head(as.data.frame(adjusted_data))

```

##### Harman and Mean Centering
Harman and Mean Centering do not require special parameters, however, using Harman on large datasets may require significant time. 


```{r}
adjusted_data <- batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat, quantnorm = FALSE),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = c("Harman", "Mean Centering"),
                                   batch.1 = "batchflex_study",
                                   variable_of_interest = "MajorCellType")
head(as.data.frame(adjusted_data[["Harman"]]))
head(as.data.frame(adjusted_data[["Mean Centering"]]))
```

#### Batch correct within Batch_FLEX
`batch_correct` can be used within the `Batch_FLEX` function and automatically pipe the output to the `batch_evaluate` functions. 

```{r}
test_example <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), mat = example_mat, meta = example_meta, correction_method = "Limma", evaluation_method = "pca", batch.1 = "batchflex_study")
head(as.data.frame(test_example$data_matrices$Limma))
test_example$batch_evaluation$Unadjusted$batch1$Plots$pca
test_example$batch_evaluation$Limma$batch1$Plots$pca
```

### Visualize

#### Batch evalute used individually with specific evaluation methods
BatchFLEX utilizes many different methods of evaluation to determine whether batch correction is necessary and to assess whether the correction effectively removed the batch effect. If a variable_of_interest is known then BatchFLEX also ensures it is maintained after correction. `batch_evaluate` can be used as a stand alone function. By default, all evaluation methods are conducted, however, if desired, individual evaluation methods can be chosen. By default, when applicable, a black and white, a colored by batch, and a colored by the variable of interest plot will be generated for all plots. If desired, the color_by parameter can be adjusted to the desired plot or plots of interest.


##### Principal component analysis
By default, pca will include all annotation types, which includes annotation by cluster and annotation by the meta file. This can be adjusted by changing the annotation parameter. If a cluster number is not provided, then BatchFLEX will use silhouette analysis to generate an ideal cluster number for each matrix file. 
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "pca", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", annotation = "meta")
test_evaluate$Plots$pca
```

##### Cluster analysis
By default, the cluster analysis will include an elbow plot, a silhouette plot and a dunn index plot. If all three are not desired, then the cluster_analysis_method parameter can be changed. 
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "cluster_analysis", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", cluster_analysis_method = "silhouette")
test_evaluate$Plots$cluster_analysis
```

##### PCA details
The PCA details analysis includes many helpful plots and matrices that can be used to assess the PCA. By default, `batch_evaluate` will use the batch.1 parameter when generating these files, however, this can be adjusted by altering the pca_factors parameter

```{r}
test_evaluate <- batch_evaluate(evaluation_method = "pca_details", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", pca_factors = "MajorCellType")
test_evaluate$Plots$pca_details
```


##### Multiple components PCA
By default, the multiple components plot will include the top 5 principal components. This can be changes by altering the ncomponents parameter.
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "mc_pca", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", color_by = "batch", ncomponents = 10)
test_evaluate$Plots$mc_pca
```

##### Relative log expression analysis
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "rle", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", color_by = "variable_of_interest")
test_evaluate$Plots$rle
```

##### Explanatory variables plot
By default, the explanatory variables plot will include the batch.1 and the variable of interest groups. If another column from the meta file is desired, then the variable_choices parameter can be altered.
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "ev", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", variable_choices = c("batchflex_study", "Major Lineage"))
test_evaluate$Plots$ev
```

##### Surrogate variable analysis
By default, SVA will use the "be" method for surrogate variable estimation. If desired, the "leek" method can be chosen by chaning the sva_nsv_method.
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "sva", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", sva_nsv_method = "be")
test_evaluate$Plots$sva
```

##### Unform manifold approximation and projection
```{r}
test_evaluate <- batch_evaluate(evaluation_method = "umap", mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType")
test_evaluate$Plots$umap
```

#### Batch evaluate used individually with all evaluation methods
By default, `batch_evaluate` will include all evaluation methods that are avaiable based on the input criteria. 

```{r}
test_evaluate <- batch_evaluate(mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType")
test_evaluate$Plots
```

#### Batch evaluate used within BatchFLEX
`batch_evaluate` can also be used in conjunction with other functions within the top level `Batch_FLEX` function.

```{r}
test_evaluate <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "Limma", evaluation_method = "ev", mat = example_mat, meta = example_meta, batch.1 =  "batchflex_study", variable_of_interest = "MajorCellType")
test_evaluate$Matrices
test_evaluate$Plots
```

### BatchFlex easy usage
By default, `Batch_FLEX` will run `batch_correct` and `batch_evaluate`. It will correct using all avaialable methods except for ComBatseq, which requires a counts based matrix file. It will also use a self contained housekeeping gene list for the RUVg evaluation. By default, this is set to human, but can be changed to mouse or a user provided housekeeping gene list can be used. All available evaluation methods will be used unless a variable of interest is not provided, which eliminates the use of sva. The output is an organized list tree separated into the unadjusted and adjusted expression matrices and the the evaluation plots and matrices. The evaluation list is organized according to the  correction method. Each evaluation list contains a list of evaluation matrices and a list of evaluation plots. An example of how to access to the plots is shown below. 
```{r}
example_mat_abbreviated <- example_mat[,1:100]
example_meta_abbreviated <- example_meta[1:100,]
test_Batch_FLEX <- Batch_FLEX(mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType")

head(as.data.frame(test_Batch_FLEX[["data_matrices"]][["Unadjusted"]]))
head(as.data.frame(test_Batch_FLEX[["data_matrices"]][["Limma"]]))

test_Batch_FLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]]
test_Batch_FLEX[["batch_evaluation"]][["Limma"]][["batch1"]][["Plots"]]
```


`ggpubr` has a functioned called `ggarrange` that allows for combining plots into plates, so can visualize both the unadjusted and the Limma adjusted expression together. This allows us to see the impact that 1) the study has on the expression profiles and 2) the 

```{r, fig.width = 20}
example_mat_abbreviated <- example_mat[,1:100]
example_meta_abbreviated <- example_meta[1:100,]
test_Batch_FLEX <- Batch_FLEX(mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType", evaluation_method = "pca")

ggpubr::ggarrange(
  test_Batch_FLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["Limma"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]], 
  test_Batch_FLEX[["batch_evaluation"]][["Harman"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["RUVg"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["SVA"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  common.legend = TRUE)
```

Finally, `BatchFLEX` has the ability to generate `ggarrange` plots for all evaluation methods and to export these plots along with all other plots and matrices to a BatchFLEX_analysis folder by evoking the `BatchFLEX_export` function.

```{r}
example_mat_abbreviated <- example_mat[,1:100]
example_meta_abbreviated <- example_meta[1:100,]
test_Batch_FLEX <- Batch_FLEX(Batch_FLEX_function = c("batch_evaluate", "batch_correct", "BatchFLEX_export"), correction_method = c("Limma", "ComBat"), evaluation_method = c("pca","mc_pca"), mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "MajorCellType")
devtools::session_info()
```





# Disclamer

Copyright 2024 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.








