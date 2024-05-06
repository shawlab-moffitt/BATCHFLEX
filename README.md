# BatchFLEX

An R package for correct batch effects in numeric data (typically RNA sequencing data). Companion to the BatchFlex Shiny Application at:

* [Shiny GitHub Package](https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp)
* [Live Shiny IO](https://shawlab-moffitt.shinyapps.io/batch_flex/)

To install package:

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("shawlab-moffitt/BATCHFLEX", dependencies = TRUE)
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
### BatchFlex easy usage
By default, `Batch_FLEX` will run `batch_correct` and `batch_evaluate`. It will correct using all available methods except ComBatseq, which requires a counts based matrix file. For simplicity, Harman was removed as it requires significant time for large datasets, however, in general, `Batch_FLEX` can be ran without indicating a correction_method. `Batch_FLEX` also uses a self contained housekeeping gene list for the RUVg evaluation. By default, this is set to human, but can be changed to mouse or a user provided housekeeping gene list can be used. All available evaluation methods will be used unless a variable of interest is not provided, which eliminates the use of sva. The output is an organized list tree separated into the unadjusted and adjusted expression matrices and the the evaluation plots and matrices. The evaluation list is organized according to the  correction method. Each evaluation list contains a list of evaluation matrices and a list of evaluation plots. An example of how to access to the plots is shown below. 
```{r, fig.width = 10, fig.height = 10, results = 'hide'}
#test_Batch_FLEX <- Batch_FLEX(mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")

test_Batch_FLEX <- Batch_FLEX(mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", correction_method = c("Limma", "ComBat", "Mean Centering", "RUVg", "SVA"))

test_Batch_FLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]]
test_Batch_FLEX[["batch_evaluation"]][["Limma_adjusted"]][["batch1"]][["Plots"]]
```

`BatchFLEX` includes a small function to run the Shiny app locally by downloading and temporarily storing the tar.gz file from the BATCH-FLEX-shinyApp GitHub page. This can be run alone or within the overall BatchFLEX function. 
```{r, eval = FALSE}
BatchFLEX_shiny()
test_Batch_FLEX <- Batch_FLEX(Batch_FLEX_function = c("preprocess_matrix", "batch_evaluate", "batch_correct", "BatchFLEX_export", "BatchFLEX_shiny"), correction_method = c("ComBat", "Limma", "Harman", "SVA", "RUVg"), mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")
```

`ggpubr` has a function called `ggarrange` that allows for combining plots into plates to visualize all correction methods simultaneously. This allows the user to easily evaluate the correction methods and to identify the most appropriate correction method for their dataset.

```{r, fig.width = 10, fig.height = 10, results = 'hide'}
example_mat_abbreviated <- BatchFLEX::preprocess_matrix(BatchFLEX::example_mat)[,1:100]
example_meta_abbreviated <- example_meta[1:100,]
test_Batch_FLEX <- Batch_FLEX(mat = example_mat_abbreviated, meta = example_meta_abbreviated, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", evaluation_method = "pca")

ggpubr::ggarrange(
  test_Batch_FLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["Limma_adjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["ComBat_adjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]], 
  test_Batch_FLEX[["batch_evaluation"]][["Harman_adjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["RUVg_adjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  test_Batch_FLEX[["batch_evaluation"]][["SVA_adjusted"]][["batch1"]][["Plots"]][["pca"]][["pca_meta_batch_colored_pca"]],
  common.legend = TRUE)
```

`BatchFLEX` has the ability to generate `ggarrange` plots for all evaluation methods and to export these plots along with all other plots and matrices to a BatchFLEX_analysis folder by evoking the `BatchFLEX_export` function. Because of the size, multiple heatmaps comparing the uncorrected matrix to each corrected matrix are generated instead of a single plot comparing all methods simultaneously. Also, because of the significant number of matrices and plots generated by `BatchFLEX`, parallel functionality was added to help improve computational time. This is especially important when using larger datasets and when comparing multiple correction methods. By default, the SOCK method is used for windows and the FORK method is used for mac and linux. Additionally, the cores setting is defaulted to use one worker per matrix for the evaluation step if allowed by the computer. Because using parallel functionality is significant impacted by the users computer configuration, this setting is set to FALSE by default. 

```{r, fig.width = 10, fig.height = 10, eval = FALSE}
test_Batch_FLEX <- Batch_FLEX(Batch_FLEX_function = c("preprocess_matrix", "batch_evaluate", "batch_correct", "BatchFLEX_export"), correction_method = c("ComBat", "Limma", "Harman", "SVA", "RUVg"), mat = example_mat, meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", parallelize = TRUE)

```

The remaining focus of the tutorial is used to highlight each function employed by `BatchFLEX` in greater detail. 

### BatchFLEX with Simulated data
#### Simulated data used individually
For quick evaluation of the BatchFLEX tool, a `simulate_data` function is included within BatchFLEX. This tool can generate simulated data with multiple treatment effects and with multiple batch effects based on user specified parameters. Then, a meta file is automatically generated based on the batch information. If `simulate_data` is used alone, a list of the matrx and meta files is exported, which can then be used as inputs in `Batch_FLEX`. Defaults for `simulate_data` can be found within the help page by using `?simulate_data`. If `simulate_data` is called within the `Batch_FLEX` function, then the sim_mat and sim_meta files are automatically assigned to the mat and meta file for `Batch_FLEX`. Additionally, batch.1 and variable_of_interest is automatically assigned. In the example below, a single correction method and the rle evaluation method are being used for simplicity. 

The simulated data is based on the linear model framework introduced by Gagnon-Bartsch and Speed.

https://academic.oup.com/biostatistics/article/17/1/16/1744198


```{r, fig.width = 10, fig.height = 10}
simulated_data = simulate_data()
head(as.data.frame(simulated_data$sim_matrix), n = c(5,5))
head(simulated_data$sim_meta, n = c(5,5))
test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "pca", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", color_by = "batch")

test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["pca"]]
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat_adjusted"]][["batch1"]][["Plots"]][["pca"]]
```

If desired, users can adjust the default simulation to include additional batches or additional treatments as well as many other parameters to test `BatchFLEX`. 

#### User adjusted simulation
```{r, fig.width = 10, fig.height = 10}
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
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat_adjusted"]][["batch1"]][["Plots"]][["umap"]]

```


#### Simulate data within BatchFLEX
```{r, fig.width = 10, fig.height = 10}
test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("simulate_data", "batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "ev", batch.1 = "batch", variable_of_interest = "treatment", color_by = "batch")

test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["ev"]]
test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat_adjusted"]][["batch1"]][["Plots"]][["ev"]]
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
head(as.data.frame(test_merge$merged_matrix), n = c(5,5))
head(test_merge$merged_meta, n = c(5,5))
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
head(as.data.frame(test_merge$merged_matrix), n = c(5,5))
head(test_merge$merged_meta, n = c(5,5))
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
head(as.data.frame(test_merge$merged_matrix), n = c(5,5))
head(test_merge$merged_meta, n = c(5,5))
```
#### Merge data within BatchFLEX

Similar to `simulate_data`, `merge_data` can be used within the `Batch_FLEX` function and automatically pipe the output to the `batch_correct` and `batch_evaluate` functions. In this step `preprocess_matrix` is also used to log2 and quantnorm the input matrix. `preprocess_matrix` can also handle raw counts and can remove duplicates.

```{r, fig.width = 10, fig.height = 10}
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
  variable_of_interest = "Major_Lineage",
  correction_method = "ComBat",
  evaluation_method = "ev",
  quantnorm = FALSE
  )
head(as.data.frame(test_merge_BatchFLEX[["merge_data"]][["merged_matrix"]]), n = c(5,5))
head(test_merge_BatchFLEX[["merge_data"]][["merged_meta"]], n = c(5,5))
head(as.data.frame(test_merge_BatchFLEX[["data_matrices"]][["Unadjusted_Log2_Norm"]]), n = c(5,5))
head(as.data.frame(test_merge_BatchFLEX[["data_matrices"]][["ComBat"]]), n = c(5,5))
test_merge_BatchFLEX[["batch_evaluation"]][["Unadjusted_Log2"]][["batch1"]][["Plots"]][["ev"]]
test_merge_BatchFLEX[["batch_evaluation"]][["ComBat_adjusted"]][["batch1"]][["Plots"]][["ev"]]
```


### BatchFLEX with example input data
#### Example Data
`BatchFLEX` contains example data from several GEO studies that can be explored in the `?example_mat` help page. The number of samples was subset to include 48 to keep processing time and data size small. The data object `example_mat` includes the expression data for these 48 samples where the samples are in columns and genes in the rows. Sample level information can be found in `example_meta` such as the cell type the sample is and the study from which it came. 

```{r, fig.width = 10, fig.height = 10}
data(example_meta)
data(example_mat)
test_example <- Batch_FLEX(mat = preprocess_matrix(BatchFLEX::example_mat), 
           meta = BatchFLEX::example_meta, 
           batch.1 = "batchflex_study", 
           variable_of_interest = "Major_Lineage",
           correction_method = "Limma",
           evaluation_method = "umap")
head(as.data.frame(BatchFLEX::example_mat), n = c(5,5))
head(as.data.frame(BatchFLEX::example_meta), n = c(5,5))
test_example$batch_evaluation$Unadjusted$batch1$Plots$umap
test_example$batch_evaluation$Limma_adjusted$batch1$Plots$umap
```
#### Example File
Here's an example of reading in a matrix and meta file for analysis.
```
library(BatchFLEX)
meta <- read.delim("2809_BladderCancerSamples_JoshTimMetaShort_20240430.tsv", fill = TRUE)
matrix <- read.delim("merged_matrix2024_04_25_tsfilterd_v2.tsv")
matrix_numeric <- as.matrix(matrix[,-1])
rownames(matrix_numeric) <- matrix[,1]
test_bladder <- Batch_FLEX(Batch_FLEX_function = c("preprocess_matrix", "batch_correct"), 
           mat = matrix_numeric, 
           meta = meta, 
           batch.1 = "batchflex_study",
           correction_method = "Limma",
           log2 = FALSE,
           quantnorm = TRUE,
           remove.duplicates = FALSE)
test_preprocess <- preprocess_matrix(matrix_numeric, log2 = FALSE, quantnorm = TRUE)
```
### Batch correct

#### Batch correct used individually
`BatchFLEX` can use multiple methods of correction either individually or simultaneously. These include "Limma", "ComBat", "ComBatseq", "Mean Centering", "Harman", "RUVg", and "SVA".

##### Limma
Limma can be used with a one batch or two batches and can include the variable of interest in the model matrix. 

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "Limma",
                                   batch.1 = "batchflex_study",
                                   batch.2 = NULL,
                                   variable_of_interest = "Major_Lineage")
head(as.data.frame(adjusted_data), n = c(5,5))
```

##### ComBat
ComBat uses a single batch, but can also be used for non parametric data by changing the par.prior parameter.

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "ComBat",
                                   batch.1 = "batchflex_study",
                                   par.prior = TRUE)
head(as.data.frame(adjusted_data), n = c(5,5))
```

##### ComBatseq
ComBatseq is unique among the correction methods as it requies raw count data instead of log transformed or normalized data.
```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "ComBatseq",
                                   batch.1 = "batchflex_study",
                                   variable_of_interest = "Major_Lineage")
head(as.data.frame(adjusted_data), n = c(5,5))
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
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "RUVg",
                                   batch.1 = "batchflex_study",
                                   housekeeping = BatchFLEX::hsiao_mouse,
                                   k = 2,
                                   drop = 0,
                                   center = FALSE,
                                   round = FALSE,
                                   tolerance = 1e-8)
head(as.data.frame(adjusted_data), n = c(5,5))
```

##### SVA
SVA is unique as it does not require batch information and generates surrogate variables that are used for adjustment. The sva_nsv_method parameter is used to choose the method of surrogate variable esimation. Default is set to "be".

```{r}
adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = "SVA",
                                   variable_of_interest = "Major_Lineage",
                                   sva_nsv_method = "leek")
head(as.data.frame(adjusted_data), n = c(5,5))
```

##### Harman and Mean Centering
Harman and Mean Centering do not require special parameters, however, using Harman on large datasets may require significant time. 


```{r, results='hide'}
adjusted_data <- batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
                                   meta = BatchFLEX::example_meta,
                                   correction_method = c("Harman", "Mean Centering"),
                                   batch.1 = "batchflex_study",
                                   variable_of_interest = "Major_Lineage")
head(as.data.frame(adjusted_data[["Harman_adjusted"]]), n = c(5,5))
head(as.data.frame(adjusted_data[["Mean Centering_adjusted"]]), n = c(5,5))
```

#### Batch correct within Batch_FLEX
`batch_correct` can be used within the `Batch_FLEX` function and automatically pipe the output to the `batch_evaluate` functions. 

```{r, fig.width = 10, fig.height = 10}
test_example <- Batch_FLEX(Batch_FLEX_function = c("preprocess_matrix","batch_correct", "batch_evaluate"), mat = example_mat, meta = example_meta, correction_method = "Limma", evaluation_method = "pca", batch.1 = "batchflex_study")
head(as.data.frame(test_example$data_matrices$Limma), n = c(5,5))
test_example$batch_evaluation$Unadjusted$batch1$Plots$pca
test_example$batch_evaluation$Limma_adjusted$batch1$Plots$pca
```

### Visualize

#### Batch evalute used individually with specific evaluation methods
BatchFLEX utilizes many different methods of evaluation to determine whether batch correction is necessary and to assess whether the correction effectively removed the batch effect. If a variable_of_interest is known then BatchFLEX also ensures it is maintained after correction. `batch_evaluate` can be used as a stand alone function. By default, all evaluation methods are conducted, however, if desired, individual evaluation methods can be chosen. By default, when applicable, a black and white, a colored by batch, and a colored by the variable of interest plot will be generated for all plots. If desired, the color_by parameter can be adjusted to the desired plot or plots of interest.


##### Principal component analysis
By default, pca will include all annotation types, which includes annotation by cluster and annotation by the meta file. This can be adjusted by changing the annotation parameter. If a cluster number is not provided, then BatchFLEX will use silhouette analysis to generate an ideal cluster number for each matrix file. 
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "pca", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", annotation = "meta")
test_evaluate$Plots$pca
```

##### Cluster analysis
By default, the cluster analysis will include an elbow plot, a silhouette plot and a dunn index plot. If all three are not desired, then the cluster_analysis_method parameter can be changed. 
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "cluster_analysis", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", cluster_analysis_method = "silhouette")
test_evaluate$Plots$cluster_analysis
```
##### Cluster heterogeneity and evenness
BatchFLEX also includes multiple methods to assess the heterogeneity and evenness of the clustering to ensure appropriate batch correction occurs. Ideally, clustering representing the variable of interest should have relatively equal representation of each batch, which would suggest that the batch effect is not driving the clustering. By default, cluster_he will use the ward.d2 method and will use the MAD variance method.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "cluster_HE", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", plot_title = "Unadjusted")
test_evaluate$Plots$cluster_HE
test_evaluate$Matrices$cluster_HE[-1]
```

##### Heatmap
BatchFLEX will generate a heatmap with clustering based on the same cluster method chosen in the heterogeneity and evenness analysis. By default, this will use the top 2000 most variable genes by MAD and will be annotated by the batch.1 and variable of interest. However, the user can adjust this if desired. Additionally, the rows and columns can be labeled, however, for dense heatmaps, this is generally not recommended.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "heatmap", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")
test_evaluate$Plots$heatmap
```

##### PCA details
The PCA details analysis includes many helpful plots and matrices that can be used to assess the PCA. By default, `batch_evaluate` will use the batch.1 parameter when generating these files, however, this can be adjusted by altering the pca_factors parameter

```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "pca_details", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", pca_factors = "Major_Lineage")
test_evaluate$Plots$pca_details
```

##### Multiple components PCA
By default, the multiple components plot will include the top 5 principal components. This can be changes by altering the ncomponents parameter.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "mc_pca", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", color_by = "batch", ncomponents = 10)
test_evaluate$Plots$mc_pca
```

##### Relative log expression analysis
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "rle", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", color_by = "variable_of_interest")
test_evaluate$Plots$rle
```

##### Explanatory variables plot
By default, the explanatory variables plot will include the batch.1 and the variable of interest groups. If another column from the meta file is desired, then the variable_choices parameter can be altered.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "ev", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", variable_choices = c("batchflex_study", "Major_Lineage"))
test_evaluate$Plots$ev
```

##### Principal variance component analysis
The PVCA plot is generated from the top 20000 most variable genes using MAD by default. Additionally, pvca_pct is set to 0.8, but can be adjusted if desired. If variable_choices are not selected, the batch.1 and variable of interest will be used by default.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "pvca", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")
test_evaluate$Plots$ev
```

##### Surrogate variable analysis
By default, SVA will use the "be" method for surrogate variable estimation. If desired, the "leek" method can be chosen by chaning the sva_nsv_method.
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "sva", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage", sva_nsv_method = "be")
test_evaluate$Plots$sva
```

##### Unform manifold approximation and projection
```{r, fig.width = 10, fig.height = 10}
test_evaluate <- batch_evaluate(evaluation_method = "umap", mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")
test_evaluate$Plots$umap
```

#### Batch evaluate used individually with all evaluation methods
By default, `batch_evaluate` will include all evaluation methods that are available based on the input criteria. 

```{r, fig.width = 10, fig.height = 10, eval = FALSE}
test_evaluate <- batch_evaluate(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat), meta = example_meta, batch.1 = "batchflex_study", variable_of_interest = "Major_Lineage")
test_evaluate$Plots
```

#### Batch evaluate used within BatchFLEX
`batch_evaluate` can also be used in conjunction with other functions within the top level `Batch_FLEX` function.
Additionally, `batch_evaluate` can be ran in parallel using the parallelize parameter if desired, which can reduce computational time at the expense of overhead. Default is set to FALSE.

```{r, fig.width = 10, fig.height = 10}
test_Batch_FLEX <- Batch_FLEX(Batch_FLEX_function = c("preprocess_matrix", "batch_correct", "batch_evaluate"), correction_method = "Limma", evaluation_method = "umap", mat = example_mat, meta = example_meta, batch.1 =  "batchflex_study", variable_of_interest = "Major_Lineage", parallelize = TRUE)
test_Batch_FLEX[["batch_evaluation"]][["Unadjusted_Log2"]][["batch1"]][["Plots"]]
test_Batch_FLEX[["batch_evaluation"]][["Limma"]][["batch1"]][["Plots"]]
```

```{r}
devtools::session_info()
```









# Disclamer

Copyright 2024 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.








