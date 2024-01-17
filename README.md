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

```{r setup}
library(BatchFLEX)
library(ggpubr)
```

### BatchFLEX with Simulated data
#### Simulated data Alone
For quick evaluation of the BatchFLEX tool, a `simulate_data` function is included within BatchFLEX. This tool can generate simulated data with multiple treatment effects and with multiple batch effects based on user specified parameters. Then, a meta file is automatically generated based on the batch information. If `simulate_data` is used alone, a list of the matrx and meta files is exported, which can then be used as inputs in `Batch_FLEX`. Defaults for `simulate_data` can be found within the help page by using `?simulate_data()`. If `simulate_data` is called within the `Batch_FLEX` function, then the sim_mat and sim_meta files are automatically assigned to the mat and meta file for `Batch_FLEX`. Additionally, batch.1 and variable_of_interest is automatically assigned. In the example below, a single correction method and the rle evaluation method are being used for simplicity. 

The simulated data is based on the linear model framework introduced by Gagnon-Bartsch and Speed.

https://academic.oup.com/biostatistics/article/17/1/16/1744198


```{r}
simulated_data = simulate_data()
head(as.data.frame(simulated_data$sim_matrix))
head(simulated_data$sim_meta)
test_simulate_BatchFLEX <- Batch_FLEX(BatchFLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "rle", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", color_by = "batch")

print(test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["rle"]])
print(test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["rle"]])
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
                          batch_sample_sd = 0,
                          batch_sample_sd_max = 2,
                          epsilon_mean = 0,
                          epsilon_sd = 1)

test_simulate_BatchFLEX <- Batch_FLEX(BatchFLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "rle", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", variable_of_interest = "treatment")

print(test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["rle"]])
print(test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["rle"]])

```


#### Simulate data within BatchFLEX
```{r}
test_simulate_BatchFLEX <- Batch_FLEX(BatchFLEX_function = c("simulate_data", "batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "rle", color_by = "batch")

print(test_simulate_BatchFLEX[["batch_evaluation"]][["Unadjusted"]][["batch1"]][["Plots"]][["rle"]])
print(test_simulate_BatchFLEX[["batch_evaluation"]][["ComBat"]][["batch1"]][["Plots"]][["rle"]])
```


### BatchFLEX with merged data
`BathFLEX` comes preloaded with multiple datasets downloaded from GEO with meta files that have been manually adjusted to enable merging. These files can be merged using `merge_data` to generate a matrix and meta file similar to the example matrix and meta file included in `BatchFLEX`. 

```{r}
matrices_to_be_merged <- list(BatchFLEX::GSE112876_matrix,
                              BatchFLEX::GSE15907_matrix,
                              BatchFLEX::GSE37448_matrix,
                              BatchFLEX::GSE60336_matrix,
                              BatchFLEX::GSE60337_matrix,
                              BatchFLEX::GSE75195_matrix,
                              BatchFLEX::GSE75202_matrix,
                              BatchFLEX::GSE75203_matrix,
                              BatchFLEX::gse
```


#### Merge data alone



### BatchFLEX with example input data
#### Example Data
`BatchFLEX` has example data from several GEO studies that can be explored in the `?example_mat` help page. The number of samples was subset to include 48 to keep processing time and data size small. The data object `example_mat` includes the expression data for these 48 samples where the samples are in columns and genes in the rows. Sample level information can be found in `example_meta` such as the cell type the sample is and the study from which it came. 

```{r}
data(example_meta)
data(example_mat)
```

### Batch correction

RUVg takes in house keeping genes to for.... BatchFlex has a few house keeping genes to select from for both mouse and human.

Mouse:

Human: 
- `data(eisenberg)`
- `data(lin500)`
- `data(hsiao)`

The important parameters for adjusted batches are the expression matrix, the meta data associated with those samples, and which methods to use for adjustment. Each method has slightly different parameters that are required. 

```{r}
adjusted_data = batch_correction(mat = example_mat,
                                 meta = example_meta,
                                 method = c("Limma", "ComBat"),
                                 batch.1 = "Study",
                                 log2_transformed = TRUE)
```

### Visualize

The function `plot_umaps` takes in the list from `batch_correction` and then creates a plot list for all of the method that are specified. The `color` parameter can be used to color the points based on a column in the meta data that you are interested in. Here we will look at "Study" since that is what we adjusted for with in the batch correction function. This function is a nice wrapper for `ggplot2`, meaning that the output plots can then be further colored/added to the same way that one would with other `ggplot2` plots.

```{r}
plots = plot_umaps(adjusted_list = adjusted_data,
                   meta = example_meta,
                   color = "Study",
                   method = "Limma")
```

`ggpubr` has a functioned called `ggarrange` that allows for combining plots into plates, so can visualize both the unadjusted and the Limma adjusted expression together. This allows us to see the impact that 1) the study has on the expression profiles and 2) the 

```{r, fig.width = 8}
ggarrange(plots$Unadjusted, plots$Limma, common.legend = TRUE)
```

# Disclamer

Copyright 2024 Moffitt Cancer Center
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.








