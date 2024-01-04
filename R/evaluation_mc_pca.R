#' Evaluation Multiple Components Plot
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param ncomponents Used by the multiple components plot to select the number of principal components that will be plotted. Default is set to 5
#' @param color_by Used by the multiple components plot to select whether "batch" or "variable_of_interest" will be used to color the individuals. Default is set to "batch"
#'
#' @return A list object of multiple components plots
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_mc_pca <- function(mat,
                              meta,
                              batch.1,
                              variable_of_interest,
                              ncomponents,
                              color_by){
  mc_pca_plot_list <- list()
  color_mc_pca <- c()
  if (color_by == "batch"){
    color_mc_pca <-  batch.1
  }else if (color_by == "variable_of_interest"){
    color_mc_pca <-  variable_of_interest
  }
  mat_pca_mc <- mat
  names(mat_pca_mc) <- NULL
  pca_mc_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(mat_pca_mc)),
    colData = meta,
    rowData = rownames(mat)
  )
  SummarizedExperiment::assay(pca_mc_sce, "logcounts") <- SingleCellExperiment::counts(pca_mc_sce)
  pca_mc_sce <- scater::runPCA(pca_mc_sce, ncomponents = 50)
  mc_pca_plot_list$mcpca <- scater::plotPCA(
    pca_mc_sce,
    ncomponents = ncomponents,
    colour_by = color_mc_pca
  )
  return(mc_pca_plot_list)
}
