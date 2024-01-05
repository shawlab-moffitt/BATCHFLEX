#' Evaluation Multiple Components Plot
#'
#' @param mat A Numeric matrix or list of matrices after pre-processing and/or batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param ncomponents Used by the multiple components plot to select the number of principal components that will be plotted. Default is set to 5
#' @param color_by Used by evaluation multiple components, pca, and rle to select which feature will be used to color the individuals. Choices are "batch", "variable_of_interest", "BnW", and "all". Default is set to "all"
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

  mat_pca_mc <- mat
  names(mat_pca_mc) <- NULL
  pca_mc_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(mat_pca_mc)),
    colData = meta,
    rowData = rownames(mat)
  )
  SummarizedExperiment::assay(pca_mc_sce, "logcounts") <- SingleCellExperiment::counts(pca_mc_sce)
  pca_mc_sce <- scater::runPCA(pca_mc_sce, ncomponents = 50)

  if ("batch" %in% color_by){
    mc_pca_plot_list$mcpca$batch_colored <- scater::plotPCA(
      pca_mc_sce,
      ncomponents = ncomponents,
      colour_by = batch.1
    )
  }
  if ("variable_of_interest" %in% color_by){
    mc_pca_plot_list$mcpca$voi_colored <- scater::plotPCA(
      pca_mc_sce,
      ncomponents = ncomponents,
      colour_by = variable_of_interest
    )
  }
  if ("BnW" %in% color_by){
    mc_pca_plot_list$mcpca$bnw_colored <- scater::plotPCA(
      pca_mc_sce,
      ncomponents = ncomponents
    )
  }


  return(mc_pca_plot_list)
}
