#' Evaluation Relative Log Expression
#'
#' @param mat A Numeric matrix or list of matrices after pre-processing and/or batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by evaluation multiple components, pca, and rle to select which feature will be used to color the individuals. Choices are "batch", "variable_of_interest", "BnW", and "all". Default is set to "all"
#'
#' @return A list object of RLE plots
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_rle <- function(mat,
                           meta,
                           batch.1,
                           variable_of_interest,
                           color_by){
  evaluation_rle_list <- list()
  rle_mat <- mat
  names(rle_mat) <- NULL
  if ("batch" %in% color_by){
    RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(rle_mat)),
      colData = meta,
      rowData = rownames(mat)
    )
    evaluation_rle_list$batch_colored <- scater::plotRLE(
      RLE_SCE,
      exprs_values = "counts",
      color_by = batch.1
    )
  }
  if ("variable_of_interest" %in% color_by){
    RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(rle_mat)),
      colData = meta,
      rowData = rownames(mat)
    )
    evaluation_rle_list$voi_colored <- scater::plotRLE(
      RLE_SCE,
      exprs_values = "counts",
      color_by = variable_of_interest
    )
  }
  if ("BnW" %in% color_by){
    RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(rle_mat)),
      colData = meta,
      rowData = rownames(mat)
    )
    evaluation_rle_list$BnW_colored <- scater::plotRLE(
      RLE_SCE,
      exprs_values = "counts"
    )
  }

  return(evaluation_rle_list)
}
