#' Evaluation Relative Log Expression
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by the RLE plot to select whether "batch" or "variable_of_interest" will be used to color the individuals. Default is set to "batch"
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
  color_rle <- c()
  if (color_by == "batch"){
    color_rle <-  batch.1
  }else if (color_by == "variable_of_interest"){
    color_rle <-  variable_of_interest
  }
  rle_mat <- mat
  names(rle_mat) <- NULL
  RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(rle_mat)),
    colData = meta,
    rowData = rownames(mat)
  )
  evaluation_rle_list$rle <- scater::plotRLE(
    RLE_SCE,
    exprs_values = "counts",
    color_by = color_rle
  )
  return(evaluation_rle_list)
}
