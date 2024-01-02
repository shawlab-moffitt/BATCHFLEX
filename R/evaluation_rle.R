#' Evaluation Relative Log Expression
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param rawmat Numeric matrix before pre-processing with features as rownames and sample names as the column names
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by the RLE plot to select whether "batch" or "variable_of_interest" will be used to color the individuals. Default is set to "batch"
#'
#' @return A list object of uncorrected and batch corrected RLE plots
#' @export
#'
#' @examples
evaluation_rle <- function(mat,
                           rawmat,
                           batch_correction,
                           meta,
                           batch,
                           variable_of_interest,
                           color_by){
  evaluation_rle_list <- list()
  color_rle <- c()
  if (color_by == "batch"){
    color_rle <-  batch
  }else if (color_by == "variable_of_interest"){
    color_rle <-  variable_of_interest
  }
  uncorrected_mat <- mat
  names(uncorrected_mat) <- NULL
  uncorrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(uncorrected_mat)),
    colData = meta,
    rowData = rawmat[,1]
  )
  evaluation_rle_list$uncrle <- scater::plotRLE(
    uncorrected_RLE_SCE,
    exprs_values = "counts",
    color_by = color_rle
  )
  if(!is.null(batch_correction)){
    corrected_mat <- batch_correction
    names(corrected_mat) <- NULL
    corrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(corrected_mat)),
      colData = meta,
      rowData = rawmat[,1]
    )
    evaluation_rle_list$corrle <- scater::plotRLE(
      corrected_RLE_SCE,
      exprs_values = "counts",
      color_by = color_rle
    )
  }
  return(evaluation_rle_list)
}
