#' Adjust Harmon
#'
#' Remove batch effects from expression data using Harman from `Harman`.
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param variable_of_interest Column name from the meta file of the column that will be used for variable_of_interest information
#' @param log2_transformed logical whether the data is alrady transformed
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
#' adjusted_data <- batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
#' meta = BatchFLEX::example_meta,
#' correction_method = c("Harman"),
#' batch.1 = "batchflex_study",
#' variable_of_interest = "Major_Lineage")
#' head(as.data.frame(adjusted_data[["Harman_adjusted"]]), n = c(5,5))
#'
#'
adjust_harman = function(mat,
                         meta,
                         batch.1,
                         variable_of_interest,
                         log2_transformed){
  harman_correction_PCA <- Harman::harman(
    if(!log2_transformed) log2(mat) else {mat},
    expt = meta[[variable_of_interest]],
    batch = meta[[batch.1]],
    limit = 0.95,
    printInfo = T
  )

  Harman::reconstructData(harman_correction_PCA)
}
