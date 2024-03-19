#' Adjust Limma
#'
#' Remove batch effects from expression data using `removeBatchEffect` from `limma`.
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_of_interest Column name from the meta file of the column that will be used for variable_of_interest information
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param batch.2 Column name from the meta file of the column that will be used for batch two information
#' @param log2_transformed whether data is log2 transformed already or not
#' @param ... other parameters to pass to `limma::removeBatchEffects`
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
#' adjusted_data <-  batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
#' meta = BatchFLEX::example_meta,
#' correction_method = "Limma",
#' batch.1 = "batchflex_study",
#' batch.2 = NULL,
#' variable_of_interest = "Major_Lineage")
#' head(as.data.frame(adjusted_data), n = c(5,5))
#'
adjust_limma = function(mat,
                        meta,
                        variable_of_interest,
                        batch.1,
                        batch.2,
                        log2_transformed = TRUE,
                        ...){
  if (is.null(variable_of_interest)){
    model_matrix <- NULL
  } else {
    total_covariates <- paste0(variable_of_interest,collapse = "+")
    model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(meta))
  }
  #get the values associated with batch 1 and batch 2
  #input for 'removeBatchEffect is vectors
  b1 = if(is.null(batch.1)) NULL else {meta[[batch.1]]}
  b2 = if(is.null(batch.2)) NULL else {meta[[batch.2]]}
  limma::removeBatchEffect(
    if(!log2_transformed) log2(mat) else {mat},
    batch = b1,
    batch2 = b2,
    covariates = model_matrix,
    ...
  )
}
