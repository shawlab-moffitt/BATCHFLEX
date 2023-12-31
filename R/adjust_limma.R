#' Adjust Limma
#'
#' Wrapper for `limma::removeBatchEffects`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param treatment Column name from the meta file of the column that will be used for treatment information
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
adjust_limma = function(mat,
                        meta,
                        treatment,
                        batch.1,
                        batch.2,
                        log2_transformed = TRUE,
                        ...){
  if (is.null(treatment)){
    model_matrix <- NULL
  } else {
    total_covariates <- paste0(treatment,collapse = "+")
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
