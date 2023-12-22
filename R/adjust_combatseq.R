#' Adjust ComBat_Seq
#'
#' Wrapper for `sva::ComBat_seq`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param treatment Column name from the meta file of the column that will be used for treatment information
#' @param log2_transformed logical whether the data is alrady transformed
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
adjust_combatseq = function(mat, meta, batch.1, treatment, log2_transformed){
  if (is.null(batch.1)) stop("Mean Centering requires a batch.1 input")
  if (is.null(treatment)){
    model_matrix <- NULL
  } else {
    total_covariates <- paste0(treatment,collapse = "+")
    model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(meta))
  }
  combatseq_corrected <- sva::ComBat_seq(
    if(!log2_transformed) log2(mat) else {mat},
    batch = meta[[batch.1]],
    covar_mod = model_matrix
  )
  log2(combatseq_corrected + 1)
}
