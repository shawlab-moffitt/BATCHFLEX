#' Adjust ComBat
#'
#' Wrapper for `sva::ComBat`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param log2_transformed whether data is log2 transformed already or not
#' @param par.prior Used in the ComBat method, TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
adjust_combat = function(mat,
                         meta,
                         batch.1,
                         log2_transformed,
                         par.prior){
  if (is.null(batch.1)) stop("ComBat requires a batch.1 input")
  batch_combat <- meta[[batch.1]]
  modcombat <-  stats::model.matrix(~1, data = as.data.frame(meta))
  sva::ComBat(
    dat = if(!log2_transformed) log2(mat) else {mat},
    batch = batch_combat,
    mod = modcombat,
    par.prior = par.prior
  )
}
