#' Adjust_sva
#'
#' Wrapper for `sva::sva`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_of_interest Column name from the meta file of the column that will be used as the variable of interest to generate the model matrix
#' @param sva_nsv_method Input method for the num.sv function in sva. Default is set to "be", but can be manually set to "leek"
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
adjust_sva = function(mat,
                      meta,
                      variable_of_interest,
                      sva_nsv_method){
  if (is.null(variable_of_interest)) {
    stop("SVA requires a variable_of_interest input")}
  if (sva_nsv_method == "be") {
    print("Method 'be' is selected. Input 'leek' to sva_nsv_method if alternative is desired")
  }
  sva_mod <- stats::model.matrix(reformulate(variable_of_interest), data = meta)
  sva_mod0 <- stats::model.matrix(~1, data = meta)
  n_sv <- sva::num.sv(mat, sva_mod, method = sva_nsv_method)
  svobj <- sva::sva(mat, sva_mod, sva_mod0, n.sv = n_sv)
  fsvaobj <- sva::fsva(mat, sva_mod, svobj, mat)
  fsvaobj$db
}
