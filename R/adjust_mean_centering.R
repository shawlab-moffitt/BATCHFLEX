#' Adjust Mean Centering
#'
#' Wrapper for `bapred::meancenter`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param log2_transformed whether data is log2 transformed already or not
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
adjust_mean_centering = function(mat,
                                 meta,
                                 batch.1,
                                 log2_transformed){
  if (is.null(batch.1)) stop("Mean Centering requires a batch.1 input")
  mean_centering_batch = meta[[batch.1]]
  mean_centering_data = t(if(!log2_transformed) log2(mat) else {mat})
  mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
  t(mean_center$xadj)
}
