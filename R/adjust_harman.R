#' Adjust Harmon
#'
#' Wrapper for `Harman::harman` and `Harman::reconstructData`
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
adjust_harman = function(mat,
                         meta,
                         batch.1,
                         treatment,
                         log2_transformed){
  harman_correction_PCA <- Harman::harman(
    if(!log2_transformed) log2(mat) else {mat},
    expt = meta[[treatment]],
    batch = meta[[batch.1]],
    limit = 0.95,
    printInfo = T
  )

  Harman::reconstructData(harman_correction_PCA)
}
