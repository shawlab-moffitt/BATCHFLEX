#' Adjust RUVg
#'
#' Wrapper for `RUVSeq::RUVg`
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param housekeeping Name of housekeeping gene set or character vector of housekeeping genes
#' @param k Used in the RUVg method, the number of factors of unwanted variation to be estimated from the data
#' @param drop Used in the RUVg method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k
#' @param center Used in the RUVg method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression
#' @param round Used in the RUVg method, if TRUE, the normalized measures are rounded to form pseudo-counts
#' @param tolerance Used in the RUVg method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive
#'
#' @return an adjusted expression matrix
#' @export
#'
#' @examples
#' set.seed(333)
adjust_ruvg = function(mat,
                       housekeeping,
                       k,
                       drop,
                       center,
                       round,
                       tolerannce){
  h_i_m = housekeeping %in% row.names(mat)
  message(sum(h_i_m), " genes of ", length(h_i_m), " found in data")
  RUVg_correction <- RUVSeq::RUVg(
    if(!log2_transformed) log2(mat) else {mat},
    cIdx = housekeeping,
    k = k,
    drop = drop,
    center = center,
    round = round,
    tolerance = tolerance,
    isLog = T #forcing true due to transformation of the matrix
  )
  RUVg_correction$normalizedCounts
}
