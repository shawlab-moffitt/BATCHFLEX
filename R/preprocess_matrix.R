#' Pre Process Matrix
#'
#' @param mat Numeric matrix or data frame. The features can be rownames or the first column followed by the sample names as the column names with numeric count data as the values.
#' @param raw.counts Logical. TRUE indicates the input data is raw counts. FALSE indicated the input data is not raw counts.
#' @param raw.norm.method Character string of what method to use for raw count normalization. Supported methods are "TMM" or "upperquartile". If raw.counts is FALSE this will be ignored.
#' @param log2 Logical. If TRUE, the input data with be logged with the method "log2+1".
#' @param quantnorm Logical. If TRUE, the input data will be quantile normalized.
#' @param remove.duplicates Logical. If TRUE, if duplicate row features are found, they will be summarized to the row with the maximum average feature value.
#'
#' @return A numeric matrix that has been processed through the users input choices.
#' @export
#'
#' @examples
#' set.seed(101)
preprocess_matrix = function(mat = NULL,
                             raw.counts = FALSE,
                             raw.norm.method = NULL,
                             log2 = TRUE,
                             quantnorm = TRUE,
                             remove.duplicates = TRUE) {
  #checks to make sure the data is in the right format
  if (!raw.counts) raw.norm.method = NULL
  if (raw.counts)
    if (!raw.norm.method %in% c("TMM","upperquartile")) stop("raw.norm.method not found")
  #if (is.data.frame(data))
    if (is.character(as.data.frame(mat)[,1]) & all(apply(mat[,2:ncol(mat)],2,is.numeric))) {
      mat <- as.data.frame(mat)
      if (TRUE %in% duplicated(mat[, 1])) {
        if (remove.duplicates) {
          cat("\tRemoving duplicate genes\n")
          data_dup <- mat %>% dplyr::group_by(Genes) %>% dplyr::filter(n() > 1) %>% as.data.frame()
          data_nodup <- mat %>% dplyr::group_by(Genes) %>% dplyr::filter(n() == 1) %>% as.data.frame()
          data_dup_summ <- data_dup %>%
            dplyr::group_by_at(1) %>%
            dplyr::summarise_all(max) %>%
            as.data.frame()
          mat <- rbind(data_nodup,data_dup_summ)
        } else stop("Duplicate row features found. Please remove duplicates")
      }
      rownames(mat) <- mat[,1]
      mat <- mat[,-1]
      mat <- as.matrix(mat)
    } else if (all(apply(mat,2,is.numeric))) {
      mat <- as.matrix(mat)
    } else {
      stop("Must provide count data in one of the following formats:\n\tA numeric matrix\n\tA numeric data frame\n\tA dataframe with the first column consisting of features and the following columns of numeric count data")
    }

  if (raw.counts) {
    cat("\tNormalizing Raw Counts\n")
    mat <- edgeR::DGEList(counts = mat)
    mat <- edgeR::calcNormFactors(mat, method = raw.norm.method)
    mat <- edgeR::cpm(mat)
  }
  if (log2) {
    cat("\tLogging Matrix\n")
    mat <- log2(mat + 1)
  }
  if (quantnorm) {
    cat("\tQuantile Normalizing\n")
    mat <- preprocessCore::normalize.quantiles(mat, keep.names = T)
  }
  return(mat)
}
