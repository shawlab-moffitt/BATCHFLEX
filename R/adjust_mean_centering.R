#' Adjust Mean Centering
#'
#' Center each feature within each batch to zero mean using base R.
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
#' adjusted_data <- batch_correct(mat = BatchFLEX::preprocess_matrix(BatchFLEX::example_mat),
#' meta = BatchFLEX::example_meta,
#' correction_method = "Mean Centering",
#' batch.1 = "batchflex_study",
#' variable_of_interest = "Major_Lineage")
#' head(as.data.frame(adjusted_data[["Mean Centering_adjusted"]]), n = c(5,5))
#'
adjust_mean_centering = function(mat,
                                 meta,
                                 batch.1,
                                 log2_transformed){
  if (is.null(batch.1)) stop("Mean Centering requires a batch.1 input")
  m <- as.matrix(if (!log2_transformed) log2(mat) else mat)
  batch <- as.factor(meta[[batch.1]])

  if (length(batch) != ncol(m) || anyNA(batch)) {
    stop("Batch metadata must contain one nonmissing value per sample.")
  }

  if (any(!is.finite(m))) {
    stop("Mean Centering requires finite expression values after transformation.")
  }

  for (idx in split(seq_len(ncol(m)), batch, drop = TRUE)) {
    block <- m[, idx, drop = FALSE]
    m[, idx] <- sweep(block, 1L, rowMeans(block), FUN = "-")
  }

  m
}
