#' Batch Correction
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param method A character vector of batch correction methods in c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg")
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param batch.2 Column name from the meta file of the column that will be used for batch two information
#' @param log2_transformed logical whether the data is already transformed
#' @param treatment Column name from the meta file of the column that will be used for treatment information
#' @param housekeeping Name of housekeeping gene set or character vector of housekeeping genes
#' @param k Used in the RUVg method, the number of factors of unwanted variation to be estimated from the data
#' @param drop Used in the RUVg method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k
#' @param center Used in the RUVg method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression
#' @param round Used in the RUVg method, if TRUE, the normalized measures are rounded to form pseudo-counts
#' @param tolerance Used in the RUVg method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive
#' @param par.prior Used in the ComBat method, TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#'
#' @return List object of length of method
#' @export
#'
#' @examples
#' set.seed(101)
batch_correction = function(mat = NULL,
                            meta = NULL,
                            method,
                            batch.1 = NULL,
                            batch.2 = NULL,
                            log2_transformed = TRUE,
                            treatment = NULL,
                            housekeeping = NULL,
                            k = 2,
                            drop = 0,
                            center = FALSE,
                            round = FALSE,
                            tolerance = 1e-8,
                            par.prior = TRUE) {
  #checks to make sure the data is in the right format
  if (is.null(mat)) stop("Must provide matrix")
  if (!all(apply(mat,2,is.numeric)) | !is(mat,"matrix")) stop("Must be numeric matrix")
  if (is.null(meta)) stop("Must provide meta data")
  if (!all(meta[,1] %in% colnames(mat))) stop("Sample names must match between matrix and meta data")
  if("all" %in% method) method = c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg")
  if (!all(method %in% c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg"))) stop("Batch correction methods not found")
  #if (!is.null(housekeeping)) stop("Must provide name of housekeeping gene set or vector of housekeeping genes")
  if(!is.null(batch.1))
    if(!(batch.1 %in% colnames(meta)))
        stop("batch.1 needs to be a column name in the metadata")
  if(length(batch.1) > 1) stop('batch.1 can only be of length 1 or NULL')
  if(!is.null(batch.2))
    if(!(batch.2 %in% colnames(meta)))
      stop("batch.2 needs to be a column name in the metadata")
  if(length(batch.2) > 1) stop('batch.2 can only be of length 1 or NULL')

  meta <- meta[match(colnames(mat), meta[[1]]),]
  batch_corrected_list <- list()
  batch_corrected_list$Unadjusted = mat
  if ("Limma" %in% method){
    cat("\tAdjusting Limma\n")
    batch_corrected_list$Limma <- adjust_limma(mat, meta, treatment, batch.1, batch.2)
  }
  if("ComBat" %in% method){
    cat("\tAdjusting ComBat\n")
    batch_corrected_list$ComBat = adjust_combat(mat, meta, batch.1, log2_transformed, par.prior)
  }
  if("Mean Centering" %in% method){
    cat("\tAdjusting Mean Centering\n")
    batch_corrected_list$`Mean Centering` = adjust_mean_centering(mat, meta, batch.1, log2_transformed)
  }
  if("ComBatseq" %in% method){
    cat("\tAdjusting ComBatseq\n")
    batch_corrected_list$ComBatseq <- adjust_combatseq(mat, meta, batch.1, treatment, log2_transformed)
  }
  if("Harman" %in% method){
    cat("\tAdjusting Harman\n")
    batch_corrected_list$Harman <- adjust_harman(mat, meta, batch.1, treatment, log2_transformed)
  }
  if("RUVg" %in% method){
    cat("\tAdjusting RUVg\n")
    batch_corrected_list$RUVg <- adjust_ruvg(mat, housekeeping, k, drop, center, round, tolerannce)
  }
  return(batch_corrected_list)
}
