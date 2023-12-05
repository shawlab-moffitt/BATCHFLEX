#' Batch Correction
#'
#' @param mat Numeric matrix with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information
#' @param batch.2 Column name from the meta file of the column that will be used for batch two information
#' @param treatment Column name from the meta file of the column that will be used for treatment information
#' @param housekeeping Name of housekeeping gene set or character vector of housekeeping genes
#' @param k Used in the RUVg method, the number of factors of unwanted variation to be estimated from the data
#' @param drop Used in the RUVg method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k
#' @param center Used in the RUVg method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression
#' @param round Used in the RUVg method, if TRUE, the normalized measures are rounded to form pseudo-counts
#' @param tolerance Used in the RUVg method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive
#' @param par.prior Used in the ComBat method, TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param method A character vector of batch correction methods
#'
#' @return List object of length of method
#' @export
#'
#' @examples
#' set.seed(101)
batch_correction = function(mat = NULL,meta = NULL,method,batch.1 = NULL,batch.2 = NULL,treatment = NULL,housekeeping = NULL,
                            k = 2,drop = 0,center = FALSE,round = FALSE,tolerance = 1e-8,par.prior = TRUE) {

  if (is.null(mat)) stop("Must provide matrix")
  if (!all(apply(mat,2,is.numeric)) | !is(mat,"matrix")) stop("Must be numeric matrix")
  if (is.null(meta)) stop("Must provide meta data")
  if (!all(method %in% c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg"))) stop("Batch correction methods not found")
  #if (!is.null(housekeeping)) stop("Must provide name of housekeeping gene set or vector of housekeeping genes")


  meta <- meta[match(colnames(mat), meta[[1]]),]
  batch_corrected_list <- list()
  if ("Limma" %in% method){
    if (is.null(treatment)){
      model_matrix <- NULL
    } else {
      total_covariates <- paste0(treatment,collapse = "+")
      model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(meta))
    }
    batch_corrected_list$limma_corrected <- limma::removeBatchEffect(
      mat,
      batch = if(is.null(batch.1)) NULL else {meta[,batch.1]},
      batch2 = if(is.null(batch.2)) NULL else {meta[,batch.2]},
      covariates = model_matrix
    )

  }
  if("ComBat" %in% method){
    if (!is.null(batch.1)) {
      batch_combat <- meta[,batch.1]
      modcombat <-  stats::model.matrix(~1, data = as.data.frame(meta))
      batch_corrected_list$ComBat_corrected <- sva::ComBat(
        dat = mat,
        batch = batch_combat,
        mod = modcombat,
        par.prior = par.prior
      )
    }
  }
  if("Mean Centering" %in% method){
    if (!is.null(batch.1)) {
      mean_centering_batch = meta[,batch.1]
      mean_centering_data = t(mat)
      mean_center <- bapred::meancenter(as.matrix(mean_centering_data), as.factor(mean_centering_batch))
      batch_corrected_list$mean_center_correction <- as.data.frame(t(mean_center$xadj))
    }
  }
  if("ComBatseq" %in% method){
    if (is.null(treatment)){
      model_matrix <- NULL
    } else {
      total_covariates <- paste0(treatment,collapse = "+")
      model_matrix <- stats::model.matrix(reformulate(total_covariates), data = as.data.frame(meta))
    }
    combatseq_corrected <- sva::ComBat_seq(
      (2^mat),
      batch = meta[,batch.1],
      covar_mod = model_matrix
    )
    batch_corrected_list$combatseq_corrected <- log2(combatseq_corrected + 1)
  }
  if("Harman" %in% method){
    harman_correction_PCA <- Harman::harman(
      mat,
      expt = meta[,treatment],
      batch = meta[,batch.1],
      limit = 0.95,
      printInfo = T
    )
    batch_corrected_list$harman_corrected <- Harman::reconstructData(harman_correction_PCA)
  }
  if("RUVg" %in% method){
    RUVg_correction <- RUVSeq::RUVg(
      mat,
      cIdx = housekeeping,
      k = k,
      drop = drop,
      center = center,
      round = round,
      tolerance = tolerance,
      isLog = T
    )
    batch_corrected_list$RUVg_corrected <- as.data.frame(RUVg_correction$normalizedCounts)
  }
  return(batch_corrected_list)
}
