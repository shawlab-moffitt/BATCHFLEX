#' Batch Correct
#'
#' @param mat A Numeric matrix after preprocessing with features as rownames and sample names as the column names.
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix.
#' @param batch.1 Column name from the meta file of the column that will be used for batch one information.
#' @param batch.2 Column name from the meta file of the column that will be used for batch two information.
#' @param log2_transformed logical whether the data is already transformed. Default is set to TRUE.
#' @param variable_of_interest Column name from the meta file of the column that will be used for variable of interest information.
#' @param housekeeping Name of housekeeping gene set or character vector of housekeeping genes.
#' @param k Used in the RUVg correction_method, the number of factors of unwanted variation to be estimated from the data.
#' @param drop Used in the RUVg correction_method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k.
#' @param center Used in the RUVg correction_method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression.
#' @param round Used in the RUVg correction_method, if TRUE, the normalized measures are rounded to form pseudo-counts.
#' @param tolerance Used in the RUVg correction_method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive.
#' @param par.prior Used in the ComBat correction_method, TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used.
#' @param correction_method A character vector of batch correction methods in c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg", "SVA).
#' @param sva_nsv_method Input correction_method for the num.sv function in sva. Default is set to "be", but can be manually set to "leek".
#'
#' @return List object of length of correction_method
#' @export
#'
#' @examples
#' set.seed(333)
batch_correct = function(mat = NULL,
                         meta = NULL,
                         correction_method = "all",
                         batch.1 = NULL,
                         batch.2 = NULL,
                         log2_transformed = TRUE,
                         variable_of_interest = NULL,
                         housekeeping = BatchFLEX::hsiao_mouse,
                         k = 2,
                         drop = 0,
                         center = FALSE,
                         round = FALSE,
                         tolerance = 1e-8,
                         par.prior = TRUE,
                         sva_nsv_method = "be") {

  #checks to make sure the data is in the right format
  if (is.null(mat)) stop("Must provide matrix")
  if (!all(apply(mat,2,is.numeric)) | !is(mat,"matrix"))
    if (!prep_matrix) stop("Must be numeric matrix")
  if (is.null(meta)) stop("Must provide meta data")

  if (is.null(batch.1) & !"SVA" %in% correction_method) {
    stop("Please select column name in the meta file for the batch information")
  }
  if (is.null(variable_of_interest)){
    #message("Missing variable of interest. Some correction methods and evaluation techniques are not available.")
    message("Missing variable of interest. SVA and Harman correction methods are not available.")

  }
  if (is.null(housekeeping) & "RUVg" %in% correction_method | is.null(housekeeping) & "all" %in% correction_method){
    stop("Please provide a list of housekeeping genes for the RUVg correction method")
  }
  if("all" %in% correction_method) correction_method = c("Limma", "ComBat", "Mean Centering", "Harman", "RUVg", "SVA")
  if(is.null(variable_of_interest)) {
    remove = c("SVA", "Harman")
    correction_method = correction_method[!correction_method %in% remove]
  }
  if (!all(correction_method %in% c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg", "SVA"))) {
    stop("Batch correction method not found \nAcceptable methods include Limma, ComBat, Mean Centering, Combatseq, Harman, RUVg, and SVA")}
  if(!is.null(batch.1))
    if(!(batch.1 %in% colnames(meta)))
        stop("batch.1 needs to be a column name in the metadata")
  if(length(batch.1) > 1) stop('batch.1 can only be of length 1 or NULL')
  if(!is.null(batch.2))
    if(!(batch.2 %in% colnames(meta)))
      stop("batch.2 needs to be a column name in the metadata")
  if(length(batch.2) > 1) stop('batch.2 can only be of length 1 or NULL')

  batch_corrected_list <- list()

  meta <- meta[match(colnames(mat), meta[[1]]),]

  if ("Limma" %in% correction_method){
    cat("\tAdjusting Limma\n")
    batch_corrected_list$Limma_adjusted <- adjust_limma(mat, meta, variable_of_interest, batch.1, batch.2)
  }
  if("ComBat" %in% correction_method){
    cat("\tAdjusting ComBat\n")
    batch_corrected_list$ComBat_adjusted = adjust_combat(mat, meta, batch.1, log2_transformed, par.prior)
  }
  if("Mean Centering" %in% correction_method){
    cat("\tAdjusting Mean Centering\n")
    batch_corrected_list$`Mean Centering_adjusted` = adjust_mean_centering(mat, meta, batch.1, log2_transformed)
  }
  if("ComBatseq" %in% correction_method){
    cat("\tAdjusting ComBatseq\n")
    batch_corrected_list$ComBatseq_adjusted <- adjust_combatseq(mat, meta, batch.1, variable_of_interest, log2_transformed)
  }
  if("Harman" %in% correction_method){
    cat("\tAdjusting Harman\n")
    batch_corrected_list$Harman_adjusted <- adjust_harman(mat, meta, batch.1, variable_of_interest, log2_transformed)
  }
  if("RUVg" %in% correction_method){
    cat("\tAdjusting RUVg\n")
    batch_corrected_list$RUVg_adjusted <- adjust_ruvg(mat, housekeeping, k, drop, center, round, tolerance, log2_transformed)
  }
  if("SVA" %in% correction_method){
    cat("\tAdjusting SVA\n")
    batch_corrected_list$SVA_adjusted <- adjust_sva(mat, meta, variable_of_interest, sva_nsv_method)
  }
  return(batch_corrected_list)
}
