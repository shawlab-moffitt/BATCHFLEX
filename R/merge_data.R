#' merge_data
#'
#' Merges multiple expression matrices and generates a simple meta file to be used by `Batch_FLEX`, if none is provided.
#' If meta files are provided, they will be merged and a simple meta file will not be generated.
#' Samples from the merged meta file and matrix file will be matched using `RecordLinkage::levenshteinSim` generated scores and will be ordered similarly.
#'
#'
#' @param merge_matrix_files A list of matrix files (as data frames) to be merged. If supplied as a named list, batchflex will add a named column to the meta file with the specified names. If not, batch flex will add a generic name per meta file.
#' @param merge_meta_files A list of meta files to be merged
#' @param keep_all_genes A TRUE/FALSE statement to indicate how unmatched genes will be handled. If TRUE, then NAs will be inserted for unmatched genes. If FALSE, unmatched genes will be deleted. Default is set to FALSE
#'
#' @return A list of merged matrix and meta files
#' @export
#'
#' @examples
#' set.seed(333)
#' matrices_to_be_merged <- list("GSE112876" = BatchFLEX::GSE112876_matrix,
#' "GSE15907" = BatchFLEX::GSE15907_matrix,
#' "GSE37448" = BatchFLEX::GSE37448_matrix,
#' "GSE60336" = BatchFLEX::GSE60336_matrix,
#' "GSE75195" = BatchFLEX::GSE75195_matrix,
#' "GSE75202" = BatchFLEX::GSE75202_matrix,
#' "GSE75203" = BatchFLEX::GSE75203_matrix)
#' meta_files_to_be_merged <- list(BatchFLEX::GSE112876_meta,
#'                                BatchFLEX::GSE15907_meta,
#'                                BatchFLEX::GSE37448_meta,
#'                                BatchFLEX::GSE60336_meta,
#'                                BatchFLEX::GSE75195_meta,
#'                                BatchFLEX::GSE75202_meta,
#'                                BatchFLEX::GSE75203_meta)
#' test_merge <- merge_data(merge_matrix_files = matrices_to_be_merged, merge_meta_files = meta_files_to_be_merged)
#' head(as.data.frame(test_merge$merged_matrix), n = c(5,5))
#' head(test_merge$merged_meta, n = c(5,5))
#'
merge_data <- function(merge_matrix_files = NULL,
                        merge_meta_files = NULL,
                        keep_all_genes = FALSE){

  if (is.null(merge_matrix_files)){
    message("Must provide a list of matrices and meta files simulatenously to match matrix and meta information")
  }
  if (!is.null(names(merge_matrix_files))){
    user_names <- names(merge_matrix_files)
  }
  num_matrix_samples <- sum(unlist(lapply(merge_matrix_files, ncol))) - length(merge_matrix_files)
  num_meta_samples <- sum(unlist(lapply(merge_meta_files, nrow)))
  if (num_matrix_samples != num_meta_samples){
    stop("Number of matrix and meta samples are not equivalent. Please check the input files.")
  }
  merge_data <- list()
  merged_matrix <- base::merge(merge_matrix_files[[1]], merge_matrix_files[[2]], all = keep_all_genes)
  if (length(merge_matrix_files) > 2){
    for (file in merge_matrix_files[3:length(merge_matrix_files)]){
      merged_matrix <- base::merge(merged_matrix, file, all = keep_all_genes)
    }
  }
  row.names(merged_matrix) <- merged_matrix[,1]
  merged_matrix <- as.matrix(merged_matrix[,-1])
  merge_data$merged_matrix <- merged_matrix
  if (is.null(merge_meta_files)){
    study_vector <- vector()
    merged_meta <- data.frame()
    merged_meta <- rbind(merged_meta, as.data.frame(colnames(merged_matrix)))
    names(merged_meta) <- "SampleID"
    for (study_number in 1:length(merge_matrix_files)){
      batchflex_study <- paste0("Study_", study_number, sep = "")
      if (!is.null(names(merge_matrix_files))){
        batchflex_study <- user_names[study_number]
      }
      study_vector <- base::append(study_vector, rep(batchflex_study, length(colnames(merge_matrix_files[[study_number]])[-1])))
    }
    merged_meta$batchflex_study <- study_vector
  }else {
    for (study_number in 1:length(merge_meta_files)){
      for (meta_file in 1:length(merge_meta_files)){
        colnames(merge_meta_files[[meta_file]])[1] <- "OrigMetaID"
      }
      batchflex_study <- paste0("Study_", study_number, sep = "")
      if(!is.null(names(merge_matrix_files))){
        batchflex_study <- user_names[study_number]
      }
      study_vector <- vector()
      study_vector <- base::append(study_vector, rep(batchflex_study, nrow(merge_meta_files[[study_number]])))
      merge_meta_files[[study_number]]$batchflex_study <- study_vector
    }
    merged_meta <- base::merge(merge_meta_files[[1]], merge_meta_files[[2]], all = TRUE)
    if (length(merge_matrix_files) > 2){
      for (file in merge_meta_files[3:length(merge_meta_files)]){
        merged_meta <- base::merge(merged_meta, file, all = TRUE)
      }
    }
  }
  matrix_colnames <- colnames(merged_matrix)
  match_scores <- sapply(merged_meta$OrigMetaID, function (x) sapply(matrix_colnames, function(y) RecordLinkage::levenshteinSim(y, x)))
  merged_meta$BestMatchID <- rownames(match_scores)[apply(match_scores, 2, which.max)]
  merged_meta$BestMatchScore <- apply(match_scores, 2, max)
  merged_meta <- merged_meta |> dplyr::relocate(BestMatchID)
  merged_mata <- merged_meta |> dplyr::relocate(BestMatchScore, .before = OrigMetaID)
  merged_meta <- merged_meta[match(colnames(merged_matrix), merged_meta$BestMatchID),]
  merge_data$merged_meta <- merged_meta
  return(merge_data)
}
