#' merge_data
#'
#' @param merge_matrix_files A list of matrix files to be merged. If supplied as a named list, batchflex will add a named column to the meta file with the specified names. If not, batch flex will add a generic name per meta file.
#' @param merge_meta_files A list of meta files to be merged
#' @param keep_all_genes A TRUE/FALSE statement to indicate how unmatched genes will be handled. If TRUE, then NAs will be inserted for unmatched genes. If FALSE, unmatched genes will be deleted. Default is set to FALSE
#'
#' @return A list of merged matrix and meta files
#' @export
#'
#' @examples
#' set.seed(333)
merge_data <- function(merge_matrix_files = NULL,
                        merge_meta_files = NULL,
                        keep_all_genes = FALSE
){
  if (is.null(merge_matrix_files)){
    message("Must provide a list of matrices and meta files simulatenously to match matrix and meta information")
  }
  if (!is.null(names(merge_matrix_files))){
    user_names <- names(merge_matrix_files)
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
  merged_meta <- merged_meta[match(colnames(merged_matrix), merged_meta[,1]),]
  merge_data$merged_meta <- merged_meta
  return(merge_data)
}
