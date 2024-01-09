#' merge_data
#'
#' @param merge_matrix_files A list of matrix files to be merged
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
  if (is.null(merge_matrix_files) | is.null(merge_meta_files)){
    message("Must provide a list of matrices and meta files simulatenously to match matrix and meta information")
  }
  merge_data <- list()
  merged_matrix <- base::merge(merge_matrix_files[[1]], merge_matrix_files[[2]], all = keep_all_genes)
  if (length(merge_matrix_files) > 2){
    for (file in merge_matrix_files[3:length(merge_matrix_files)]){
      merged_matrix <- base::merge(merged_matrix, file)
    }
  }
  merge_data$merged_matrix <- merged_matrix

  merged_meta <- base::merge(merge_meta_files[[1]], merge_meta_files[[2]], all = TRUE)
  if (length(merge_matrix_files) > 2){
    for (file in merge_meta_files[3:length(merge_meta_files)]){
      merged_meta <- base::merge(merged_meta, file, all = TRUE)
    }
  }
  merged_meta <- merged_meta[match(colnames(merged_matrix)[-1], merged_meta[,1]),]
  merge_data$merged_meta <- merged_meta

  return(merge_data)
}
