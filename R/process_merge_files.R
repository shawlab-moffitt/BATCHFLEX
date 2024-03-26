#' Process Merge Files
#'
#' @param files_tobe_merged_path A path to the folder containing the matrix and meta files.
#' @param meta_files_included A TRUE/FALSE to determine if the function should attempt to group a matrix and meta file together.
#' @param all_transnormed <- A TRUE/FALSE to determine if all matrices have been transformed and normalized equally.
#'
#' @return A list of matrix and meta files in appropriate format to be merged using `merge_data`
#' @export
#'
#' @examples
#' set.seed(333)
#' test_process_merge <- Batch_FLEX(Batch_FLEX_function = c("process_merge_files", "merge_data"), files_tobe_merged_path = "path_to_folder_containing_matrix_and_meta_files")
#'
process_merge_files <- function(files_tobe_merged_path,
                                meta_files_included = TRUE,
                                all_transnormed = TRUE){
  merge_files <- base::list.files(files_tobe_merged_path, recursive = TRUE)

  if (meta_files_included == TRUE){
    if(!(length(merge_files) %% 2 == 0)){
      stop("Number of matrix and meta files do not match. Please check the containing folder")
    }
  }
  if (meta_files_included == TRUE){
    names_to_match <- tolower(basename(merge_files))
    names_to_match <- gsub("meta", "", names_to_match)
    names_to_match <- gsub("matrix", "", names_to_match)
    names_to_match <- gsub("[[:punct:]]", "", names_to_match)
    names_to_match <- gsub("txt", "", names_to_match)
    distance_measurements  <- adist(names_to_match)
    grouped_names <- data.frame(names_to_match, "grouping" = anticlust::matching(distance_measurements, p = 2))
    grouped_names$OrigID <- merge_files
    grouped_names <- grouped_names |> dplyr::relocate(OrigID)

    merge_files_data <- list()
    for (matched_files in 1:length(unique(grouped_names$grouping))){
      grouped_files <- grouped_names$OrigID[which(grouped_names$grouping == matched_files)]
      grouped_files_data <- sapply(grouped_files, function(end_path) read.delim(paste(files_tobe_merged_path, "/", end_path, sep = "")))
      names(grouped_files_data) <- basename(names(grouped_files_data))
      file_sizes <- data.frame("file_name" = names(grouped_files_data), "object_size" = sapply(grouped_files_data, function (object) object.size(object)))
      probable_matrix <- grouped_files_data[which(names(grouped_files_data) == file_sizes$file_name[which.max(file_sizes$object_size)])]
      probable_meta <- grouped_files_data[which(names(grouped_files_data) == file_sizes$file_name[which.min(file_sizes$object_size)])]
      if (ncol(probable_matrix[[1]][-1]) != nrow(probable_meta[[1]])){
        message(paste("The number of samples in",
                      names(probable_matrix),
                      "are",
                      ncol(probable_matrix[[1]][-1]),
                      "and the number of samples in",
                      names(probable_meta),
                      "are",
                      nrow(probable_meta[[1]]),
                      "and are not equal."))
        message("Attempting to find common direct matches between files.")
        probable_matrix_genes <- probable_matrix[[1]][[1]]
        probable_meta[[1]] <- probable_meta[[1]][which(probable_meta[[1]][[1]] %in% colnames(probable_matrix[[1]])),]
        probable_matrix[[1]] <- probable_matrix[[1]][,which(colnames(probable_matrix[[1]]) %in% probable_meta[[1]][[1]])]
        probable_matrix[[1]]$Genes <- probable_matrix_genes
        probable_matrix[[1]] <- probable_matrix[[1]] |> dplyr::relocate(Genes)
        message(paste("Common samples between files are ", ncol(probable_matrix[[1]][-1]), ".", sep = ""))
      }
      if (ncol(probable_matrix[[1]][-1]) == 0){
        message("File is being skipped as no direct matches were found.")
        message("To try flexable matching, please ensure sample numbers are equal between the matrix and meta file")
        next
      }
      merge_files_data$matrix <- append(merge_files_data$matrix, probable_matrix)
      merge_files_data$meta <- append(merge_files_data$meta, probable_meta)
    }
  }else{
    merge_files_data$matrix <- sapply(merge_files, function(end_path) read.delim(paste(files_tobe_merged_path, "/", end_path, sep = "")))
  }
  for (raw_matrix in 1:length(merge_files_data$matrix)){
    raw_matrix_name <- names(merge_files_data$matrix[raw_matrix])
    raw_matrix_data <- merge_files_data$matrix[[raw_matrix]]
    raw_matrix_final_rownames <- as.vector(raw_matrix_data[,1])
    raw_matrix_final <- as.matrix(raw_matrix_data[,-1])
    raw_matrix_final <- apply(raw_matrix_final, 2, as.numeric)
    if (all_transnormed == FALSE){
      mean_expression <- apply(raw_matrix_final, 2, function(x) mean(abs(x)))
      percent_above10 <- length(median_expression[which(median_expression > 10)])/length(median_expression)
      if (percent_above10 > 0.5){
        raw_matrix_final <- log2(raw_matrix_final + 1)
        raw_matrix_final <- as.matrix(raw_matrix_final)
      }
    }
    rownames(raw_matrix_final) <- raw_matrix_final_rownames
    merge_files_data$matrix[[raw_matrix_name]] <- raw_matrix_final
  }
  for (matrix in 1:length(merge_files_data$matrix)){
    matrix_name <- names(merge_files_data$matrix[matrix])
    new_data_frame <- as.data.frame(merge_files_data$matrix[[matrix_name]])
    new_data_frame$Genes <- row.names(merge_files_data$matrix[[matrix_name]])
    new_data_frame <- new_data_frame |> dplyr::relocate(Genes)
    merge_files_data$matrix[[matrix_name]] <- new_data_frame
  }
  merge_files_data$individual_matrices <- names(merge_files_data$matrix)
  merge_files_data$individual_meta_files <- names(merge_files_data$meta)
  return(merge_files_data)
}


