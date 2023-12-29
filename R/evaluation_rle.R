evaluation_rle <- function(mat,
                           rawmat,
                           batch_correction,
                           meta,
                           batch,
                           variable_of_interest,
                           color_by){
  evaluation_rle_list <- list()
  color_rle <- c()
  if (color_by == "batch"){
    color_rle <-  batch
  }else if (color_by == "variable_of_interest"){
    color_rle <-  variable_of_interest
  }
  uncorrected_mat <- mat
  names(uncorrected_mat) <- NULL
  uncorrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(uncorrected_mat)),
    colData = meta,
    rowData = rawmat[,1]
  )
  evaluation_rle_list$uncrle <- scater::plotRLE(
    uncorrected_RLE_SCE,
    exprs_values = "counts",
    color_by = color_rle
  )
  if(!is.null(batch_correction)){
    corrected_mat <- batch_correction
    names(corrected_mat) <- NULL
    corrected_RLE_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(corrected_mat)),
      colData = meta,
      rowData = rawmat[,1]
    )
    evaluation_rle_list$corrle <- scater::plotRLE(
      corrected_RLE_SCE,
      exprs_values = "counts",
      color_by = color_rle
    )
  }
  return(evaluation_rle_list)
}
