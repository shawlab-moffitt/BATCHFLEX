#' Evaluation Explanatory Variables
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_choices Used by the explanatory variables function to select which variables to plot. Default is a combination of the batch and variable of interest.
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#'
#' @return A list object of uncorrected and batch corrected explanatory variables plots
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_ev <- function(mat,
                          batch_correction,
                          meta,
                          variable_choices,
                          batch.1,
                          variable_of_interest){
  evaluation_ev_list <- list()
  if (is.null(variable_choices)){
    variable_choices = c(batch.1, variable_of_interest)
  }
  my_colors <- metafolio::gg_color_hue(length(variable_choices))
  uncorrected_mat <- mat
  names(uncorrected_mat) <- NULL
  uncorrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(uncorrected_mat)),
    colData = meta,
    rowData = rownames(mat)
  )
  SummarizedExperiment::assay(uncorrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(uncorrected_EV_SCE)
  uncorrected_EV_SCE_PCA <- scater::runPCA(uncorrected_EV_SCE)
  evaluation_ev_list$uncev <- scater::plotExplanatoryVariables(
    uncorrected_EV_SCE_PCA,
    exprs_values = "logcounts",
    variables = variable_choices
  ) +
    ggplot2::scale_color_manual(values = my_colors)
  if(!is.null(batch_correction)){
    corrected_mat <- batch_correction
    names(corrected_mat) <- NULL
    corrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(corrected_mat)),
      colData = meta,
      rowData = rownames(mat)
    )
    SummarizedExperiment::assay(corrected_EV_SCE, "logcounts") <- SingleCellExperiment::counts(corrected_EV_SCE)
    corrected_EV_SCE_PCA <- scater::runPCA(corrected_EV_SCE)
    evaluation_ev_list$corev <- scater::plotExplanatoryVariables(
      corrected_EV_SCE_PCA,
      exprs_values = "logcounts",
      variables = variable_choices
    ) +
      ggplot2::scale_color_manual(values = my_colors)
  }
  return(evaluation_ev_list)
}
