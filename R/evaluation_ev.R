evaluation_ev <- function(mat,
                          rawmat,
                          batch_correction,
                          meta,
                          variable_choices,
                          batch,
                          variable_of_interest){
  evaluation_ev_list <- list()
  if (is.null(variable_choices)){
    variable_choices = c(batch, variable_of_interest)
  }
  my_colors <- metafolio::gg_color_hue(length(variable_choices))
  uncorrected_mat <- mat
  names(uncorrected_mat) <- NULL
  uncorrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(uncorrected_mat)),
    colData = meta,
    rowData = rawmat[,1]
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
    corrected_mat <- mat
    names(corrected_mat) <- NULL
    corrected_EV_SCE <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(corrected_mat)),
      colData = meta,
      rowData = rawmat[,1]
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
