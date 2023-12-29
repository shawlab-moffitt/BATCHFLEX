evaluation_mc_pca <- function(mat,
                              rawmat,
                              batch_correction,
                              meta,
                              batch,
                              variable_of_interest,
                              ncomponents,
                              color_by){
  mc_pca_plot_list <- list()
  color_mc_pca <- c()
  if (color_by == "batch"){
    color_mc_pca <-  batch
  }else if (color_by == "variable_of_interest"){
    color_mc_pca <-  variable_of_interest
  }
  unc_mat_pca_mc <- mat
  names(unc_mat_pca_mc) <- NULL
  uncorrected_pca_mc_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(unc_mat_pca_mc)),
    colData = meta,
    rowData = rawmat[,1]
  )
  SummarizedExperiment::assay(uncorrected_pca_mc_sce, "logcounts") <- SingleCellExperiment::counts(uncorrected_pca_mc_sce)
  uncorrected_pca_mc_sce <- scater::runPCA(uncorrected_pca_mc_sce, ncomponents = 50)
  mc_pca_plot_list$uncmcpca <- scater::plotPCA(
    uncorrected_pca_mc_sce,
    ncomponents = ncomponents,
    colour_by = color_mc_pca
  )
  if (!is.null(batch_correction)){
    cor_mat_pca_mc <- batch_correction
    names(cor_mat_pca_mc) <- NULL
    corrected_pca_mc_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as.matrix(cor_mat_pca_mc)),
      colData = meta,
      rowData = rawmat[,1]
    )
    SummarizedExperiment::assay(corrected_pca_mc_sce, "logcounts") <- SingleCellExperiment::counts(corrected_pca_mc_sce)
    corrected_pca_mc_sce <- scater::runPCA(corrected_pca_mc_sce, ncomponents = 50)
    mc_pca_plot_list$cormcpca <- scater::plotPCA(
      corrected_pca_mc_sce,
      ncomponents = ncomponents,
      colour_by = color_mc_pca)
  }
}
