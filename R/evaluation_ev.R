#' Evaluation Explanatory Variables
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_choices Used by the explanatory variables function to select which variables to plot. Default is a combination of the batch and variable of interest
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#'
#' @return A list object of explanatory variables plots
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_ev <- function(mat,
                          meta,
                          variable_choices,
                          batch.1,
                          variable_of_interest){
  evaluation_ev_list <- list()
  if (is.null(variable_choices)){
    variable_choices = c(batch.1, variable_of_interest)
  }
  mat <- mat
  names(mat) <- NULL
  EV_SCE <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(mat)),
    colData = meta,
    rowData = rownames(mat)
  )
  SummarizedExperiment::assay(EV_SCE, "logcounts") <- SingleCellExperiment::counts(EV_SCE)
  EV_SCE_PCA <- scater::runPCA(EV_SCE)
  EV_SCE_PCA_Vars <- getVarianceExplained(EV_SCE_PCA,exprs_values = "logcounts",variables = variable_choices)
  EV_SCE_PCA_Vars_melt <- reshape2::melt(EV_SCE_PCA_Vars)
  evaluation_ev_list$ev <- ggplot(EV_SCE_PCA_Vars_melt, aes(x = value,color = Var2)) +
    geom_density(linewidth = 1) +
    scale_x_log10(limit = c(0.0001,100),labels = ~ format(.x, scientific = FALSE), breaks = c(0.001,0.01,0.1,1,10,100)) +
      geom_vline(xintercept = 1, linetype="dashed") +
    theme_classic() +
    labs(x = "% Variance Explained", y = "Density", color = "Variable")

  return(evaluation_ev_list)
}


