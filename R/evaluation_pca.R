#' Evaluation Pincipal Component Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param annotation Used to select whether a cluster or meta annotated PCA plot is generated. cluster = "cluster", meta = "meta", all = c("cluster", "meta")
#' @param cluster_number Used to select the number of kmeans generated clusters to display in the uncorrected plot. If NULL is selected, a Silhouette generated cluster number is used.
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by the PCA plot to select which feature will be used to color the individuals
#'
#' @return A list object of PCA plots for each annotation method selected
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_pca = function(mat,
                          meta,
                          annotation,
                          cluster_number,
                          batch.1,
                          variable_of_interest,
                          color_by){
  pca_plot_list <- list()
  row.names(mat) <- NULL
  if("all" %in% annotation) annotation = c("cluster", "meta")
  if (!all(annotation %in% c("cluster", "meta"))){
    stop("Annotation method not found")
  }
  if ("cluster" %in% annotation & is.null(cluster_number) ){
    silhouette <- factoextra::fviz_nbclust(x = t(mat), FUNcluster = kmeans, method = "silhouette", verbose = T)
    cluster_number = silhouette$data$clusters[which(max(silhouette$data$y) == silhouette$data$y)]
  }
  if ("meta" %in% annotation & is.null(batch.1)){
    print("Missing batch information")
  }
  if ("meta" %in% annotation & is.null(variable_of_interest)){
    print("Missing variable of interest")
  }
  if ("cluster" %in% annotation){
    PCA <- cluster::pam(as.data.frame(t(mat)), cluster_number)
    pca_plot_list$pca_clust <- ggplot2::autoplot(
      PCA,
      frame = TRUE,
      frame.type = "norm")
  }
  if ("meta" %in% annotation){
    PCA <- stats::prcomp(t(mat), scale. = TRUE)
    PCA_data <- cbind(t(mat), meta)
    if (color_by == "batch"){
      pca_plot_list$pcameta <- ggplot2::autoplot(
        PCA,
        PCA_data,
        color = batch.1,
        shape = variable_of_interest
      )+
        ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
    }else if (color_by == "voi"){
      pca_plot_list$pcameta <- ggplot2::autoplot(
        PCA,
        PCA_data,
        color = variable_of_interest,
        shape = batch.1
      )+
        ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
    }
  }
  return(pca_plot_list)
}

