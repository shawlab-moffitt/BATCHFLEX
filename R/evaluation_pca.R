#' Evaluation Pincipal Component Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param annotation Used to select whether a cluster or meta annotated PCA plot is generated. cluster = "cluster", meta = "meta", all = c("cluster", "meta")
#' @param uncorrected_cluster_number Used to select the number of kmeans generated clusters to display in the uncorrected plot. If NULL is selected, a Dunn generated cluster number is used.
#' @param corrected_cluster_number Used to select the number of kmeans generated clusters to display in the corrected plot. If NULL is selected, a Dunn generated cluster number is used.
#' @param batch Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by the PCA plot to select which feature will be used to color the individuals
#'
#' @return A list object of uncorrected and batch corrected PCA plots for each annotation method selected
#' @export
#'
#' @examples
evaluation_pca = function(mat,
                          batch_correction,
                          meta,
                          annotation,
                          uncorrected_cluster_number,
                          corrected_cluster_number,
                          batch,
                          variable_of_interest,
                          color_by){
  pca_plot_list <- list()
  if("all" %in% annotation) annotation = c("cluster", "meta")
  if (!all(annotation %in% c("cluster", "meta"))){
    stop("Annotation method not found")
  }
  if ("cluster" %in% annotation & is.null(uncorrected_cluster_number) ){
    uncorrected_dunn_k <- c(2:10)
    uncorrected_dunnin <- c()
    for (i in uncorrected_dunn_k){
      uncorrected_dunnin[i] <- clValid::dunn(
        distance = dist(t(mat)),
        clusters = kmeans(t(mat), i)$cluster
      )
    }
    uncorrected_dunn_index_analysis <- as.data.frame(cbind(uncorrected_dunn_k, uncorrected_dunnin[-1]))
    colnames(uncorrected_dunn_index_analysis) <- c("cluster_number", "dunn_index")
    uncorrected_cluster_number = uncorrected_dunn_index_analysis$cluster_number[which(max(uncorrected_dunn_index_analysis$dunn_index) == uncorrected_dunn_index_analysis$dunn_index)]
  }
  if ("cluster" %in% annotation & !is.null(batch_correction) & is.null(corrected_cluster_number) ){
    corrected_dunn_k <- c(2:10)
    corrected_dunnin <- c()
    for (i in corrected_dunn_k){
      corrected_dunnin[i] <- clValid::dunn(
        distance = dist(t(batch_correction)),
        clusters = kmeans(t(batch_correction), i)$cluster
      )
    }
    corrected_dunn_index_analysis <- as.data.frame(cbind(corrected_dunn_k, corrected_dunnin[-1]))
    colnames(corrected_dunn_index_analysis) <- c("cluster_number", "dunn_index")
    corrected_cluster_number = corrected_dunn_index_analysis$cluster_number[which(max(corrected_dunn_index_analysis$dunn_index) == corrected_dunn_index_analysis$dunn_index)]
  }
  if ("meta" %in% annotation & is.null(batch)){
    print("Missing batch information")
  }
  if ("meta" %in% annotation & is.null(variable_of_interest)){
    print("Missing variable of interest")
  }
  if ("cluster" %in% annotation){
    uncorrected_PCA <- cluster::pam(as.data.frame(t(mat)), uncorrected_cluster_number)
    pca_plot_list$uncpcaclust <- ggplot2::autoplot(
      uncorrected_PCA,
      frame = TRUE,
      frame.type = "norm")
    if (!is.null(batch_correction)){
      corrected_PCA <- cluster::pam(as.data.frame(t(batch_correction)), corrected_cluster_number)
      pca_plot_list$corpcaclust <- ggplot2::autoplot(
        corrected_PCA,
        frame = TRUE,
        frame.type = "norm")
    }
  }
  if ("meta" %in% annotation){
    uncorrected_PCA <- stats::prcomp(t(mat), scale. = TRUE)
    uncorrected_PCA_data <- cbind(t(mat), meta)
    if (color_by == "batch"){
      pca_plot_list$uncpcameta <- ggplot2::autoplot(
        uncorrected_PCA,
        uncorrected_PCA_data,
        color = batch,
        shape = variable_of_interest
      )+
        ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
    }else if (color_by == "voi"){
      pca_plot_list$uncpcameta <- ggplot2::autoplot(
        uncorrected_PCA,
        uncorrected_PCA_data,
        color = variable_of_interest,
        shape = batch
      )+
        ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
    }
    if (!is.null(batch_correction)){
      corrected_PCA <- stats::prcomp(t(batch_correction), scale. = TRUE)
      corrected_PCA_data <- cbind(t(batch_correction), meta)
      if (color_by == "batch"){
        pca_plot_list$corpcameta <- ggplot2::autoplot(
          corrected_PCA,
          corrected_PCA_data,
          color = batch,
          shape = variable_of_interest
        )+
          ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
      }else if (color_by == "voi"){
        pca_plot_list$corpcameta <- ggplot2::autoplot(
          corrected_PCA,
          corrected_PCA_data,
          color = variable_of_interest,
          shape = batch
        )+
          ggplot2::scale_shape_manual(values = seq(0, length(meta[,variable_of_interest])))
      }
    }
  }
  return(pca_plot_list)
}

