#' Evaluation Cluster Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param cluster_analysis_method Used to select cluster analysis method. Elbow = "wss", Silhouette = "silhouette", Dunn = "dunn', all generates plots from each method
#'
#' @return A list object of uncorrected and batch corrected cluster plots for each selected cluster analysis method
#' @export
#'
#' @examples
evaluation_cluster_analysis = function(mat,
                                       batch_correction,
                                       cluster_analysis_method){
  if("all" %in% cluster_analysis_method){
    cluster_analysis_method = c("wss", "silhouette", "dunn")
  }
  if(!all(cluster_analysis_method %in% c("wss", "silhouette", "dunn"))){
    stop("Cluster analysis method not found")
  }
  cluster_plot_list <- list()
  if ("wss" %in% cluster_analysis_method){
    cluster_plot_list$uncwss <- factoextra::fviz_nbclust(x = t(mat), FUNcluster = kmeans, method = "wss", verbose = T)
    if (!is.null(batch_correction)){
      cluster_plot_list$corwss <- factoextra::fviz_nbclust(x = t(batch_correction), FUNcluster = kmeans, method = "wss", verbose = T)
    }
  }
  if ("silhouette" %in% cluster_analysis_method){
    cluster_plot_list$uncsil <- factoextra::fviz_nbclust(x = t(mat), FUNcluster = kmeans, method = "silhouette", verbose = T)
    if (!is.null(batch_correction)){
      cluster_plot_list$corsil <- factoextra::fviz_nbclust(x = t(batch_correction), FUNcluster = kmeans, method = "silhouette", verbose = T)
    }
  }
  if ("dunn" %in% cluster_analysis_method){
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
    cluster_plot_list$uncdunn <- ggplot(data = uncorrected_dunn_index_analysis, mapping = aes(x = cluster_number, y = dunn_index))+
      geom_point(color = "dodgerblue1")+
      geom_line(color = "dodgerblue1")+
      geom_vline(
        xintercept = uncorrected_dunn_index_analysis$cluster_number[which(max(uncorrected_dunn_index_analysis$dunn_index) == uncorrected_dunn_index_analysis$dunn_index)],
        color = "dodgerblue1",
        linetype = 2
      ) +
      theme_classic()
    if (!is.null(batch_correction)){
      corrected_dunn_k <- c(2:10)
      corrected_dunnin <- c()
      for (i in corrected_dunn_k){
        corrected_dunnin[i] <- clValid::dunn(
          distance = dist(t(mat)),
          clusters = kmeans(t(mat), i)$cluster
        )
      }
      corrected_dunn_index_analysis <- as.data.frame(cbind(corrected_dunn_k, corrected_dunnin[-1]))
      colnames(corrected_dunn_index_analysis) <- c("cluster_number", "dunn_index")
      cluster_plot_list$cordunn <- ggplot(data = corrected_dunn_index_analysis, mapping = aes(x = cluster_number, y = dunn_index))+
        geom_point(color = "dodgerblue1")+
        geom_line(color = "dodgerblue1")+
        geom_vline(
          xintercept = corrected_dunn_index_analysis$cluster_number[which(max(corrected_dunn_index_analysis$dunn_index) == corrected_dunn_index_analysis$dunn_index)],
          color = "dodgerblue1",
          linetype = 2
        ) +
        theme_classic()
    }
  }
  return(cluster_plot_list)
}
