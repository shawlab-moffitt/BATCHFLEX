#' Evaluation Cluster Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param cluster_analysis_method Used to select cluster analysis method. Elbow = "wss", Silhouette = "silhouette", Dunn = "dunn', all generates plots from each method
#'
#' @return A list object of cluster plots for each selected cluster analysis method
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_cluster_analysis = function(mat,
                                       cluster_analysis_method){
  if("all" %in% cluster_analysis_method){
    cluster_analysis_method = c("wss", "silhouette", "dunn")
  }
  if(!all(cluster_analysis_method %in% c("wss", "silhouette", "dunn"))){
    stop("Cluster analysis method not found")
  }
  cluster_plot_list <- list()
  if ("wss" %in% cluster_analysis_method){
    cluster_plot_list$wss <- factoextra::fviz_nbclust(x = t(mat), FUNcluster = kmeans, method = "wss", verbose = T)
  }
  if ("silhouette" %in% cluster_analysis_method){
    cluster_plot_list$silhouette <- factoextra::fviz_nbclust(x = t(mat), FUNcluster = kmeans, method = "silhouette", verbose = T)
  }
  if ("dunn" %in% cluster_analysis_method){
    set.seed(333)
    dunn_k <- c(2:10)
    dunnin <- c()
    for (i in dunn_k){
      dunnin[i] <- clValid::dunn(
        distance = dist(t(mat)),
        clusters = kmeans(t(mat), i)$cluster
      )
    }
    dunn_index_analysis <- as.data.frame(cbind(dunn_k, dunnin[-1]))
    colnames(dunn_index_analysis) <- c("cluster_number", "dunn_index")
    cluster_plot_list$dunn <- ggplot(data = dunn_index_analysis, mapping = aes(x = cluster_number, y = dunn_index))+
      geom_point(color = "dodgerblue1")+
      geom_line(color = "dodgerblue1")+
      geom_vline(
        xintercept = dunn_index_analysis$cluster_number[which(max(dunn_index_analysis$dunn_index) == dunn_index_analysis$dunn_index)],
        color = "dodgerblue1",
        linetype = 2
      ) +
      theme_classic()
  }
  return(cluster_plot_list)
}
