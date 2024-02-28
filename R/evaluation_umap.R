
#' Evaluation umap
#'
#' @param mat A Numeric matrix or list of matrices after pre-processing and/or batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch.1 Column name from the meta file of the column that will be used for batch information
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param color_by Used by evaluation multiple components, pca, and rle to select which feature will be used to color the individuals. Choices are "batch", "variable_of_interest", and "all". Default is set to "all"
#'
#' @return A list object of UMAP plots for each annotation method selected
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_umap = function(mat,
                           meta,
                           batch.1,
                           variable_of_interest,
                           color_by,
                           plot_title){
  evaluation_umap_list <- list()
  #run umaps
  umap_res = umap::umap(t(mat))

  if ("batch" %in% color_by){
    evaluation_umap_list$batch_colored_umap <- dplyr::bind_cols(meta,
                                             as.data.frame(umap_res$layout) %>%
                                               dplyr::rename("x" = 1, "y" = 2)) %>%
      ggplot2::ggplot(aes(x = x, y = y)) +
      geom_point(aes(color = get(batch.1))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(plot_title) +
      theme(axis.text.x = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 18),
            plot.title = element_text(size = 18))
      #labs(title = paste0(x, " Expression"), color = batch.1)
  }
  if ("variable_of_interest" %in% color_by){
    evaluation_umap_list$voi_colored_umap <- dplyr::bind_cols(meta,
                                                              as.data.frame(umap_res$layout) %>%
                                                                dplyr::rename("x" = 1, "y" = 2)) %>%
      ggplot2::ggplot(aes(x = x, y = y)) +
      geom_point(aes(color = get(variable_of_interest))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(plot_title) +
      theme(axis.text.x = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 18),
            plot.title = element_text(size = 18))
      #labs(title = paste0(x, " Expression"), color = variable_of_interest)
  }
  if ("BnW" %in% color_by){
    evaluation_umap_list$BnW_colored_umap <- dplyr::bind_cols(meta,
                                                              as.data.frame(umap_res$layout) %>%
                                                                dplyr::rename("x" = 1, "y" = 2)) %>%
      ggplot2::ggplot(aes(x = x, y = y)) +
      geom_point() +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(plot_title) +
      theme(axis.text.x = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 18),
            plot.title = element_text(size = 18))
    #labs(title = paste0(x, " Expression"), color = variable_of_interest)
  }
  return(evaluation_umap_list)
}
