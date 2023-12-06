#' Plot UMAPs
#'
#' @param adjusted_list list of matrices from `batch_correction()`
#' @param meta data frame of the sample level information
#' @param color column in `meta` which to color points - usually batch
#' @param method character string for which batch corrections to plot
#'
#' @import ggplot2
#'
#' @return list of plots including unadjusted data and adjusted data
#' @export
#'
#' @examples
#' set.seed(333)
plot_umaps = function(adjusted_list, meta, color = NULL, method = NULL){
  #checks and things
  if(!is(adjusted_list, "list")) stop('adjusted_list must be a list')
  if(!all(method %in% c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg"))) stop("Batch correction method not found")
  if(!all(method %in% names(adjusted_list))) stop("method not in adjusted_list")
  if(!is.null(color))
    if(length(color) > 1 | !(color %in% colnames(meta))) stop("color either length > 0 or not column in meta data")

  #run umaps
  umap_res = lapply(adjusted_list[c("Unadjusted", method)], function(x){
    umap::umap(t(x))
  })

  plotlist = lapply(c("Unadjusted", method), function(x){
    dplyr::bind_cols(meta,
              as.data.frame(umap_res[[x]]$layout) %>%
                dplyr::rename("x" = 1, "y" = 2)) %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(aes(color = get(color))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = paste0(x, " Expression"),
           color = color)
  })
  names(plotlist) = c("Unadjusted", method)

  return(plotlist)
}
