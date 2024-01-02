#' Evaluation PCA Details
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param batch Column name from the meta file of the column that will be used for batch information
#' @param pca_factors Column name from the meta file of the column that will be used to group the summary details of the PCA plot. If NULL, batch is selected
#'
#' @return A list object of uncorrected and batch corrected PCA details including a scree plot, a matrix of pca components, a matrix of contributions grouped by factor, and a matrix of contribution counts
#' @export
#'
#' @examples
evaluation_pca_details <- function(mat,
                                   meta,
                                   batch_correction,
                                   batch,
                                   pca_factors
){
  if (if.null(pca_factors)){
    pca_factors <- batch
  }

  pca_details_list <- list()

  uncorrected_PCA_details <- FactoMineR::PCA(t(mat), graph = F)

  # uncorrected components matrix
  pca_details_list$Matrices$unccomponents <- uncorrected_PCA_details$x

  #uncorrected scree plot
  uncorrected_scree_eig <- as.data.frame(factoextra::get_eig(uncorrected_PCA_details))
  uncorrected_scree_eig$Dimensions <- gsub("Dim\\.","", rownames(uncorrected_scree_eig))
  uncorrected_scree_eig$`Variance Percent` <- paste0(round(uncorrected_scree_eig$variance.percent,1),"%")
  uncorrected_scree_eig_top10 <- uncorrected_scree_eig[1:10,]
  pca_details_list$Plots$uncscree <- ggplot(data=uncorrected_scree_eig_top10, aes(x=reorder(Dimensions,-variance.percent), y=variance.percent)) +
    geom_bar(stat="identity", fill="steelblue")+
    theme_minimal() +
    geom_text(aes(label=`Variance Percent`), vjust=-0.3,size=4.5) +
    labs(x = "Dimensions",y = "Variance Percent") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))

  #uncorrected contributions matrix
  uncorrected_PCA_individuals <- uncorrected_PCA_details
  uncorrected_PCA_factors <- as.data.frame(uncorrected_PCA_individuals$ind$contrib)
  uncorrected_PCA_factors$factors <- as.vector(meta[,pca_factors])
  uncorrected_PCA_factors_final_sum <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, sum)
  uncorrected_PCA_factors_final_mean <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, mean)
  uncorrected_PCA_factors_final_sdv <- stats::aggregate(. ~ factors, uncorrected_PCA_factors, sd)
  uncorrected_PCA_factors_final_longer_sum <- tidyr::pivot_longer(uncorrected_PCA_factors_final_sum, !factors, names_to = "PC_Components")
  uncorrected_PCA_factors_final_longer_mean <- tidyr::pivot_longer(uncorrected_PCA_factors_final_mean, !factors, names_to = "PC_Components")
  uncorrected_PCA_factors_final_longer_sdv <- tidyr::pivot_longer(uncorrected_PCA_factors_final_sdv, !factors, names_to = "PC_Components")
  colnames(uncorrected_PCA_factors_final_longer_sum)[3] <- "Sum"
  colnames(uncorrected_PCA_factors_final_longer_mean)[3] <- "Mean"
  colnames(uncorrected_PCA_factors_final_longer_sdv)[3] <- "SDV"
  uncorrected_PCA_factors_final_longer <- uncorrected_PCA_factors_final_longer_sum[,1:2]
  uncorrected_PCA_factors_final_longer <- cbind(
    uncorrected_PCA_factors_final_longer,
    uncorrected_PCA_factors_final_longer_sum[,3],
    uncorrected_PCA_factors_final_longer_mean[,3],
    uncorrected_PCA_factors_final_longer_sdv[,3]
  )
  colnames(uncorrected_PCA_factors_final_longer)[1] <- pca_factors
  pca_details_list$Matrices$unccontribution <- uncorrected_PCA_factors_final_longer

  #uncorrected counts plot
  uncorrected_PCA_factors_final_count <- uncorrected_PCA_factors %>% dplyr::group_by(factors) %>% dplyr::summarise(individuals = length(factors))
  colnames(uncorrected_PCA_factors_final_count)[1] <- pca_factors
  pca_details_list$Matrices$unccount <- uncorrected_PCA_factors_final_count


  if(!is.null(batch_correction)){
    corrected_PCA_details <- FactoMineR::PCA(t(batch_correction), graph = F)

    #corrected components matrix
    pca_details_list$Matrices$corcomponents <- corrected_PCA_details$x

    #corrected scree plot
    corrected_scree_eig <- as.data.frame(factoextra::get_eig(corrected_PCA_details))
    corrected_scree_eig$Dimensions <- gsub("Dim\\.","", rownames(corrected_scree_eig))
    corrected_scree_eig$`Variance Percent` <- paste0(round(corrected_scree_eig$variance.percent,1),"%")
    corrected_scree_eig_top10 <- corrected_scree_eig[1:10,]
    pca_details_list$Plots$corscree <- ggplot(data=corrected_scree_eig_top10, aes(x=reorder(Dimensions,-variance.percent), y=variance.percent)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_minimal() +
      geom_text(aes(label=`Variance Percent`), vjust=-0.3,size=4.5) +
      labs(x = "Dimensions",y = "Variance Percent") +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16))

    #corrected contributions matrix
    corrected_PCA_individuals <- corrected_PCA_details
    corrected_PCA_factors <- as.data.frame(corrected_PCA_individuals$ind$contrib)
    corrected_PCA_factors$factors <- as.vector(meta[,pca_factors])
    corrected_PCA_factors_final_sum <- stats::aggregate(. ~ factors, corrected_PCA_factors, sum)
    corrected_PCA_factors_final_mean <- stats::aggregate(. ~ factors, corrected_PCA_factors, mean)
    corrected_PCA_factors_final_sdv <- stats::aggregate(. ~ factors, corrected_PCA_factors, sd)
    corrected_PCA_factors_final_longer_sum <- tidyr::pivot_longer(corrected_PCA_factors_final_sum, !factors, names_to = "PC_Components")
    corrected_PCA_factors_final_longer_mean <- tidyr::pivot_longer(corrected_PCA_factors_final_mean, !factors, names_to = "PC_Components")
    corrected_PCA_factors_final_longer_sdv <- tidyr::pivot_longer(corrected_PCA_factors_final_sdv, !factors, names_to = "PC_Components")
    colnames(corrected_PCA_factors_final_longer_sum)[3] <- "Sum"
    colnames(corrected_PCA_factors_final_longer_mean)[3] <- "Mean"
    colnames(corrected_PCA_factors_final_longer_sdv)[3] <- "SDV"
    corrected_PCA_factors_final_longer <- corrected_PCA_factors_final_longer_sum[,1:2]
    corrected_PCA_factors_final_longer <- cbind(
      corrected_PCA_factors_final_longer,
      corrected_PCA_factors_final_longer_sum[,3],
      corrected_PCA_factors_final_longer_mean[,3],
      corrected_PCA_factors_final_longer_sdv[,3]
    )
    colnames(corrected_PCA_factors_final_longer)[1] <- pca_factors
    pca_details_list$Matrices$corcontribution <- corrected_PCA_factors_final_longer

    #corrected counts plot
    corrected_PCA_factors_final_count <- corrected_PCA_factors %>% dplyr::group_by(factors) %>% dplyr::summarise(individuals = length(factors))
    colnames(corrected_PCA_factors_final_count)[1] <- pca_factors
    pca_details_list$matrices$corcount <- corrected_PCA_factors_final_count
  }
  return(pca_details_list)
}
