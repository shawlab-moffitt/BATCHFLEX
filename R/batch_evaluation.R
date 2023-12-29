batch_evaluation = function(mat = NULL,
                            rawmat = NULL,
                            batch_correction = NULL,
                            meta = NULL,
                            evaluation_method,
                            batch = NULL,
                            annotation = NULL,
                            uncorrected_cluster_number = NULL,
                            corrected_cluster_number = NULL,
                            variable_of_interest = NULL,
                            cluster_analysis_method = NULL,
                            color_by = "batch",
                            ncomponents = 5,
                            pca_factors = NULL,
                            variable_choices = NULL,
                            sva_method = "be"){

  batch_evaluation_list <- list()
  if("all" %in% evaluation_method) method = c("pca", "cluster_analysis", "mc_pca", "pca_details", "rle", "ev", "sva")
  if (!all(evaluation_method %in% c("pca", "cluster_analysis", "mc_pca", "pca_details", "rle", "ev", "sva"))) stop("Evaluation correction method not found")
  if ("pca" %in% evaluation_method){
    cat("\tConducting principal component analysis\n")
    batch_evaluation_list$Graphs$pca <-  evaluation_pca(mat,
                                                 batch_correction,
                                                 meta,
                                                 annotation,
                                                 uncorrected_cluster_number,
                                                 corrected_cluster_number,
                                                 batch,
                                                 variable_of_interest,
                                                 color_by)
  }
  if ("cluster_analysis" %in% evaluation_method){
    cat("\tConducting cluster analysis\n")
    batch_evaluation_list$Graphs$cluster_analysis <- evaluation_cluster_analysis(mat,
                                                                          batch_correction,
                                                                          cluster_analysis_method)
  }
  if ("mc_pca" %in% evaluation_method){
    cat("\tConducting multiple components PCA analysis\n")
    batch_evaluation_list$Graphs$mc_pca <- evaluation_mc_pca(mat,
                                                             rawmat,
                                                             batch_correction,
                                                             meta,
                                                             batch,
                                                             variable_of_interest,
                                                             ncomponents,
                                                             color_by)
  }
  if ("pca_details" %in% evaluation_method){
    cat("\tConducting PCA details analysis\n")
    evaluation_pca_details_list <- evaluation_pca_details(mat,
                                                          meta,
                                                          batch_correction,
                                                          pca_factors)
    batch_evaluation_list$Graphs$pca_details <- evaluation_pca_details_list$Graphs
    batch_evaluation_list$Matrices$pca_details <- evaluation_pca_details_list$Matrices
  }
  if ("rle" %in% evaluation_method){
    cat("\tConducting relative log expression analysis\n")
    batch_evaluation_list$Graphs$rle <- evaluation_rle(mat,
                                                       rawmat,
                                                       batch_correction,
                                                       meta,
                                                       batch,
                                                       variable_of_interest,
                                                       color_by)
  }
  if ("ev" %in% evaluation_method){
    cat("\tConducting explanatory variables analysis\n")
    batch_evaluation_list$Graphs$ev <- evaluation_ev(mat,
                                                     rawmat,
                                                     batch_correction,
                                                     meta,
                                                     variable_choices,
                                                     batch,
                                                     variable_of_interest)
  }
  if ("sva" %in% evaluation_method){
    cat("\tConducting surrogate variable analysis\n")
    evaluation_sva_list <- evaluation_sva(mat,
                                          batch_correction,
                                          meta,
                                          variable_of_interest,
                                          sva_method)
    batch_evaluation_list$Graphs$sva <- evaluation_sva_list$Graphs
    batch_evaluation_list$Matrices$sva <- evaluation_sva_list$Matrices
  }
  return(batch_evaluation_list)
}
