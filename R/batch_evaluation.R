batch_evaluation = function(mat = NULL,
                            batch_correction = NULL,
                            meta = NULL,
                            evaluation_method,
                            batch = "Study",
                            annotation = NULL,
                            uncorrected_cluster_number = NULL,
                            corrected_cluster_number = NULL,
                            variable_of_interest = "CellType",
                            cluster_analysis = NULL,
                            color_by = "batch"){
  batch_evaluation_list <- list("Graphs", "Matrices")
  if ("pca" %in% evaluation_method){
    cat("\tConducting PCA analysis\n")
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
  if("cluster_analysis" %in% evaluation_method){
    cat("\tConducting cluster analysis\n")
    batch_evaluation_list$Graphs$cluster_analysis <- evaluation_cluster_analysis(mat,
                                                                          batch_correction,
                                                                          cluster_analysis)
  }
  return(batch_evaluation_list)
}
