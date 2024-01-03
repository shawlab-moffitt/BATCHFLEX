Batch_FLEX = function(Batch_FLEX_function = NULL,
                      mat = NULL,
                      meta = NULL,
                      correction_method = NULL,
                      batch.1 = NULL,
                      batch.2 = NULL,
                      log2_transformed = TRUE,
                      variable_of_interest = NULL,
                      housekeeping = NULL,
                      k = 2,
                      drop = 0,
                      center = FALSE,
                      round = FALSE,
                      tolerance = 1e-8,
                      par.prior = TRUE,
                      sva_nsv_method = "be",
                      batch_correction = NULL,
                      evaluation_method = NULL,
                      annotation = NULL,
                      uncorrected_cluster_number = NULL,
                      corrected_cluster_number = NULL,
                      cluster_analysis_method = NULL,
                      color_by = "batch",
                      ncomponents = 5,
                      pca_factors = NULL,
                      variable_choices = NULL){
  set.seed(333)
  Batch_FLEX_list <- list()
  if (!is.null(batch_correction)){
    Batch_FLEX_list$batch_correction$User_provided <- batch_correction
  }
  if (is.null(mat)){
    stop("Please provide a matrix file or use retrieve_data or generate_data to generate a matrix file")
  }
  if (is.null(meta)){
    stop("please provide a meta file or use retrieve_data or generate_data to generate a meta file")
  }
  if (is.null(batch.1) & "batch_correct" %in% Batch_FLEX_function | "batch_evaluate" %in% Batch_FLEX_function){
    stop("Please select column for batch information")
  }
  if (is.null(variable_of_interest) & "batch_correct" %in% Batch_FLEX_function | "batch_evaluate" %in% Batch_FLEX_function){
    stop("Please select the variable of interest")
  }
  if (is.null(housekeeping) & "batch_correct" %in% Batch_FLEX_function & "RUVg" %in% correction_method | "all" %in% correction_method){
    stop("Please provide a list of housekeeping genes")
  }
  if (is.null(Batch_FLEX_function)){
    Batch_FLEX_function = c("batch_correct", "batch_evaluate")
  }
  if (!all(Batch_FLEX_function %in% c("retrieve_data", "generate_data", "preprocess_data", "batch_correct", "batch_evaluate"))){
    stop("BatchFLEX function not found")
  }
  if (is.null(correction_method)){
    correction_method = "all"
  }
  if(is.null(evaluation_method)){
    evaluation_method = "all"
  }
  if(is.null(annotation)){
    annotation = "all"
  }
  if(is.null(cluster_analysis_method)){
    cluster_analysis_method = "all"
  }
  if ("batch_correct" %in% Batch_FLEX_function){
    Batch_FLEX_list$batch_correction <- batch_correct(mat, meta, correction_method, batch.1, batch.2, log2_transformed, variable_of_interest, housekeeping,
                                          k, drop, center, round, tolerance, par.prior, sva_nsv_method)
  }
  if ("batch_evaluate" %in% Batch_FLEX_function){
    for (correction in 1:length(Batch_FLEX_list$batch_correction)){
      correction_name <- names(Batch_FLEX_list$batch_correction)[[correction]]
      batch_correction = Batch_FLEX_list$batch_correction[[correction]]
      Batch_FLEX_list$batch_evaluation[[correction_name]]$batch1 <- batch_evaluate(mat, batch_correction, meta, evaluation_method, batch.1, annotation, uncorrected_cluster_number,
                                                         corrected_cluster_number, variable_of_interest, cluster_analysis_method, color_by, ncomponents,
                                                         pca_factors, variable_choices, sva_nsv_method)
      if (!is.null(batch.2)){
        Batch_FLEX_list$batch_evaluation[[correction_name]]$batch2 <- batch_evaluate(mat, batch_correction, meta, evaluation_method, batch.1 = batch.2, annotation, uncorrected_cluster_number,
                                                                  corrected_cluster_number, variable_of_interest, cluster_analysis_method, color_by, ncomponents,
                                                                  pca_factors, variable_choices, sva_nsv_method)
      }
    }
  }
  return(Batch_FLEX_list)
}
