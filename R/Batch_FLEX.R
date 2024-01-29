#' BatchFLEX
#'
#' @param Batch_FLEX_function A Character vector of BatchFLEX functions in c("retrieve_data", "merge_data", "simulate_data", "preprocess_data", "batch_correct", "batch_evaluate").
#' @param num_samples An integer setting the number of samples to simulate.
#' @param num_genes An integer setting the number of genes to simulate.
#' @param num_control An integer setting the number of genes that do not respond to treatment.
#' @param num_treatments An integer setting the number of total treatments.
#' @param treatment_effect_multiplier A number or numbers setting how each treatment impacts the treatment effect. Must equal the number of treatments.
#' @param treatment_effect A number setting average impact of the treatments.
#' @param treatment_sd A number setting the standard deviation of the treatments.
#' @param additional_batch An integer setting the number of additional batches. e.g. setting to 1 indicates 2 total batches in the simulation.
#' @param batch_proportion A number setting the proportion of the total sample that is a different batch. Each additional batch is given an equal proportion of the total batch proportion.
#' @param batch_effect_multiplier A number or numbers setting how each batch impacts the treatment effect.
#' @param batch_effect A number setting the initial expression of each gene for a single sample for the batch effect. e.g. a single column of num_genes length that will be used to set the average expression for that gene across a set of batch samples
#' @param batch_effect_sd A number setting the deviation for the initial expression of each gene for a single sample for the batch effect.
#' @param batch_sample_sd_mean A number setting average for the standard deviation for each sample derived from the initial gene expression.
#' @param batch_sample_sd_sd A number setting the deviation for the average standard deviation used in the batch effect calculations.
#' @param epsilon_mean A number setting the average noise level
#' @param epsilon_sd A number setting the deviation for the average noise level
#' @param merge_matrix_files A list of matrix files to be merged. If supplied as a named list, batchflex will add a named column to the meta file with the specified names. If not, batch flex will add a generic name per meta file.
#' @param merge_meta_files A list of meta files to be merged.
#' @param keep_all_genes A TRUE/FALSE statement to indicate how unmatched genes will be handled. If TRUE, then NAs will be inserted for unmatched genes. If FALSE, unmatched genes will be deleted. Default is set to FALSE.
#' @param mat A Numeric matrix or list of matrices after pre-processing and/or batch correction with features as rownames and sample names as the column names.
#' @param raw.counts Logical. TRUE indicates the input data is raw counts. FALSE indicated the input data is not raw counts.
#' @param raw.norm.method Character string of what method to use for raw count normalization. Supported methods are "TMM" or "upperquartile". If raw.counts is FALSE this will be ignored.
#' @param log2 Logical. If TRUE, the input data with be logged with the method "log2+1".
#' @param quantnorm Logical. If TRUE, the input data will be quantile normalized.
#' @param remove.duplicates Logical. If TRUE, if duplicate row features are found, they will be summarized to the row with the maximum average feature value.
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix.
#' @param correction_method A character vector of batch correction methods in c("Limma", "ComBat", "Mean Centering", "ComBatseq", "Harman", "RUVg", "SVA). Default is set to "all", which excludes ComBatseq as it requires a count matrix.
#' @param batch.1 Column name from the meta file of the column that will be used for batch information.
#' @param batch.2 Column name from the meta file of the column that will be used for batch two information.
#' @param log2_transformed logical whether the data is already transformed. Default is set to TRUE.
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information.
#' @param housekeeping Name of housekeeping gene set or character vector of housekeeping genes. Default is set to the preloaded HSIAO human gene list.
#' @param k Used in the RUVg correction_method, the number of factors of unwanted variation to be estimated from the data.
#' @param drop Used in the RUVg correction_method, the number of singular values to drop in the estimation of the factors of unwanted variation. This number is usually zero, but might be set to one if the first singular value captures the effect of interest. It must be less than k.
#' @param center Used in the RUVg correction_method, if TRUE, the counts are centered, for each gene, to have mean zero across samples. This is important to ensure that the first singular value does not capture the average gene expression.
#' @param round Used in the RUVg correction_method, if TRUE, the normalized measures are rounded to form pseudo-counts.
#' @param tolerance Used in the RUVg correction_method, tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than tolerance to be considered positive.
#' @param par.prior Used in the ComBat correction_method, TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used.
#' @param sva_nsv_method Input correction_method for the num.sv function in sva. Default is set to "be", but can be manually set to "leek".
#' @param evaluation_method A character vector of batch correction methods in c("pca", "cluster_analysis", "mc_pca", "pca_details", "rle", "ev", "sva", "umap").
#' @param annotation Used by evaluation pca to select whether a cluster or meta annotated PCA plot is generated. cluster = "cluster", meta = "meta", all = c("cluster", "meta").
#' @param cluster_number Used by evaluation pca to select the number of kmeans generated clusters to display in the uncorrected plot. If NULL is selected, a Dunn generated cluster number is used..
#' @param cluster_analysis_method Used to select cluster analysis method. Elbow = "wss", Silhouette = "silhouette", Dunn = "dunn', all generates plots from each method.
#' @param color_by Used by evaluation multiple components, pca, and rle to select which feature will be used to color the individuals. Choices are "batch", "variable_of_interest", "BnW", and "all". Default is set to "all".
#' @param ncomponents Used by evaluation multiple components to select the number of principal components that will be plotted. Default is set to 5.
#' @param pca_factors Used by evaluation pca details to select the Column name from the meta file of the column that will be used to group the summary details of the PCA plot. If NULL, batch is selected.
#' @param variable_choices Used by the explanatory variables function to select which variables to plot. Default is a combination of the batch and variable of interest.
#'
#' @return A list object containing all plots and matrices generated by BatchFLEX.
#' @export
#'
#' @examples
#' set.seed(333)
Batch_FLEX = function(Batch_FLEX_function = c("batch_correct", "batch_evaluate"),
                      num_samples = 1000,
                      num_genes = 20000,
                      num_control = 18000,
                      num_treatments = 2,
                      treatment_effect_multiplier = c(0, 1),
                      treatment_effect = 5,
                      treatment_sd = 1,
                      additional_batch = 1,
                      batch_proportion = 0.5,
                      batch_effect_multiplier = c(1),
                      batch_effect = 5,
                      batch_effect_sd = 1,
                      batch_sample_sd_mean = 0,
                      batch_sample_sd_sd = 2,
                      epsilon_mean = 0,
                      epsilon_sd = 1,
                      merge_matrix_files = NULL,
                      merge_meta_files = NULL,
                      keep_all_genes = FALSE,
                      mat = NULL,
                      raw.counts = FALSE,
                      raw.norm.method = NULL,
                      log2 = TRUE,
                      quantnorm = TRUE,
                      remove.duplicates = TRUE,
                      meta = NULL,
                      correction_method = "all",
                      batch.1 = NULL,
                      batch.2 = NULL,
                      log2_transformed = TRUE,
                      variable_of_interest = NULL,
                      housekeeping = BatchFLEX::hsiao_mouse,
                      k = 2,
                      drop = 0,
                      center = FALSE,
                      round = FALSE,
                      tolerance = 1e-8,
                      par.prior = TRUE,
                      sva_nsv_method = "be",
                      evaluation_method = "all",
                      annotation = "all",
                      cluster_number = NULL,
                      cluster_analysis_method = "all",
                      color_by = "all",
                      ncomponents = 5,
                      pca_factors = NULL,
                      variable_choices = NULL,
                      plot_title = NULL,
                      large_list = NULL){
  Batch_FLEX_list <- list()
  if ("simulate_data" %in% Batch_FLEX_function | "merge_data" %in% Batch_FLEX_function){
    message("Generating mat from a BatchFLEX function")
  }else {
    if (is.null(mat)){
      stop("Please provide a matrix file or use retrieve_data or simulate_data to generate a matrix file")
    }
  }
  if ("simulate_data" %in% Batch_FLEX_function | "merge_data" %in% Batch_FLEX_function){
    message("Generating meta from a BatchFLEX function")
  }else {
    if (is.null(meta)){
      stop("Please provide a meta file or use retrieve_data or simulate_data to generate a matrix file")
    }
  }
  if (is.null(batch.1) & "batch_correct" %in% Batch_FLEX_function & !"simulate_data" %in% Batch_FLEX_function | is.null(batch.1) & "batch_evaluate" %in% Batch_FLEX_function & !"simulate_data" %in% Batch_FLEX_function){
    stop("Please select column name in the meta file for the batch information")
  }
  if (is.null(variable_of_interest) & "batch_evaluate" %in% Batch_FLEX_function & !"simulate_data" %in% Batch_FLEX_function | is.null(variable_of_interest) & "batch_correct" %in% Batch_FLEX_function & !"simulate_data" %in% Batch_FLEX_function){
    message("Missing variable of interest. Some correction methods and evaluation techniques are not available.")
  }
  if (is.null(housekeeping) & "batch_correct" %in% Batch_FLEX_function & "RUVg" %in% correction_method | is.null(housekeeping) & "batch_correct" %in% Batch_FLEX_function & "all" %in% correction_method){
    stop("Please provide a list of housekeeping genes for the RUVg correction method")
  }
  if (!all(Batch_FLEX_function %in% c("retrieve_data", "simulate_data", "merge_data", "preprocess_matrix", "batch_correct", "batch_evaluate", "BatchFLEX_export"))){
    stop("BatchFLEX function not found")
  }
  if ("simulate_data" %in% Batch_FLEX_function){
    Batch_FLEX_list$simulate_data <- simulate_data(num_samples, num_genes, num_control, num_treatments, treatment_effect_multiplier, treatment_effect,
                                                   treatment_sd, additional_batch, batch_proportion, batch_effect_multiplier, batch_effect, batch_effect_sd,
                                                   batch_sample_sd_mean, batch_sample_sd_sd, epsilon_mean, epsilon_sd)
  }
  if ("simulate_data" %in% Batch_FLEX_function & "batch_correct" %in% Batch_FLEX_function | "simulate_data" %in% Batch_FLEX_function & "batch_evaluate" %in% Batch_FLEX_function){
    mat <- Batch_FLEX_list$simulate_data$sim_matrix
    meta <- Batch_FLEX_list$simulate_data$sim_meta
    batch.1 = "batch"
    variable_of_interest = "treatment"
  }
  if ("merge_data" %in% Batch_FLEX_function){
    Batch_FLEX_list$merge_data <- merge_data(merge_matrix_files, merge_meta_files, keep_all_genes)
  }
  if ("merge_data" %in% Batch_FLEX_function & "batch_correct" %in% Batch_FLEX_function | "merge_data" %in% Batch_FLEX_function & "batch_evaluate" %in% Batch_FLEX_function | "merge_data" %in% Batch_FLEX_function & "preprocess_data" %in% Batch_FLEX_function){
    mat <- Batch_FLEX_list$merge_data$merged_mat
    meta <- Batch_FLEX_list$merge_data$merged_meta
  }
  if (!all(apply(mat,2,is.numeric)) | !is(mat,"matrix")) stop("Must be numeric matrix")
  if ("preprocess_matrix" %in% Batch_FLEX_function){
    cat("\tPre-processing input matix\n")
    mat <- preprocess_matrix(mat, raw.counts, raw.norm.method, log2, quantnorm, remove.duplicates)
    Batch_FLEX_list$data_matrices[[paste0("Unadjusted_", ifelse(log2, "Log2", ""), ifelse(quantnorm, "_Norm", ""))]] <-  as.matrix(mat)
  }
  if ("batch_correct" %in% Batch_FLEX_function){
    if (!"preprocess_matrix" %in% Batch_FLEX_function){
      Batch_FLEX_list$data_matrices$Unadjusted <- mat
    }
    Batch_FLEX_list$data_matrices <- base::append(Batch_FLEX_list$data_matrices, batch_correct(mat, meta, correction_method, batch.1, batch.2, log2_transformed, variable_of_interest, housekeeping,
                                          k, drop, center, round, tolerance, par.prior, sva_nsv_method))
  }
  if (!"batch_correct" %in% Batch_FLEX_function){
    if (is.matrix(mat) & !"preprocess_matrix" %in% Batch_FLEX_function){
      Batch_FLEX_list$data_matrices$Unadjusted <- mat
    }
    if (is.list(mat)){
      Batch_FLEX_list$data_matrices <- mat
    }
  }
  if (is.null(plot_title)){
    plot_titles <- list()
    for (name in names(Batch_FLEX_list$data_matrices)){
      if(is.null(batch.2)){
        if (grepl("Unadjusted", name)){
          plot_titles[[name]] <- name
        } else {
          plot_titles[[name]] <- paste0(name, "_", "batch_1_is_", batch.1, "_", "voi_is_", if (is.null(variable_of_interest)) {"NULL"} else {variable_of_interest}, sep = "")
        }
      } else {
        if (grepl("Unadjusted", name)){
          plot_titles[[name]] <- name
        } else {
          plot_titles[[name]] <- paste0(name, "_", "batch_1_is_", batch.1, "_", "batch_2_is_", batch.2, if (is.null(variable_of_interest)) {"NULL"} else {variable_of_interest}, sep = "")
        }
      }
    }

  }else {
    plot_titles <- plot_title
  }
  if ("batch_evaluate" %in% Batch_FLEX_function){
    for (matrix in 1:length(Batch_FLEX_list$data_matrices)){
      plot_title <- plot_titles[[matrix]]
      matrix_name <- names(Batch_FLEX_list$data_matrices)[[matrix]]
      all_matrices = Batch_FLEX_list$data_matrices[[matrix]]
      Batch_FLEX_list$batch_evaluation[[matrix_name]]$batch1 <- batch_evaluate(mat = all_matrices, meta, evaluation_method, batch.1, annotation, cluster_number,
                                                                               variable_of_interest, cluster_analysis_method, color_by, ncomponents,
                                                                               pca_factors, variable_choices, sva_nsv_method, plot_title)
      if (!is.null(batch.2)){
        Batch_FLEX_list$batch_evaluation[[matrix_name]]$batch2 <- batch_evaluate(mat = all_matrices, meta, evaluation_method, batch.1 = batch.2, annotation, cluster_number,
                                                                                     variable_of_interest, cluster_analysis_method, color_by, ncomponents,
                                                                                     pca_factors, variable_choices, sva_nsv_method, plot_title)
      }
    }
  }
  if (is.null(large_list)){
    large_list <- Batch_FLEX_list
  }
  if ("BatchFLEX_export" %in% Batch_FLEX_function){
    BatchFLEX_export(large_list)
  }
  return(Batch_FLEX_list)
}
