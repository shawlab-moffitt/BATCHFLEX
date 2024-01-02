#' Evaluation Surrogate Variable Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param batch_correction Numeric matrix following batch correction with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param sva_method Used by SVA to generate surrogate variables. Choices are "be" and "leek". Default is set to "be"
#'
#' @return A list object of uncorrected and batch corrected probability association plots and surrogate variable matrices
#' @export
#'
#' @examples
evaluation_sva <- function(mat,
                           batch_correction,
                           meta,
                           variable_of_interest,
                           sva_method){
  evaluation_sva_list <- list()
  uncorrected_mod <- stats::model.matrix(reformulate(variable_of_interest), data = meta)
  uncorrected_mod_null <- stats::model.matrix(~1, data = meta)
  uncorrected_SVA_nsv <- sva::num.sv(mat, uncorrected_mod, method = sva_method)
  uncorrected_SVA_object <- sva::sva(mat,
                                     uncorrected_mod, uncorrected_mod_null,
                                     n.sv = uncorrected_SVA_nsv,
                                     numSVmethod = sva_method
  )
  uncorrected_SVA_probability_df <- data.frame(Genes = 1:length(uncorrected_SVA_object$pprob.gam),
                                               latent_variable = uncorrected_SVA_object$pprob.gam,
                                               variable_of_intrest = uncorrected_SVA_object$pprob.b
  )
  uncorrected_SVA_probability_df_longer <- pivot_longer(uncorrected_SVA_probability_df, !Genes, names_to = "variable_type")
  uncorrected_SVA_probability_df_longer <- dplyr::rename(uncorrected_SVA_probability_df_longer, "probability_association_of_each_gene" = "value")
  evaluation_sva_list$Plots$uncsvaprob <- ggplot(uncorrected_SVA_probability_df_longer,
                                                  mapping = aes(x = probability_association_of_each_gene, fill = variable_type, alpha = 0.5)
                                                  )+ geom_density()
  evaluation_sva_list$Matrices$uncsv <- uncorrected_SVA_object$sv
  if(!is.null(batch_correction)){
    corrected_mod <- stats::model.matrix(reformulate(variable_of_interest), data = meta)
    corrected_mod_null <- stats::model.matrix(~1, data = meta)
    corrected_SVA_nsv <- sva::num.sv(batch_correction, uncorrected_mod, method = sva_method)
    corrected_SVA_object <- sva::sva(mat,
                                       corrected_mod, uncorrected_mod_null,
                                       n.sv = corrected_SVA_nsv,
                                       numSVmethod = sva_method
    )
    corrected_SVA_probability_df <- data.frame(Genes = 1:length(corrected_SVA_object$pprob.gam),
                                                 latent_variable = corrected_SVA_object$pprob.gam,
                                                 variable_of_intrest = corrected_SVA_object$pprob.b
    )
    corrected_SVA_probability_df_longer <- pivot_longer(corrected_SVA_probability_df, !Genes, names_to = "variable_type")
    corrected_SVA_probability_df_longer <- dplyr::rename(corrected_SVA_probability_df_longer, "probability_association_of_each_gene" = "value")
    evaluation_sva_list$Plots$corsvaprob <- ggplot(corrected_SVA_probability_df_longer,
                                                    mapping = aes(x = probability_association_of_each_gene, fill = variable_type, alpha = 0.5)
    )+ geom_density()
    evaluation_sva_list$Matrices$corsv <- corrected_SVA_object$sv
  }
  return(evaluation_sva_list)
}


