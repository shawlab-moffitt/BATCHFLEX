#' Evaluation Surrogate Variable Analysis
#'
#' @param mat Numeric matrix after pre-processing with features as rownames and sample names as the column names
#' @param meta Data frame of sample data with the first column being sample names that match the column names of the matrix
#' @param variable_of_interest Column name from the meta file of the column that will be used for the variable of interest information
#' @param sva_nsv_method Used by SVA to generate surrogate variables. Choices are "be" and "leek". Default is set to "be"
#'
#' @return A list object of probability association plots and surrogate variable matrices
#' @export
#'
#' @examples
#' set.seed(333)
evaluation_sva <- function(mat,
                           meta,
                           variable_of_interest,
                           sva_nsv_method){
  evaluation_sva_list <- list()
  mod <- stats::model.matrix(reformulate(variable_of_interest), data = meta)
  mod_null <- stats::model.matrix(~1, data = meta)
  SVA_nsv <- sva::num.sv(mat, mod, method = sva_nsv_method)
  SVA_object <- sva::sva(mat,
                         mod, mod_null,
                         n.sv = SVA_nsv,
                         numSVmethod = sva_nsv_method
  )
  SVA_probability_df <- data.frame(Genes = 1:length(SVA_object$pprob.gam),
                                               latent_variable = SVA_object$pprob.gam,
                                               variable_of_intrest = SVA_object$pprob.b
  )
  SVA_probability_df_longer <- pivot_longer(SVA_probability_df, !Genes, names_to = "variable_type")
  SVA_probability_df_longer <- dplyr::rename(SVA_probability_df_longer, "probability_association_of_each_gene" = "value")
  evaluation_sva_list$Plots$svaprob <- ggplot(SVA_probability_df_longer,
                                                  mapping = aes(x = probability_association_of_each_gene, fill = variable_type, alpha = 0.5)
                                                  )+ geom_density()
  evaluation_sva_list$Matrices$sv <- SVA_object$sv

  return(evaluation_sva_list)
}


