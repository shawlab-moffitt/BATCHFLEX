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
  evaluation_sva_list$Graphs$uncsvaprob <- ggplot(uncorrected_SVA_probability_df_longer,
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
    evaluation_sva_list$Graphs$corsvaprob <- ggplot(corrected_SVA_probability_df_longer,
                                                    mapping = aes(x = probability_association_of_each_gene, fill = variable_type, alpha = 0.5)
    )+ geom_density()
    evaluation_sva_list$Matrices$corsv <- corrected_SVA_object$sv
  }
  return(evaluation_sva_list)
}


