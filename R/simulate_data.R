#' Simulate Data
#'
#' Generates simulated data based on the linear model framework introduced by Gagnon-Bartsch and Speed.
#'
#' @param BatchFLEX_function A Character vector of BatchFLEX functions in c("retrieve_data", "simulate_data", "preprocess_data", "batch_correct", "batch_evaluate").
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
#'
#' @return A list of a simulated matrix and a simulated meta file
#' @export
#'
#' @examples
#' set.seed(333)
#' simulated_data = simulate_data()
#' head(as.data.frame(simulated_data$sim_matrix), n = c(5,5))
#' head(simulated_data$sim_meta, n = c(5,5))
#' test_simulate_BatchFLEX <- Batch_FLEX(Batch_FLEX_function = c("batch_correct", "batch_evaluate"), correction_method = "ComBat", evaluation_method = "pca", mat = simulated_data$sim_matrix, meta = simulated_data$sim_meta, batch.1 = "batch", color_by = "batch")

simulate_data <- function(num_samples = 100,
                          num_genes = 10000,
                          num_control = 9000,
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
                          epsilon_sd = 1
){
  if (num_treatments != length(treatment_effect_multiplier)){
    stop("The treatment effect multiplier must equal the number of treatments!")
  }
  if (additional_batch != length(batch_effect_multiplier)){
    stop("The batch effect multiplier must equal the number of additional batches!")
  }
  simulate_data_list <- list()
  ctl <- rep(FALSE, num_genes)
  ctl[1:num_control] <- TRUE
  # treatment effect
  X <- vector()
  for (effect in treatment_effect_multiplier) {
    X <- base::append(X, c(rep(effect, num_samples/num_treatments)))
  }
  X <- base::matrix(X, num_samples, 1)
  beta <- base::matrix(rnorm(1*num_genes, treatment_effect, treatment_sd), 1, num_genes) #treatment coefficients
  beta[, ctl] <- 0

  # batch effect
  W <- vector()
  segment_length <- num_samples/num_treatments
  samples_per_batch <- num_samples/num_treatments*batch_proportion/additional_batch
  for (segments in 1:num_treatments){
    for (beffect in batch_effect_multiplier){
      W <- base::append(W, rep(beffect, samples_per_batch))
    }
    W <- base::append(W, rep(0, segment_length - num_samples/num_treatments*batch_proportion))
  }
  W <- base::matrix(W, num_samples, 1)
  alpha <- base::matrix(rnorm(1*num_genes, batch_effect, batch_effect_sd), 1, num_genes)
  Y_alpha <- base::sapply(alpha, function(alpha){rnorm(num_samples, mean =  alpha,
                                                       abs(rnorm(1, mean = batch_sample_sd_mean, sd = batch_sample_sd_sd)))})
  YY_alpha <- base::apply(Y_alpha, 2, function(x){x*W})

  epsilon <- base::matrix(rnorm(num_samples*num_genes, epsilon_mean, epsilon_sd), num_samples, num_genes)
  sim_matrix <- X%*%beta + YY_alpha + epsilon
  sim_matrix <- t(sim_matrix)
  colnames(sim_matrix) <- c(paste0("Sample", 1:num_samples))
  rownames(sim_matrix) <- rownames(BatchFLEX::example_mat)[1:num_genes]
  simulate_data_list$sim_matrix <- sim_matrix
  sim_meta <- as.data.frame(colnames(sim_matrix))
  sim_meta <- as.data.frame(sim_meta)
  sim_meta <- cbind(sim_meta, W)
  sim_meta <- cbind(sim_meta, X)
  names(sim_meta) <- c("sampleID", "batch", "treatment")
  simulate_data_list$sim_meta <- sim_meta
  return(simulate_data_list)
}

