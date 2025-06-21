# R/functions.R

#' Install and Load Required R Packages
#'
#' Checks if a list of packages is installed, installs them if they are not,
#' and then loads them into the session.
#'
#' @param packages_to_install A character vector of package names.
#' @return Invisible. This function is called for its side effect of loading packages.
#' @importFrom BiocManager install
install_and_load_packages <- function(packages_to_install) {
  # Ensure BiocManager is installed
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages_to_install) {
    # Check for BioConductor/CRAN packages and install if not present
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
    # Load the package
    library(pkg, character.only = TRUE)
  }
}


#' Create a Base Single-Cell Dataset
#'
#' Generates a synthetic single-cell dataset using the SimBu package structure.
#'
#' @param n_genes An integer, the number of genes.
#' @param n_cells An integer, the number of cells.
#' @return A `SimBu_dataset` S4 object containing synthetic data.
#' @importFrom Matrix Matrix t
#' @importFrom stats rpois
create_base_sc_dataset <- function(n_genes = 1000, n_cells = 300) {
  gene_names <- paste0("gene_", 1:n_genes)
  cell_names <- paste0("cell_", 1:n_cells)
  
  # Create sparse matrices for efficiency
  counts <- Matrix::Matrix(matrix(stats::rpois(n_genes * n_cells, 5), ncol = n_cells), sparse = TRUE)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  tpm <- Matrix::Matrix(matrix(stats::rpois(n_genes * n_cells, 5), ncol = n_cells), sparse = TRUE)
  tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
  rownames(tpm) <- gene_names
  colnames(tpm) <- cell_names
  
  # Create cell annotations
  annotation <- data.frame(
    "ID" = cell_names,
    "cell_type" = c(
      rep("T cells CD4", 50), rep("T cells CD8", 50),
      rep("Macrophages", 100), rep("NK cells", 10),
      rep("B cells", 70), rep("Monocytes", 20)
    )
  )
  
  # Create the SimBu dataset object
  simbu_ds <- SimBu::dataset(
    annotation = annotation,
    count_matrix = counts,
    tpm_matrix = tpm,
    name = "test_dataset"
  )
  return(simbu_ds)
}


#' Generate Bulk RNA-seq Simulations
#'
#' Generates multiple batches of simulated bulk RNA-seq data from a single-cell dataset.
#'
#' @param simbu_dataset A `SimBu_dataset` object.
#' @param n_sims An integer, the number of simulations (batches) to generate.
#' @param n_samples_per_sim An integer, the number of samples in each simulation.
#' @param n_cells_per_sample An integer, the number of cells to pool for each bulk sample.
#' @return A list of `SimBu_simulation` objects.
generate_simulations <- function(simbu_dataset, n_sims = 3, n_samples_per_sim = 30, n_cells_per_sample = 100) {
  sim_list <- lapply(1:n_sims, function(i) {
    SimBu::simulate_bulk(
      data = simbu_dataset,
      scenario = "random",
      scaling_factor = "NONE",
      ncells = n_cells_per_sample,
      nsamples = n_samples_per_sim,
      run_parallel = FALSE
    )
  })
  return(sim_list)
}


#' Merge Simulations and Introduce a Batch Effect
#'
#' Merges a list of simulations and manually introduces a batch effect
#' to the count matrix for demonstration purposes.
#'
#' @param sim_list A list of `SimBu_simulation` objects.
#' @param n_genes_affected An integer, the number of genes to be affected by the batch effect.
#' @param effect_multiplier A numeric value for the batch effect strength.
#' @return A list containing the modified `SummarizedExperiment` object (`se_with_batch`),
#'         a character vector of batch information (`batch_info`), and the names of
#'         the genes affected by the batch effect (`affected_genes`).
#' @importFrom SummarizedExperiment assay `assay<-`
introduce_batch_effect <- function(sim_list, n_genes_affected = 150, effect_multiplier = 2.5) {
  merged_list <- SimBu::merge_simulations(sim_list)
  se_object <- merged_list$bulk
  colnames(se_object) <- make.names(colnames(se_object), unique = TRUE)
  
  counts_matrix <- SummarizedExperiment::assay(se_object, "bulk_counts")
  n_samples_per_batch <- ncol(sim_list[[1]]$bulk)
  batch_info <- rep(paste("Batch", 1:length(sim_list)), each = n_samples_per_batch)
  
  set.seed(456)
  affected_genes <- sample(rownames(counts_matrix), n_genes_affected)
  
  batch2_cols <- which(batch_info == "Batch 2")
  batch3_cols <- which(batch_info == "Batch 3")
  
  counts_matrix[affected_genes, batch2_cols] <- round(counts_matrix[affected_genes, batch2_cols] * effect_multiplier)
  counts_matrix[affected_genes, batch3_cols] <- round(counts_matrix[affected_genes, batch3_cols] / effect_multiplier)
  
  SummarizedExperiment::assay(se_object, "bulk_counts") <- counts_matrix
  
  # *** CHANGE: Return affected_genes as well ***
  return(list(
    se_with_batch = se_object, 
    batch_info = batch_info, 
    affected_genes = affected_genes
  ))
}


#' Perform PCA and Generate a Plot
#'
#' Performs Principal Component Analysis (PCA) on log-transformed counts
#' and creates a ggplot object to visualize the results.
#'
#' @param log_counts A matrix of log-transformed counts (genes in rows, samples in columns).
#' @param batch_info A character or factor vector with batch information for each sample.
#' @param plot_title A string for the plot title.
#' @return A `ggplot` object.
#' @importFrom edgeR cpm
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw scale_color_brewer coord_fixed
perform_pca_and_plot <- function(log_counts, batch_info, plot_title) {
  pca_result <- stats::prcomp(t(log_counts))
  pca_summary <- summary(pca_result)
  
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    batch = batch_info
  )
  
  pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
  pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)
  
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = batch)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::labs(
      title = plot_title,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::coord_fixed()
  
  return(p)
}

#' Apply ComBat-Seq for Batch Correction
#'
#' Applies the ComBat-Seq algorithm from the sva package to correct for batch
#' effects in raw count data. This wrapper handles the conversion from a potential
#' sparse matrix to the dense matrix required by ComBat-Seq.
#'
#' @param raw_counts A matrix of non-normalized, integer counts (genes in rows, samples in columns).
#'   Can be a standard or sparse matrix.
#' @param batch_info A character or factor vector with batch information for each sample.
#' @return A matrix of batch-corrected counts.
#' @importFrom sva ComBat_seq
run_combat_seq <- function(raw_counts, batch_info) {
  # *** THE FIX: ComBat-Seq requires a standard dense matrix. ***
  # Coerce the input to a standard matrix to ensure compatibility.
  counts_matrix <- as.matrix(raw_counts)
  
  # ComBat-Seq requires a factor for the batch argument.
  batch_factor <- as.factor(batch_info)
  
  # ComBat-Seq works on raw, integer counts.
  # It returns a matrix of corrected counts.
  corrected_counts <- sva::ComBat_seq(counts = counts_matrix, batch = batch_factor)
  
  return(corrected_counts)
}

#' Apply original ComBat for Batch Correction
#'
#' Applies the original ComBat algorithm, which expects log-transformed,
#' normalized data (e.g., from microarrays or log-CPM from RNA-seq).
#'
#' @param log_counts A matrix of log-transformed counts.
#' @param batch_info A character or factor vector with batch information.
#' @return A matrix of batch-corrected log-transformed counts.
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
run_combat <- function(log_counts, batch_info) {
  # ComBat expects a standard dense matrix.
  log_counts_matrix <- as.matrix(log_counts)
  batch_factor <- as.factor(batch_info)
  
  # Use an intercept-only model as there are no covariates to preserve.
  # This is the 'mod' argument in ComBat.
  mod_combat <- stats::model.matrix(~1, data = data.frame(batch = batch_factor))
  
  # Run ComBat on the log-transformed data
  corrected_log_counts <- sva::ComBat(
    dat = log_counts_matrix,
    batch = batch_factor,
    mod = mod_combat
  )
  
  return(corrected_log_counts)
}


#' Apply RUVg from RUVSeq for Batch Correction
#'
#' Applies the RUVg algorithm, which uses control genes to estimate and remove
#' unwanted variation. It operates on raw, integer counts.
#'
#' @param raw_counts A matrix of non-normalized, integer counts.
#' @param control_genes A character vector of gene names to be used as controls.
#' @param k An integer, the number of factors of unwanted variation to remove.
#' @return A matrix of RUVg-normalized counts.
#' @importFrom RUVSeq RUVg
#' @importFrom EDASeq newSeqExpressionSet
#' @importFrom BiocGenerics counts
run_ruvg <- function(raw_counts, control_genes, k = 1) {
  # RUVg requires a SeqExpressionSet object, which is created by a function
  # in the EDASeq package.
  counts_matrix <- round(as.matrix(raw_counts))
  
  # Filter for control genes that are present in the count matrix
  valid_control_genes <- intersect(control_genes, rownames(counts_matrix))
  if (length(valid_control_genes) == 0) {
    stop("None of the provided control genes were found in the count matrix.")
  }
  
  # *** THE FINAL FIX: Call newSeqExpressionSet from its source package, EDASeq. ***
  set <- EDASeq::newSeqExpressionSet(counts = counts_matrix)
  
  # Run RUVg using the control genes to estimate factors of unwanted variation
  set_normalized <- RUVSeq::RUVg(set, valid_control_genes, k = k)
  
  # Return the normalized count matrix
  return(counts(set_normalized))
}
