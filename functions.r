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

#' Generate Faceted Box Plots of Log2-CPM Values
#'
#' Creates a ggplot object showing box plots of Log2-CPM distributions,
#' faceted by correction method and coloured by batch.
#'
#' @param tidy_logcpm_df A tidy data frame with columns: Method, Batch, and Log2CPM.
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_boxplot facet_wrap labs theme_bw theme element_text
generate_log2cpm_boxplot <- function(tidy_logcpm_df) {
  
  # Ensure the order of methods is consistent for plotting
  method_order <- c("Pre-Correction", "Limma", "ComBat", "ComBat-Seq", "RUVg", "PCA Correction")
  plot_df <- tidy_logcpm_df
  plot_df$Method <- factor(plot_df$Method, levels = intersect(method_order, unique(plot_df$Method)))
  
  p <- ggplot(plot_df, aes(x = Batch, y = Log2CPM, fill = Batch)) +
    # Use outlier.shape = NA to avoid a messy plot with tens of thousands of points
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~Method, scales = "free_y") +
    labs(
      title = "Distribution of Log2-CPM Values by Batch and Method",
      subtitle = "Median-centered distributions indicate successful batch effect removal",
      x = "Batch",
      y = "Log2-CPM Value"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "none" # The fill is redundant with the x-axis
    )
  
  return(p)
}

#' Perform PCA and Generate a Plot
#'
#' Performs PCA on log-counts OR uses pre-computed PCA data to create a ggplot object.
#'
#' @param log_counts A matrix of log-transformed counts. Required if `pca_data` is NULL.
#' @param batch_info A character or factor vector with batch information.
#' @param plot_title A string for the plot title.
#' @param pca_data A data frame with columns PC1, PC2, and batch. If provided, `log_counts` is ignored.
#' @return A `ggplot` object.
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw scale_color_brewer coord_fixed
perform_pca_and_plot <- function(log_counts = NULL, batch_info, plot_title, pca_data = NULL) {
  if (is.null(pca_data)) {
    # If no pre-computed PCs are given, calculate them.
    pca_result <- stats::prcomp(t(log_counts))
    pca_summary <- summary(pca_result)
    
    pc1_var <- round(pca_summary$importance[2, 1] * 100, 1)
    pc2_var <- round(pca_summary$importance[2, 2] * 100, 1)
    
    plot_data <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      batch = batch_info
    )
    
    x_lab <- paste0("PC1 (", pc1_var, "%)")
    y_lab <- paste0("PC2 (", pc2_var, "%)")
  } else {
    # Use the pre-computed PCs.
    plot_data <- pca_data
    x_lab <- "PC1" # We don't know the variance explained, so use generic labels
    y_lab <- "PC2"
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = batch)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(title = plot_title, x = x_lab, y = y_lab) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::coord_fixed()
  
  return(p)
}

#' Apply ComBat-ref for Batch Correction using a Reference Batch
#'
#' Applies the ComBat-ref algorithm, which adjusts data to a specified
#' reference batch. Expects log-transformed, normalized data.
#'
#' @param log_counts A matrix of log-transformed counts.
#' @param batch_info A character or factor vector with batch information.
#' @param ref_batch A string specifying the name of the reference batch.
#' @return A matrix of batch-corrected log-transformed counts.
#' @importFrom Combat_ref Combat_ref
run_combat_ref <- function(log_counts, batch_info, ref_batch) {
  log_counts_matrix <- as.matrix(log_counts)
  
  # Combat_ref expects a factor for the batch argument.
  batch_factor <- as.factor(batch_info)
  
  # Run Combat_ref on the log-transformed data
  corrected_log_counts <- Combat_ref::Combat_ref(
    dat = log_counts_matrix,
    batch = batch_factor,
    ref.batch = ref_batch
  )
  
  return(corrected_log_counts)
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
  
  # *** FIX: Call newSeqExpressionSet from its source package, EDASeq. ***
  set <- EDASeq::newSeqExpressionSet(counts = counts_matrix)
  
  # Run RUVg using the control genes to estimate factors of unwanted variation
  set_normalized <- RUVSeq::RUVg(set, valid_control_genes, k = k)
  
  # Return the normalized count matrix
  return(counts(set_normalized))
}

#' Apply fastMNN for Batch Correction
#'
#' Applies the Mutual Nearest Neighbors (MNN) correction from the batchelor package.
#' It returns the corrected principal components, not a corrected count matrix.
#'
#' @param log_counts A matrix of log-transformed counts.
#' @param batch_info A character or factor vector with batch information.
#' @return A data frame containing the corrected principal components for each sample.
#' @importFrom batchelor fastMNN
#' @importFrom SingleCellExperiment reducedDim
run_fastmnn <- function(log_counts, batch_info) {
  # fastMNN returns a SingleCellExperiment object
  sce_corrected <- batchelor::fastMNN(log_counts, batch = as.factor(batch_info))
  
  # *** FIX: The 'reducedDim' accessor is in the SingleCellExperiment package ***
  corrected_pcs <- SingleCellExperiment::reducedDim(sce_corrected, "corrected")
  
  # Return as a data frame for easier handling
  return(as.data.frame(corrected_pcs))
}

#' Apply a PCA-based Batch Correction
#'
#' This method identifies principal components (PCs) that are significantly
#' associated with batch and regresses them out of the expression data.
#'
#' @param log_counts A matrix of log-transformed counts.
#' @param batch_info A character or factor vector with batch information.
#' @param sig_threshold A numeric p-value threshold to identify significant PCs.
#' @return A matrix of batch-corrected log-transformed counts.
#' @importFrom stats prcomp aov
#' @importFrom limma removeBatchEffect
run_pca_correction <- function(log_counts, batch_info, sig_threshold = 0.01) {
  # Perform PCA on the data
  pca <- stats::prcomp(t(log_counts), scale. = TRUE)
  
  # Identify PCs significantly associated with batch
  significant_pcs <- c()
  for (i in 1:ncol(pca$x)) {
    # Fit an ANOVA model: PC_i ~ batch
    fit <- stats::aov(pca$x[, i] ~ as.factor(batch_info))
    p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    # If the p-value is below the threshold, consider this PC batch-related
    if (!is.na(p_val) && p_val < sig_threshold) {
      significant_pcs <- c(significant_pcs, i)
    }
  }
  
  if (length(significant_pcs) == 0) {
    warning("No principal components were found to be significantly associated with batch. Returning original data.")
    return(log_counts)
  }
  
  cat(paste0("Found ", length(significant_pcs), " PCs significantly associated with batch: ", paste(significant_pcs, collapse=", "), "\n"))
  
  # Extract the values of the significant PCs
  batch_pcs <- pca$x[, significant_pcs, drop = FALSE]
  
  # Regress out the effect of these PCs from the original log-counts matrix.
  # We use limma's function with the 'covariates' argument for this.
  corrected_log_counts <- limma::removeBatchEffect(log_counts, covariates = batch_pcs)
  
  return(corrected_log_counts)
}

#' Calculate PERMANOVA R-squared for Batch Effect
#' @param log_counts A matrix of log-transformed counts. Optional if `pca_data` is provided.
#' @param pca_data A data frame of principal components (samples in rows, PCs in columns).
#' @return The R-squared value from the adonis test.
#' @importFrom vegan vegdist adonis2
calculate_permanova_r2 <- function(log_counts = NULL, batch_info, pca_data = NULL) {
  if (is.null(pca_data)) {
    dist_matrix <- vegan::vegdist(t(log_counts), method = "euclidean")
  } else {
    dist_matrix <- stats::dist(pca_data, method = "euclidean")
  }
  df <- data.frame(batch = as.factor(batch_info))
  permanova_res <- vegan::adonis2(dist_matrix ~ batch, data = df)
  return(permanova_res$R2[1])
}

#' Calculate Average Silhouette Width for Batches
#' @param log_counts A matrix of log-transformed counts. Optional if `pca_data` is provided.
#' @param pca_data A data frame of principal components (samples in rows, PCs in columns).
#' @return The average silhouette width.
#' @importFrom cluster silhouette
#' @importFrom stats dist prcomp
calculate_silhouette_width <- function(log_counts = NULL, batch_info, n_pcs = 5, pca_data = NULL) {
  if (is.null(pca_data)) {
    pca <- stats::prcomp(t(log_counts), scale. = TRUE)
    pca_data_subset <- pca$x[, 1:n_pcs]
  } else {
    # Ensure we don't try to select more PCs than are available
    max_pcs <- min(n_pcs, ncol(pca_data))
    pca_data_subset <- pca_data[, 1:max_pcs]
  }
  
  dist_matrix <- stats::dist(pca_data_subset)
  sil <- cluster::silhouette(x = as.numeric(as.factor(batch_info)), dist = dist_matrix)
  return(mean(sil[, "sil_width"]))
}