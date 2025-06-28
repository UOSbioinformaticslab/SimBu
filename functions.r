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

#' Sample from an empirical distribution
sample_empirical <- function(x, n) {
  sample(x, size = n, replace = TRUE)
}


#' Clinical-calibrated Base Single-Cell Dataset
#'
#' @param n_genes Number of genes to simulate.
#' @param n_cells Number of cells to simulate.
#' @param empirical_means Numeric vector of real gene mean expressions.
#' @param empirical_libsizes Numeric vector of real total UMIs per cell.
#' @param dropout_mid Midpoint for logistic dropout curve.
#' @param donor_ids Vector of donor identifiers for patient variability.
#' @return A SimBu_dataset with aligned count and TPM matrices and clinical realism.
create_base_sc_dataset <- function(
    n_genes            = 1000,
    n_cells            = 300,
    empirical_means,           # from real clinical data
    empirical_libsizes,        # from real clinical data
    dropout_mid       = 1,
    donors            = paste0("Donor", 1:5)
) {
  # 1) Calibrate gene means & size factors from empirical distributions
  gene_means   <- sample_empirical(empirical_means, n_genes)
  size_factors <- sample_empirical(empirical_libsizes / mean(empirical_libsizes), n_cells)
  
  gene_names <- paste0("gene_", seq_len(n_genes))
  cell_names <- paste0("cell_", seq_len(n_cells))
  
  # 2) Introduce donor-level effects
  n_donors <- length(donors)
  cells_per_donor <- ceiling(n_cells / n_donors)
  donor_assign <- rep(donors, each = cells_per_donor)[seq_len(n_cells)]
  # multiplicative effects per donor
  donor_effects <- matrix(rnorm(n_genes * n_donors, mean = 1, sd = 0.2), nrow = n_genes)
  colnames(donor_effects) <- donors
  
  # 3) Simulate counts with NB and dropout
  lam <- outer(gene_means, size_factors)
  counts <- matrix(0, nrow = n_genes, ncol = n_cells)
  for (i in seq_len(n_cells)) {
    mu_i <- lam[, i] * donor_effects[, donor_assign[i]]
    raw_counts <- rnbinom(n_genes, mu = mu_i, size = 0.5)
    # dropout probability logistic function
    p_drop <- 1 / (1 + exp((log(mu_i + 1) - log(dropout_mid)) * 1.5))
    kept <- rbinom(n_genes, 1, 1 - p_drop)
    counts[, i] <- raw_counts * kept
  }
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- gene_names; colnames(counts) <- cell_names
  
  # 4) Add ambient RNA contamination (~5% of reads)
  ambient_prof <- gene_means / sum(gene_means)
  lib_sizes <- Matrix::colSums(counts)
  ambient_counts <- matrix(
    stats::rmultinom(n_cells, size = round(0.05 * lib_sizes), prob = ambient_prof),
    nrow = n_genes
  )
  counts <- counts + Matrix::Matrix(ambient_counts, sparse = TRUE)
  
  # 5) Compute TPM from counts
  lib_sizes2 <- Matrix::colSums(counts)
  raw_tpm <- Matrix::t(1e6 * Matrix::t(counts) / lib_sizes2)
  tpm <- Matrix::Matrix(raw_tpm, sparse = TRUE)
  rownames(tpm) <- gene_names; colnames(tpm) <- cell_names
  
  # 6) Annotate cells
  annotation <- data.frame(
    ID         = cell_names,
    donor      = donor_assign,
    cell_type  = sample(
      c("T cells CD4", "T cells CD8", "Macrophages", "NK cells", "B cells", "Monocytes"),
      size = n_cells, replace = TRUE,
      prob = c(0.15, 0.15, 0.30, 0.05, 0.25, 0.10)
    ),
    stringsAsFactors = FALSE
  )
  
  SimBu::dataset(
    annotation   = annotation,
    count_matrix = counts,
    tpm_matrix   = tpm,
    name         = "clinical_sc_dataset"
  )
}

#' Filter to the top N cell-type–specific genes
filter_genes_by_specificity <- function(simbu_ds, top_n = 200) {
  counts_mat <- SummarizedExperiment::assay(simbu_ds, "counts")
  cell_types <- SummarizedExperiment::colData(simbu_ds)$cell_type
  types      <- unique(cell_types)
  
  mean_mat <- sapply(types, function(t) rowMeans(counts_mat[, cell_types == t, drop=FALSE]))
  rownames(mean_mat) <- rownames(counts_mat)
  spec_score <- apply(mean_mat, 1, function(x) max(x) / sum(x))
  top_genes  <- names(sort(spec_score, decreasing = TRUE))[seq_len(min(top_n, length(spec_score)))]
  
  return(simbu_ds[top_genes, ])
}

#' Generate Bulk RNA-seq Simulations with Patient-level Heterogeneity
generate_simulations <- function(
    simbu_dataset,
    n_sims            = 3,
    n_samples_per_sim = 30,
    n_cells_per_sample= 100
) {
  types <- unique(colData(simbu_dataset)$cell_type)
  sim_list <- vector("list", n_sims)
  
  for (i in seq_len(n_sims)) {
    # sample patient-level Dirichlet concentration
    alpha_p <- pmax(rnorm(1, mean = 1, sd = 0.5), 0.1)
    fracs   <- gtools::rdirichlet(n_samples_per_sim, rep(alpha_p, length(types)))
    colnames(fracs) <- types
    
    sim_list[[i]] <- SimBu::simulate_bulk(
      data                 = simbu_dataset,
      scenario             = "custom",
      custom_scenario_data = fracs,
      scaling_factor       = "NONE",
      ncells               = n_cells_per_sample,
      nsamples             = n_samples_per_sim,
      run_parallel         = FALSE
    )
  }
  return(sim_list)
}

#' Merge Simulations and Introduce a Batch Effect
introduce_batch_effect <- function(sim_list,
                                   n_genes_affected  = 50,
                                   effect_multiplier = 7.0) {
  merged_list <- SimBu::merge_simulations(sim_list)
  se_object   <- merged_list$bulk
  colnames(se_object) <- make.names(colnames(se_object), unique = TRUE)
  
  counts_matrix <- assay(se_object, "bulk_counts")
  n_per_batch   <- ncol(sim_list[[1]]$bulk)
  batch_info    <- rep(paste0("Batch", seq_along(sim_list)), each = n_per_batch)
  
  set.seed(456)
  affected_genes <- sample(rownames(counts_matrix), n_genes_affected)
  
  b2 <- which(batch_info == "Batch2"); b3 <- which(batch_info == "Batch3")
  counts_matrix[affected_genes, b2] <- round(counts_matrix[affected_genes, b2] * effect_multiplier)
  counts_matrix[affected_genes, b3] <- round(counts_matrix[affected_genes, b3] / effect_multiplier)
  
  assay(se_object, "bulk_counts") <- counts_matrix
  
  list(
    se_with_batch  = se_object,
    batch_info     = batch_info,
    affected_genes = affected_genes
  )
}

#' Native rarefaction (multinomial)
rarefy_counts <- function(counts, depth) {
  p <- counts / sum(counts)
  as.numeric(rmultinom(1, size = depth, prob = p))
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
#' @param log_counts A matrix of log-transformed counts. Ignored if `pca_data` is provided.
#' @param batch_info A character or factor vector with batch information.
#' @param plot_title A string for the plot title.
#' @param pca_data A matrix or data.frame of precomputed PCs (samples × components). If provided, `log_counts` is ignored.
#' @return A `ggplot` object.
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw scale_color_brewer coord_fixed
perform_pca_and_plot <- function(log_counts = NULL,
                                 batch_info,
                                 plot_title,
                                 pca_data = NULL) {
  if (is.null(pca_data)) {
    # 1) compute PCA from log-counts
    pca_res <- stats::prcomp(t(log_counts), scale. = TRUE)
    # percent variance for axes
    pct <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
    plot_df <- data.frame(
      PC1   = pca_res$x[, 1],
      PC2   = pca_res$x[, 2],
      batch = batch_info,
      check.names = FALSE
    )
    x_lab <- paste0("PC1 (", pct[1], "%)")
    y_lab <- paste0("PC2 (", pct[2], "%)")
    
  } else {
    # 2) use the first two columns of the provided PC matrix/frame
    plot_df <- data.frame(
      PC1   = pca_data[, 1],
      PC2   = pca_data[, 2],
      batch = batch_info,
      check.names = FALSE
    )
    x_lab <- "PC1"
    y_lab <- "PC2"
  }
  
  # now the scatter
  ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = batch)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(title = plot_title, x = x_lab, y = y_lab) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::coord_fixed()
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
#' @param raw_counts A matrix of non-normalized, integer counts (genes in rows, samples in columns).
#' Can be a standard or sparse matrix.
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
#' @param log_counts A matrix of log-transformed counts.
#' @param batch_info A character or factor vector with batch information.
#' @param ref_batch Optional string specifying a reference batch. If provided, all other batches will be adjusted to match this reference to help avoid overcorrection.
#' @return A matrix of batch-corrected log-transformed counts.
#' @importFrom sva ComBat
#' @importFrom stats model.matrix
run_combat <- function(log_counts, batch_info, ref_batch = NULL) {
  # ComBat expects standard dense matrix.
  log_counts_matrix <- as.matrix(log_counts)
  batch_factor <- as.factor(batch_info)
  
  # Use an intercept-only model as there are no covariates to preserve.
  # This is the 'mod' argument in ComBat.
  mod_combat <- stats::model.matrix(~1, data = data.frame(batch = batch_factor))
  
  # Run ComBat on the log-transformed data
  corrected_log_counts <- sva::ComBat(
    dat = log_counts_matrix,
    batch = batch_factor,
    mod = mod_combat,
    ref.batch = ref_batch
  )
  
  return(corrected_log_counts)
}

#' Apply RUVg from RUVSeq for Batch Correction
#'
#' Applies the RUVg algorithm, which uses control genes to estimate and remove
#' unwanted variation.  It operates on raw, integer counts.
#'
#' @param raw_counts A matrix of non-normalized, integer counts.
#' @param control_genes A character vector of gene names to be used as controls.
#' @param k An integer, the number of factors of unwanted variation to remove.
#' @return A matrix of RUVg‐corrected counts.
#' @importFrom EDASeq newSeqExpressionSet
#' @importFrom RUVSeq RUVg
#' @importFrom Biobase assayDataElement
run_ruvg <- function(raw_counts, control_genes, k = 1) {
  # 1) round to integer
  counts_matrix <- round(as.matrix(raw_counts))
  
  # 2) pick only the control genes that exist
  valid_ctrl <- intersect(control_genes, rownames(counts_matrix))
  if (length(valid_ctrl) == 0) {
    stop("None of the provided control genes were found in the count matrix.")
  }
  
  # 3) build the SeqExpressionSet (counts go in the 'counts' slot)
  seqset <- EDASeq::newSeqExpressionSet(counts = counts_matrix)
  
  # 4) run RUVg – this returns a SeqExpressionSet with a 'normalizedCounts' assay
  seqset_ruv <- RUVSeq::RUVg(x    = seqset,
                             cIdx = valid_ctrl,
                             k    = k)
  
  # 5) extract the corrected counts matrix from the normalizedCounts slot
  corrected <- Biobase::assayDataElement(seqset_ruv, "normalizedCounts")
  
  return(corrected)
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

#' Calculate PERMANOVA R-squared for a given variable (batch or biological group).
#'
#' @param log_counts Matrix of log-counts. Required if `pca_data` is NULL.
#' @param variable_info A character or factor vector with the variable of interest (e.g., batch_info, group_info).
#' @param pca_data Data frame with PCs. If provided, `log_counts` is ignored.
#' @return The R-squared value from the adonis2 test.
#' @importFrom vegan vegdist adonis2
calculate_permanova_r2 <- function(log_counts = NULL, variable_info, pca_data = NULL) {
  if (is.null(pca_data)) {
    dist_matrix <- vegan::vegdist(t(log_counts), method = "euclidean")
  } else {
    dist_matrix <- stats::dist(pca_data, method = "euclidean")
  }
  df <- data.frame(variable = as.factor(variable_info))
  # PERMANOVA can fail if a group has only one member, so wrap in tryCatch
  permanova_res <- tryCatch({
    vegan::adonis2(dist_matrix ~ variable, data = df)
  }, error = function(e) {
    warning("PERMANOVA failed, likely due to a group having only one sample. Returning NA.")
    return(NULL)
  })
  
  if (is.null(permanova_res)) return(NA)
  
  return(permanova_res$R2[1])
}

#' Calculate Average Silhouette Width for a given variable (batch or biological group).
#'
#' @param log_counts Matrix of log-counts. Required if `pca_data` is NULL.
#' @param variable_info A character or factor vector with the variable of interest (e.g., batch_info, group_info).
#' @param n_pcs Number of principal components to use for the calculation.
#' @param pca_data Data frame with PCs. If provided, `log_counts` is ignored.
#' @return The average silhouette width. For batch, lower is better. For biology, higher is better.
#' @importFrom cluster silhouette
#' @importFrom stats dist
calculate_silhouette_width <- function(log_counts = NULL, variable_info, n_pcs = 5, pca_data = NULL) {
  if (is.null(pca_data)) {
    # Ensure there are enough samples to compute PCs
    if (ncol(log_counts) <= n_pcs) {
      n_pcs <- ncol(log_counts) - 1
    }
    if (n_pcs < 2) {
      warning("Not enough samples to compute multiple PCs for Silhouette Width. Returning NA.")
      return(NA)
    }
    pca <- stats::prcomp(t(log_counts), scale. = TRUE)
    pca_data_subset <- pca$x[, 1:n_pcs, drop = FALSE]
  } else {
    max_pcs <- min(n_pcs, ncol(pca_data))
    pca_data_subset <- pca_data[, 1:max_pcs, drop = FALSE]
  }
  
  dist_matrix <- stats::dist(pca_data_subset)
  # Ensure there's more than one unique group
  if (length(unique(variable_info)) < 2) {
    warning("Silhouette width cannot be calculated with only one group. Returning NA.")
    return(NA)
  }
  sil <- cluster::silhouette(x = as.numeric(as.factor(variable_info)), dist = dist_matrix)
  
  return(mean(sil[, "sil_width"]))
}

#' Run kBET (k-Nearest Neighbor Batch Effect Test).
#'
#' Calculates the rejection rate of a chi-squared test for local batch mixing.
#' A lower rejection rate indicates better batch integration.
#'
#' @param log_counts Matrix of log-counts. Required if `pca_data` is NULL.
#' @param batch_info A character or factor vector with batch information.
#' @param n_pcs Number of principal components to use for the calculation.
#' @param pca_data Data frame with PCs. If provided, `log_counts` is ignored.
#' @return The kBET rejection rate. Lower is better.
#' @importFrom kbet kbet
run_kbet <- function(log_counts = NULL, batch_info, n_pcs = 10, pca_data = NULL) {
  if (is.null(pca_data)) {
    # kBET works best on PCs, so we compute them if not provided
    if (ncol(log_counts) <= n_pcs) {
      n_pcs <- ncol(log_counts) - 1
    }
    if (n_pcs < 2) {
      warning("Not enough samples to compute multiple PCs for kBET. Returning NA.")
      return(NA)
    }
    pca <- prcomp(t(log_counts), scale. = TRUE)
    data_for_kbet <- pca$x[, 1:n_pcs]
  } else {
    max_pcs <- min(n_pcs, ncol(pca_data))
    data_for_kbet <- pca_data[, 1:max_pcs]
  }
  
  # kBET can fail if batches are too small, so we wrap it in a tryCatch
  kbet_results <- tryCatch({
    kbet::kbet(df = data_for_kbet, batch = batch_info, plot = FALSE)
  }, error = function(e) {
    warning("kBET failed, likely due to small batch size. Returning NA.")
    return(NULL)
  })
  
  if (is.null(kbet_results)) return(NA)
  
  # Return the observed rejection rate
  return(kbet_results$summary$kbet.observed)
}