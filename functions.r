R_REMOTES_STANDALONE=TRUE

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
    empirical_means,
    empirical_libsizes,
    dropout_mid        = 1,
    donors             = paste0("Donor", 1:5),
    dispersion         = 0.2
) {
  # 1) Calibrate gene means & size factors from empirical distributions
  gene_means   <- sample_empirical(empirical_means, n_genes)
  size_factors <- sample_empirical(empirical_libsizes / mean(empirical_libsizes), n_cells)
  
  gene_names <- paste0("gene_", seq_len(n_genes))
  cell_names <- paste0("cell_", seq_len(n_cells))
  
  # --- FIX: Re-introduce the definition of 'lam' ---
  # 'lam' is the base expression rate (gene_mean * size_factor)
  lam <- outer(gene_means, size_factors, "*")
  
  # 2) Introduce donor-level effects
  n_donors <- length(donors)
  cells_per_donor <- ceiling(n_cells / n_donors)
  donor_assign <- rep(donors, each = cells_per_donor)[seq_len(n_cells)]
  donor_effects <- matrix(rnorm(n_genes * n_donors, mean = 1, sd = 0.2), nrow = n_genes)
  colnames(donor_effects) <- donors
  
  # 3) Simulate counts with NB and dropout (Vectorized)
  # Create the full mu matrix by multiplying base rates by donor effects
  mu_matrix <- lam * donor_effects[, donor_assign]
  
  # Simulate all counts at once
  raw_counts_matrix <- matrix(
    rnbinom(n_genes * n_cells, mu = mu_matrix, size = 1 / dispersion),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Calculate all dropout probabilities at once
  p_drop_matrix <- 1 / (1 + exp((log(mu_matrix + 1) - log(dropout_mid)) * 1.5))
  
  # Simulate all dropout events at once
  kept_matrix <- matrix(
    rbinom(n_genes * n_cells, 1, 1 - p_drop_matrix),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Apply dropout
  counts <- raw_counts_matrix * kept_matrix
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- gene_names; colnames(counts) <- cell_names
  
  # 4) Add ambient RNA contamination (~5% of reads)
  ambient_prof <- gene_means / sum(gene_means)
  lib_sizes <- Matrix::colSums(counts)
  # Prevent errors if some cells have zero counts
  valid_libs <- lib_sizes > 0
  if(any(valid_libs)){
    ambient_counts <- matrix(0, nrow = n_genes, ncol = n_cells)
    # Ensure rmultinom gets a single size value or a vector of the correct length
    ambient_sizes <- round(0.05 * lib_sizes[valid_libs])
    if(length(ambient_sizes) > 0) {
      ambient_counts[, valid_libs] <- stats::rmultinom(
        n = length(ambient_sizes),
        size = ambient_sizes,
        prob = ambient_prof
      )
    }
    counts <- counts + Matrix::Matrix(ambient_counts, sparse = TRUE)
  }
  
  
  # 5) Compute TPM from counts
  lib_sizes2 <- Matrix::colSums(counts)
  # Avoid division by zero for cells that might have zero counts after all steps
  lib_sizes2[lib_sizes2 == 0] <- 1
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
  method_order <- c("Pre-Correction", "Limma", "ComBat", "ComBat-Seq", "RUVs", "RUVg", "PCA Correction", "fastMNN")
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
  # ComBat-Seq requires a standard dense matrix.
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

#' Apply RUVs from RUVSeq with Upper‐Quartile Normalization
#'
#' Uses between‐lane UQ normalization then RUVs to remove unwanted variation using sample replicates.
#'
#' @param raw_counts Integer matrix (genes × samples).
#' @param control_genes Character vector of negative‐control gene names.
#' @param scIdx Integer matrix of replicate‐set sample indices (rows = sets, cols = max size; pad with -1).
#' @param k Numeric, number of unwanted factors to estimate.
#' @return A list with:
#'   - `corrected`: numeric matrix of RUV‐normalized counts
#'   - `k`: the k used
#' @importFrom EDASeq newSeqExpressionSet betweenLaneNormalization
#' @importFrom RUVSeq RUVs
#' @importFrom Biobase assayDataElement
run_ruvs <- function(raw_counts, control_genes, scIdx, k = 1) {
  # 1) Round & coerce
  counts_mat <- round(as.matrix(raw_counts))
  # 2) Build SeqExpressionSet
  seqset     <- EDASeq::newSeqExpressionSet(counts = counts_mat)
  # 3) UQ normalization
  seqset_uq  <- EDASeq::betweenLaneNormalization(seqset, which = "upper")
  # 4) Subset controls
  valid_ctrl <- intersect(control_genes, rownames(counts_mat))
  if (!length(valid_ctrl)) stop("No control genes found.")
  # 5) Run RUVs
  out        <- RUVSeq::RUVs(
    x     = seqset_uq,
    cIdx  = valid_ctrl,
    k     = as.numeric(k),
    scIdx = scIdx
  )
  # 6) Extract the 'normalizedCounts' assay
  corrected <- Biobase::assayDataElement(out, "normalizedCounts")
  list(corrected = corrected, k = k)
}

#' Apply RUVg from RUVSeq with Upper‐Quartile Normalization
#'
#' Uses between‐lane UQ normalization then RUVg to remove unwanted variation
#' using control genes. It does not use replicate information.
#'
#' @param raw_counts Integer matrix (genes × samples).
#' @param control_genes Character vector of negative‐control gene names.
#' @param k Numeric, number of unwanted factors to estimate.
#' @return A numeric matrix of RUVg‐normalized counts.
#' @importFrom EDASeq newSeqExpressionSet betweenLaneNormalization
#' @importFrom RUVSeq RUVg
#' @importFrom Biobase assayDataElement pData 'pData<-' AnnotatedDataFrame
run_ruvg <- function(raw_counts, control_genes, k = 1) {
  # 1) Round & coerce to a standard matrix
  counts_mat <- round(as.matrix(raw_counts))
  
  # 2) Build SeqExpressionSet. RUVg needs pData to exist.
  pheno_data <- data.frame(row.names = colnames(counts_mat),
                           dummy_var = rep(1, ncol(counts_mat)))
  seqset     <- EDASeq::newSeqExpressionSet(counts = counts_mat,
                                            phenoData = as(pheno_data, "AnnotatedDataFrame"))
  
  # 3) UQ normalization
  seqset_uq  <- EDASeq::betweenLaneNormalization(seqset, which = "upper")
  
  # 4) Ensure control genes are present in the data
  valid_ctrl <- intersect(control_genes, rownames(counts_mat))
  if (length(valid_ctrl) == 0) {
    stop("No control genes provided were found in the dataset for RUVg.")
  }
  
  # 5) Run RUVg using control genes
  out <- RUVSeq::RUVg(
    x    = seqset_uq,
    cIdx = valid_ctrl,
    k    = as.numeric(k)
  )
  
  # 6) Extract the 'normalizedCounts' assay
  corrected <- Biobase::assayDataElement(out, "normalizedCounts")
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
  
  # The 'reducedDim' accessor is in the SingleCellExperiment package
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

#' Run kBET (k-Nearest Neighbour Batch Effect Test).
#' @param log_counts A matrix of log-transformed data, used if pca_data is NULL.
#' @param batch_info A character or factor vector with batch information.
#' @param n_pcs Number of principal components to use.
#' @param pca_data A pre-computed matrix of PCs (samples x PCs).
#' @return The observed kBET rejection rate (0 to 1). Lower is better.
#' @importFrom kBET kBET
run_kbet <- function(log_counts = NULL, batch_info, n_pcs = 15, pca_data = NULL) {
  if (is.null(pca_data)) {
    if (ncol(log_counts) <= n_pcs) n_pcs <- ncol(log_counts) - 1
    if (n_pcs < 2 || nrow(log_counts) < n_pcs) {
      warning("Not enough samples/features for kBET. Returning NA.")
      return(NA)
    }
    pca <- stats::prcomp(t(log_counts), scale. = TRUE)
    data_for_kbet <- pca$x[, 1:n_pcs, drop = FALSE]
  } else {
    max_pcs <- min(n_pcs, ncol(pca_data))
    if (max_pcs < 2) {
      warning("Not enough PCs for kBET. Returning NA.")
      return(NA)
    }
    data_for_kbet <- pca_data[, 1:max_pcs, drop = FALSE]
  }
  min_batch_size <- min(table(batch_info))
  if (min_batch_size <= 10) {
    warning(paste("kBET cannot run: smallest batch is", min_batch_size, "(<=10). Returning NA."))
    return(NA)
  }
  k0 <- max(10, floor(min_batch_size * 0.25))
  if (k0 >= min_batch_size) k0 <- min_batch_size - 1
  
  kbet_results <- tryCatch({
    kBET::kBET(df = data_for_kbet, batch = batch_info, k0 = k0, plot = FALSE, verbose = FALSE)
  }, error = function(e) {
    warning(paste("kBET failed with k0 =", k0, ". Error:", e$message, ". Returning NA."))
    return(NULL)
  })
  
  if (is.null(kbet_results)) return(NA)
  
  return(kbet_results$summary$kBET.observed)
}



# Function to map a bibentry object to RIS lines
bibentry_to_ris <- function(b) {
  # helper to format authors
  format_authors <- function(auth) {
    names <- paste(auth$family, auth$given, sep = ", ")
    unlist(lapply(names, function(x) paste0("AU  - ", x)))
  }
  
  # ensure we have a single, character type string
  type_key <- as.character(b$bibtype)[1]
  
  # determine RIS type from that single key
  ris_type <- switch(type_key,
                     Article   = "JOUR",
                     Book      = "BOOK",
                     Manual    = "MANUAL",
                     Software  = "COMP",
                     Online    = "ELEC",
                     "GEN"     # generic fallback
  )
  
  out <- c(sprintf("TY  - %s", ris_type))
  
  if (!is.null(b$title))
    out <- c(out, sprintf("TI  - %s", b$title))
  if (!is.null(b$author))
    out <- c(out, format_authors(b$author))
  if (!is.null(b$year))
    out <- c(out, sprintf("PY  - %s", b$year))
  if (!is.null(b$journal))
    out <- c(out, sprintf("JO  - %s", b$journal))
  if (!is.null(b$publisher))
    out <- c(out, sprintf("PB  - %s", b$publisher))
  if (!is.null(b$volume))
    out <- c(out, sprintf("VL  - %s", b$volume))
  if (!is.null(b$issue))
    out <- c(out, sprintf("IS  - %s", b$issue))
  if (!is.null(b$pages))
    out <- c(out, sprintf("SP  - %s", b$pages))
  if (!is.null(b$url))
    out <- c(out, sprintf("UR  - %s", b$url))
  if (!is.null(b$doi))
    out <- c(out, sprintf("DO  - %s", b$doi))
  if (!is.null(b$note))
    out <- c(out, sprintf("N1  - %s", b$note))
  
  c(out, "ER  - ")
}

# Main function to export citations
export_citations_to_ris <- function(pkgs = NULL, file = "citations.ris") {
  if (is.null(pkgs)) {
    cits <- citation()
  } else {
    cits <- unlist(lapply(pkgs, citation), recursive = FALSE)
  }
  
  ris_lines <- unlist(lapply(cits, bibentry_to_ris))
  writeLines(ris_lines, con = file)
  message(sprintf("Written %d entries to '%s'", length(cits), file))
}

