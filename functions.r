#' Install and Load Required R Packages
#'
#' Handles installation from CRAN, Bioconductor, and GitHub.
#' @param packages_to_install A character vector of package names.
install_and_load_packages <- function(packages_to_install) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  
  for (pkg in packages_to_install) {
    is_github <- grepl("/", pkg)
    pkg_name <- if(is_github) basename(pkg) else pkg
    
    if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
      if (is_github) {
        remotes::install_github(pkg, force = TRUE) # Force install to get latest
      } else {
        BiocManager::install(pkg, update = FALSE)
      }
    }
    library(pkg_name, character.only = TRUE)
  }
}

#' Sample from an empirical distribution
sample_empirical <- function(x, n) sample(x, size = n, replace = TRUE)

#' Create a base single-cell dataset
create_base_sc_dataset <- function(n_genes=1000, n_cells=300, empirical_means,
                                   empirical_libsizes, dropout_mid=1,
                                   donors=paste0("Donor", 1:5), dispersion=0.2) {
  gene_means <- sample_empirical(empirical_means, n_genes)
  size_factors <- sample_empirical(empirical_libsizes / mean(empirical_libsizes), n_cells)
  gene_names <- paste0("gene_", seq_len(n_genes))
  cell_names <- paste0("cell_", seq_len(n_cells))
  lam <- outer(gene_means, size_factors, "*")
  n_donors <- length(donors)
  cells_per_donor <- ceiling(n_cells / n_donors)
  donor_assign <- rep(donors, each = cells_per_donor)[seq_len(n_cells)]
  donor_effects <- matrix(rnorm(n_genes * n_donors, mean = 1, sd = 0.2), nrow = n_genes)
  colnames(donor_effects) <- donors
  mu_matrix <- lam * donor_effects[, donor_assign]
  raw_counts <- matrix(rnbinom(n_genes * n_cells, mu = mu_matrix, size = 1 / dispersion), nrow = n_genes)
  p_drop <- 1 / (1 + exp((log(mu_matrix + 1) - log(dropout_mid)) * 1.5))
  kept <- matrix(rbinom(n_genes * n_cells, 1, 1 - p_drop), nrow = n_genes)
  counts <- raw_counts * kept
  
  # Add ambient RNA contamination
  ambient_prof <- gene_means / sum(gene_means)
  lib_sizes <- colSums(counts)
  valid_libs <- lib_sizes > 0
  if (any(valid_libs)) {
    ambient_counts <- matrix(0, nrow = n_genes, ncol = n_cells)
    ambient_sizes <- round(0.05 * lib_sizes[valid_libs])
    
    if(length(ambient_sizes) > 0) {
      generated_ambient_counts <- sapply(ambient_sizes, function(s) {
        stats::rmultinom(n = 1, size = s, prob = ambient_prof)
      })
      ambient_counts[, valid_libs] <- generated_ambient_counts
    }
    counts <- counts + ambient_counts
  }
  
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- gene_names; colnames(counts) <- cell_names
  lib_sizes2 <- Matrix::colSums(counts)
  lib_sizes2[lib_sizes2 == 0] <- 1
  tpm <- Matrix::t(1e6 * Matrix::t(counts) / lib_sizes2)
  annotation <- data.frame(
    ID = cell_names, donor = donor_assign,
    cell_type = sample(c("T cells", "B cells", "Monocytes"), size = n_cells, replace = TRUE)
  )
  SimBu::dataset(annotation=annotation, count_matrix=counts, tpm_matrix=tpm, name="sc_dataset")
}

#' Filter genes by specificity
filter_genes_by_specificity <- function(simbu_ds, top_n=200) {
  counts_mat <- SummarizedExperiment::assay(simbu_ds, "counts")
  cell_types <- SummarizedExperiment::colData(simbu_ds)$cell_type
  types <- unique(cell_types)
  mean_mat <- sapply(types, function(t) rowMeans(counts_mat[, cell_types == t, drop=FALSE]))
  
  spec_score <- matrixStats::rowMaxs(mean_mat, na.rm = TRUE) / matrixStats::rowSums2(mean_mat, na.rm = TRUE)
  
  top_genes <- names(sort(spec_score, decreasing = TRUE))[1:min(top_n, length(spec_score))]
  simbu_ds[top_genes, ]
}

#' Generate bulk RNA-seq simulations
generate_simulations <- function(simbu_dataset, n_sims=3, n_samples_per_sim=30, n_cells_per_sample=100) {
  types <- unique(colData(simbu_dataset)$cell_type)
  lapply(1:n_sims, function(i) {
    alpha_p <- pmax(rnorm(1, mean=1, sd=0.5), 0.1)
    fracs <- gtools::rdirichlet(n_samples_per_sim, rep(alpha_p, length(types)))
    colnames(fracs) <- types
    
    SimBu::simulate_bulk(
      data                 = simbu_dataset,
      scenario             = "custom",
      custom_scenario_data = fracs,
      scaling_factor       = "NONE",
      ncells               = n_cells_per_sample,
      nsamples             = n_samples_per_sim
    )
  })
}

#' Create a Non-Confounded Dataset with Batch Effect
create_confound_free_dataset <- function(sim_list, n_genes_affected, effect_multiplier) {
  for (i in seq_along(sim_list)) {
    colData(sim_list[[i]]$bulk)$biological_group <- paste0("Sim", i)
  }
  merged_se <- SimBu::merge_simulations(sim_list)$bulk
  colnames(merged_se) <- make.names(colnames(merged_se), unique = TRUE)
  n_samples <- ncol(merged_se)
  n_batches <- length(sim_list)
  samples_per_batch <- n_samples / n_batches
  set.seed(456)
  batch_info <- sample(rep(paste0("Batch", 1:n_batches), each = samples_per_batch))
  colData(merged_se)$batch <- batch_info
  biological_vars <- merged_se$biological_group
  counts_matrix <- assay(merged_se, "bulk_counts")
  set.seed(789)
  affected_genes <- sample(rownames(counts_matrix), n_genes_affected)
  b2_indices <- which(batch_info == "Batch2")
  b3_indices <- which(batch_info == "Batch3")
  counts_matrix[affected_genes, b2_indices] <- round(counts_matrix[affected_genes, b2_indices] * effect_multiplier)
  counts_matrix[affected_genes, b3_indices] <- round(counts_matrix[affected_genes, b3_indices] / effect_multiplier)
  assay(merged_se, "bulk_counts") <- counts_matrix
  list(se_with_batch=merged_se, batch_info=batch_info, biological_vars=biological_vars, affected_genes=affected_genes)
}

# --- METRIC AND CORRECTION FUNCTIONS ---

#' Perform PCA and plot
perform_pca_and_plot <- function(log_counts=NULL, batch_info, plot_title, pca_data=NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    log_counts_filtered <- log_counts[gene_vars > 1e-8, ]
    if(nrow(log_counts_filtered) < 2) {
      warning(paste("PCA for", plot_title, "skipped: <2 genes with variance.")); return(ggplot())
    }
    pca_res <- stats::prcomp(t(log_counts_filtered), scale. = TRUE)
    pca_data <- pca_res$x
  }
  plot_df <- data.frame(PC1=pca_data[,1], PC2=pca_data[,2], batch=batch_info)
  ggplot2::ggplot(plot_df, ggplot2::aes(x=PC1, y=PC2, color=batch)) +
    ggplot2::geom_point(size=2, alpha=0.8) + ggplot2::labs(title=plot_title) +
    ggplot2::theme_bw() + ggplot2::scale_color_brewer(palette="Set1")
}

#' Run ComBat
run_combat <- function(log_counts, batch_info, mod, ref_batch=NULL) {
  sva::ComBat(dat=as.matrix(log_counts), batch=as.factor(batch_info), mod=mod, ref.batch=ref_batch)
}

#' Run ComBat-Seq
run_combat_seq <- function(raw_counts, batch_info, biological_vars=NULL) {
  sva::ComBat_seq(counts=as.matrix(raw_counts), batch=as.factor(batch_info), group=biological_vars)
}

#' Run RUVs
run_ruvs <- function(seqset_uq, control_genes, scIdx, k=1) {
  valid_ctrl <- intersect(control_genes, rownames(seqset_uq))
  out <- RUVSeq::RUVs(x=seqset_uq, cIdx=valid_ctrl, k=k, scIdx=scIdx)
  list(corrected=Biobase::assayDataElement(out, "normalizedCounts"))
}

#' Run RUVg
run_ruvg <- function(seqset_uq, control_genes, k=1) {
  valid_ctrl <- intersect(control_genes, rownames(seqset_uq))
  out <- RUVSeq::RUVg(x=seqset_uq, cIdx=valid_ctrl, k=k)
  Biobase::assayDataElement(out, "normalizedCounts")
}

#' Run fastMNN
run_fastmnn <- function(log_counts, batch_info) {
  as.data.frame(SingleCellExperiment::reducedDim(batchelor::fastMNN(log_counts, batch=as.factor(batch_info))))
}

#' Run PCA-based correction
run_pca_correction <- function(log_counts, batch_info, sig_threshold=0.01) {
  pca <- stats::prcomp(t(log_counts), scale.=TRUE)
  p_vals <- sapply(1:ncol(pca$x), function(i) summary(aov(pca$x[,i]~as.factor(batch_info)))[[1]][["Pr(>F)"]][1])
  sig_pcs <- which(p_vals < sig_threshold)
  if (length(sig_pcs) == 0) { warning("No significant PCs for batch. Returning original data."); return(log_counts) }
  limma::removeBatchEffect(log_counts, covariates=pca$x[, sig_pcs, drop=FALSE])
}

#' Run SVA
run_sva <- function(log_counts, mod, mod0) {
  sv_obj <- sva::sva(dat = as.matrix(log_counts), mod = mod, mod0 = mod0)
  if (sv_obj$n.sv > 0) {
    return(limma::removeBatchEffect(log_counts, covariates = sv_obj$sv))
  } else {
    warning("SVA found 0 surrogate variables. Returning original data.")
    return(log_counts)
  }
}

#' Run SVA-Seq
run_svaseq <- function(raw_counts, mod, mod0) {
  sv_obj <- sva::svaseq(dat = as.matrix(raw_counts), mod = mod, mod0 = mod0)
  if (sv_obj$n.sv > 0) {
    return(limma::removeBatchEffect(edgeR::cpm(raw_counts, log = TRUE, prior.count = 1), covariates = sv_obj$sv))
  } else {
    warning("SVA-Seq found 0 surrogate variables. Returning original log-counts.")
    return(edgeR::cpm(raw_counts, log = TRUE, prior.count = 1))
  }
}

# --- METRIC CALCULATION FUNCTIONS ---

#' Calculate PERMANOVA R-squared
calculate_permanova_r2 <- function(log_counts=NULL, variable_info, pca_data=NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    if(sum(gene_vars > 1e-8) < 2) return(NA)
    pca_data <- stats::prcomp(t(log_counts[gene_vars > 1e-8, ]), scale.=TRUE)$x
  }
  res <- tryCatch(vegan::adonis2(dist(pca_data) ~ as.factor(variable_info)), error=function(e) NULL)
  if (is.null(res)) NA else res$R2[1]
}

#' Calculate Silhouette Width
calculate_silhouette_width <- function(log_counts=NULL, variable_info, n_pcs=5, pca_data=NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    if(sum(gene_vars > 1e-8) < 2) return(NA)
    pca_data <- stats::prcomp(t(log_counts[gene_vars > 1e-8, ]), scale.=TRUE)$x
  }
  if (length(unique(variable_info)) < 2) return(NA)
  max_pcs <- min(n_pcs, ncol(pca_data))
  sil <- cluster::silhouette(x=as.numeric(as.factor(variable_info)), dist=dist(pca_data[, 1:max_pcs]))
  mean(sil[, "sil_width"])
}

#' Calculate PC Regression R-squared
calculate_pc_regression <- function(log_counts=NULL, variable_info, n_pcs=10, pca_data=NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    if(sum(gene_vars > 1e-8) < 2) return(NA)
    pca_data <- stats::prcomp(t(log_counts[gene_vars > 1e-8, ]), scale.=TRUE)$x
  }
  max_pcs <- min(n_pcs, ncol(pca_data))
  r_sqs <- sapply(1:max_pcs, function(i) summary(lm(pca_data[,i]~as.factor(variable_info)))$r.squared)
  mean(r_sqs)
}

#' Calculate LISI
run_lisi <- function(log_counts=NULL, meta_data, var_to_check, n_pcs=15, pca_data=NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    if(sum(gene_vars > 1e-8) < 2) return(NA)
    pca_data <- stats::prcomp(t(log_counts[gene_vars > 1e-8, ]), scale.=TRUE)$x
  }
  max_pcs <- min(n_pcs, ncol(pca_data))
  res <- tryCatch(lisi::compute_lisi(pca_data[,1:max_pcs, drop=FALSE], meta_data, var_to_check), error=function(e) NULL)
  if(is.null(res)) NA else mean(res[[var_to_check]], na.rm=TRUE)
}

#' Calculate Ground-Truth Gene R-squared
calculate_gene_metric <- function(log_counts, batch_info, affected_genes) {
  genes_to_test <- intersect(affected_genes, rownames(log_counts))
  if (length(genes_to_test) == 0) return(NA)
  r_sqs <- sapply(genes_to_test, function(g) summary(lm(log_counts[g,]~as.factor(batch_info)))$r.squared)
  mean(r_sqs, na.rm=TRUE)
}

#' Run kBET (Robust Version)
run_kbet <- function(log_counts = NULL, batch_info, n_pcs = 15, pca_data = NULL) {
  if (is.null(pca_data)) {
    gene_vars <- matrixStats::rowVars(as.matrix(log_counts))
    if(sum(gene_vars > 1e-8) < n_pcs + 1) {
      warning("kBET skipped: Not enough variable genes for PCA.")
      return(NA)
    }
    pca_res <- stats::prcomp(t(log_counts[gene_vars > 1e-8, ]), scale. = TRUE)
    pca_data <- pca_res$x
  }
  
  max_pcs <- min(n_pcs, ncol(pca_data))
  if (max_pcs < 2) {
    warning("kBET skipped: Fewer than 2 PCs available.")
    return(NA)
  }
  data_for_kbet <- pca_data[, 1:max_pcs, drop = FALSE]
  
  min_batch_size <- min(table(batch_info))
  k0 <- floor(min_batch_size * 0.25)
  if (k0 < 10) {
    warning(paste("kBET skipped: Smallest batch size is", min_batch_size, "which is too small."))
    return(NA)
  }
  
  kbet_results <- tryCatch({
    kBET::kBET(df = data_for_kbet, batch = batch_info, k0 = k0, plot = FALSE, do.pca=FALSE, verbose=FALSE)
  }, error = function(e) {
    warning(paste("kBET failed with error:", e$message))
    return(NULL)
  })
  
  if (is.null(kbet_results)) return(NA)
  return(kbet_results$summary$kBET.observed)
}
