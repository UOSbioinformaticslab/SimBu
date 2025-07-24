# --- 0. Setup ---
# source("D:/rstudio/SimBu/functions.r")

# Define high-level parameters
PACKAGES <- c(
  "SimBu", "limma",   "edgeR", "dplyr",  "ggplot2",
  "patchwork", "SummarizedExperiment", "Matrix",
  "sva",  "RUVSeq", "EDASeq", "vegan",
  "cluster", "tidyr", "knitr", "batchelor",
  "SingleCellExperiment", "gtools", "R.cache",
  "matrixStats", "lisi", "remotes", "welch-lab/kBET"
)


NSAMPLES_PER_BATCH <- 50
N_GENES            <- 1000
N_CELLS            <- 300
TOP_GENES          <- 200
# --- Parameter Adjustments for Realistic Batch Effects ---
# We are increasing the complexity and strength of the simulated batch effects
# to better mimic challenging real-world scenarios.

# 1. Increased number of batches to add complexity.
NUM_SIMS           <- 4      # Original: 4
# 2. Baseline gene means remain the same.
MEAN_SHAPE   <- 1.5
MEAN_SCALE   <- 0.5
# 3. Increased library size variability for more realistic technical noise.
SIZE_MEANLOG <- 0
SIZE_SDLOG   <- 0.5    # Original: 0.5
# 4. Increased dispersion to model higher biological/technical noise.
DISPERSION   <- 0.4    # Original: 0.4
# 5. Increased number of affected genes to make the batch effect more pervasive.
N_AFFECTED_GENES  <- 150    # Original: 50
# 6. Increased effect multiplier for a stronger, more challenging batch effect.
EFFECT_MULTIPLIER <- 3.5    # Original: 3.5


# --- 1. Environment Setup ---
install_and_load_packages(PACKAGES)
set.seed(123) # Set seed for reproducibility

# --- 2. Data Generation and Simulation ---

# <<< CACHING START for Data Generation >>>
cat("--- Checking cache for prepared dataset... ---\n")
data_params_key <- list(
  N_GENES=N_GENES, N_CELLS=N_CELLS, TOP_GENES=TOP_GENES, NUM_SIMS=NUM_SIMS,
  NSAMPLES_PER_BATCH=NSAMPLES_PER_BATCH, EFFECT_MULTIPLIER=EFFECT_MULTIPLIER,
  N_AFFECTED_GENES=N_AFFECTED_GENES, DISPERSION=DISPERSION, script_version="v1"
)
prepared_data <- R.cache::loadCache(key = data_params_key)

if (is.null(prepared_data)) {
  cat("--- No cache found. Generating and preparing dataset... ---\n")
  empirical_gene_means <- rgamma(N_GENES, shape = MEAN_SHAPE, scale = MEAN_SCALE)
  empirical_lib_sizes  <- rlnorm(N_CELLS, meanlog = SIZE_MEANLOG, sdlog = SIZE_SDLOG)
  base_ds <- create_base_sc_dataset(
    n_genes=N_GENES, n_cells=N_CELLS, empirical_means=empirical_gene_means,
    empirical_libsizes=empirical_lib_sizes, dispersion=DISPERSION
  )
  filtered_ds <- filter_genes_by_specificity(base_ds, top_n = TOP_GENES)
  sim_list <- generate_simulations(
    simbu_dataset=filtered_ds, n_sims=NUM_SIMS, n_samples_per_sim=NSAMPLES_PER_BATCH
  )
  prepared_data <- create_confound_free_dataset(
    sim_list, n_genes_affected=N_AFFECTED_GENES, effect_multiplier=EFFECT_MULTIPLIER
  )
  R.cache::saveCache(prepared_data, key = data_params_key)
  cat("--- Dataset generation complete. Saved to cache. ---\n")
} else {
  cat("--- Loaded prepared dataset from cache. ---\n")
}
# <<< CACHING END for Data Generation >>>


# --- 3. Prepare Non-Confounded Dataset ---
counts_with_effect <- assay(prepared_data$se_with_batch, "bulk_counts")
batch_info         <- prepared_data$batch_info
biological_vars    <- prepared_data$biological_vars
affected_genes     <- prepared_data$affected_genes
meta_data <- data.frame(
  batch = batch_info, biology = biological_vars, row.names = colnames(counts_with_effect)
)

# --- 4. Pre-Correction Data & Pre-computation ---
log_counts_before  <- edgeR::cpm(as.matrix(counts_with_effect), log=TRUE, prior.count=1)
mod_bio <- model.matrix(~as.factor(biological_vars))
mod_null <- model.matrix(~1, data = as.data.frame(biological_vars))
mod_combat <- model.matrix(~1, data=data.frame(batch=batch_info))


# --- 5. Run Correction Methods ---
cat("--- Running Batch Correction Methods ---\n")

log_limma <- limma::removeBatchEffect(log_counts_before, batch = batch_info)
log_combat <- run_combat(log_counts_before, batch_info, mod=mod_combat, ref_batch = "Batch1")
counts_seq <- run_combat_seq(counts_with_effect, batch_info, biological_vars = biological_vars)
log_seq    <- edgeR::cpm(counts_seq, log = TRUE, prior.count = 1)

cat("Tuning RUVs and RUVg...\n")
logcpm     <- edgeR::cpm(counts_with_effect, log = TRUE, prior.count = 1)
gene_var   <- matrixStats::rowVars(logcpm)
ctrl_cands <- names(sort(gene_var))[1:50]
ctrl_genes <- setdiff(ctrl_cands, affected_genes)

# --- ROBUSTNESS FIX 2: Redefine control gene selection using ground truth ---
# The previous method of using low-variance genes failed because strong, pervasive
# batch effects make variance an unreliable indicator of a stable gene.
# In a simulation, we know the ground truth, so we can select controls robustly.

# 1. Get the list of all genes that were NOT affected by the simulated batch effect.
all_genes <- rownames(counts_with_effect)
true_unaffected_genes <- setdiff(all_genes, affected_genes)

# 2. Check if we have enough unaffected genes to choose from.
if (length(true_unaffected_genes) < 50) {
  stop("The number of unaffected genes is less than 50, cannot select a control set.")
}

# 3. Select a random sample of 50 of these true unaffected genes as our controls.
#    This is a stable and reliable method in a simulation context.
set.seed(456) # for reproducibility of control gene selection
ctrl_genes <- sample(true_unaffected_genes, 50)
# --- End of Control Gene Fix ---


members    <- split(seq_along(biological_vars), biological_vars)
max_size   <- max(lengths(members))
scIdx      <- t(sapply(members, function(idxs) c(idxs, rep(-1, max_size - length(idxs)))))

seqset_uq_ruvs <- EDASeq::betweenLaneNormalization(EDASeq::newSeqExpressionSet(round(as.matrix(counts_with_effect))), which="upper")
pheno_ruvg <- data.frame(row.names=colnames(counts_with_effect), x=rep(1, ncol(counts_with_effect)))
seqset_uq_ruvg <- EDASeq::betweenLaneNormalization(EDASeq::newSeqExpressionSet(round(as.matrix(counts_with_effect)), phenoData=Biobase::AnnotatedDataFrame(pheno_ruvg)), which="upper")

# RUV parameter tuning range increased to handle more complex effects
k_range_ruv <- 1:5 
ruvs_metrics_list <- lapply(k_range_ruv, function(k) {
  result <- try({
    res <- run_ruvs(seqset_uq_ruvs, ctrl_genes, scIdx, k = k)
    log_ds <- edgeR::cpm(res$corrected, log=TRUE, prior.count=1)
    if (!is.matrix(log_ds) || any(!is.finite(log_ds))) {
      stop("RUVs produced invalid matrix")
    }
    data.frame(k=k, R2=calculate_permanova_r2(log_ds, batch_info))
  }, silent = TRUE) 
  
  if (inherits(result, "try-error")) {
    cat("Warning: RUVs failed for k =", k, "\n")
    return(NULL)
  }
  return(result)
})
ruvs_metrics_df <- do.call(rbind, ruvs_metrics_list)

if (is.null(ruvs_metrics_df) || nrow(ruvs_metrics_df) == 0) {
  stop("RUVs correction failed for all values of k.")
}

best_k_ruvs <- ruvs_metrics_df$k[which.min(ruvs_metrics_df$R2)]
ruvs_method_name <- paste0("RUVs (k=", best_k_ruvs, ")")
log_ruvs <- edgeR::cpm(run_ruvs(seqset_uq_ruvs, ctrl_genes, scIdx, k = best_k_ruvs)$corrected, log=TRUE, prior.count=1)

ruvg_metrics_list <- lapply(k_range_ruv, function(k) {
  result <- try({
    res <- run_ruvg(seqset_uq_ruvg, ctrl_genes, k = k)
    log_ds <- edgeR::cpm(res, log=TRUE, prior.count=1)
    if (!is.matrix(log_ds) || any(!is.finite(log_ds))) {
      stop("RUVg produced invalid matrix")
    }
    data.frame(k=k, R2=calculate_permanova_r2(log_ds, batch_info))
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    cat("Warning: RUVg failed for k =", k, "\n")
    return(NULL)
  }
  return(result)
})
ruvg_metrics_df <- do.call(rbind, ruvg_metrics_list)

if (is.null(ruvg_metrics_df) || nrow(ruvg_metrics_df) == 0) {
  stop("RUVg correction failed for all values of k.")
}

best_k_ruvg <- ruvg_metrics_df$k[which.min(ruvg_metrics_df$R2)]
ruvg_method_name <- paste0("RUVg (k=", best_k_ruvg, ")")
log_ruvg <- edgeR::cpm(run_ruvg(seqset_uq_ruvg, ctrl_genes, k = best_k_ruvg), log=TRUE, prior.count=1)

pcs_mnn <- run_fastmnn(log_counts = log_counts_before, batch_info = batch_info)
log_pca <- run_pca_correction(log_counts_before, batch_info, sig_threshold = 1e-4)
log_sva <- run_sva(log_counts_before, mod=mod_bio, mod0=mod_null)
log_svaseq <- run_svaseq(counts_with_effect, mod=mod_bio, mod0=mod_null)

cat("--- All correction methods complete. ---\n\n")


# --- 6. Calculate Metrics and Compare All Methods ---
all_method_data <- list(
  list(name = "Pre-Correction", type = "log", data = log_counts_before),
  list(name = "Limma",            type = "log", data = log_limma),
  list(name = "ComBat",           type = "log", data = log_combat),
  list(name = "ComBat-Seq",       type = "log", data = log_seq),
  list(name = ruvs_method_name,   type = "log", data = log_ruvs),
  list(name = ruvg_method_name,   type = "log", data = log_ruvg),
  list(name = "fastMNN",          type = "pca", data = pcs_mnn),
  list(name = "PCA Correction",   type = "log", data = log_pca),
  list(name = "SVA",              type = "log", data = log_sva),
  list(name = "SVA-Seq",          type = "log", data = log_svaseq)
)

# <<< CACHING START for Metrics Calculation >>>
metrics_key <- list(data_key = data_params_key, corrected_data = all_method_data, metrics_version = "v1")
all_metrics <- R.cache::loadCache(key = metrics_key)

if(is.null(all_metrics)) {
  cat("--- No cache found for metrics. Calculating all metrics... ---\n")
  all_metrics_rows <- lapply(all_method_data, function(method) {
    cat("Calculating metrics for:", method$name, "\n")
    is_pca_method <- method$type == "pca"
    
    pca_data_for_metrics <- NULL
    if (is_pca_method) {
      pca_data_for_metrics <- method$data
    } else {
      log_data <- method$data
      gene_vars <- matrixStats::rowVars(as.matrix(log_data))
      log_counts_filtered <- log_data[gene_vars > 1e-8, ]
      if(nrow(log_counts_filtered) >= 2) {
        pca_data_for_metrics <- stats::prcomp(t(log_counts_filtered), scale. = TRUE)$x
      }
    }
    
    if (is.null(pca_data_for_metrics)) {
      return(data.frame(Method = method$name, PERMANOVA_R2 = NA, Batch_Silhouette = NA, Bio_Silhouette = NA, kBET_Rejection = NA, PCR_R2 = NA, iLISI = NA, cLISI = NA, Gene_R2 = NA))
    }
    
    data.frame(
      Method           = method$name,
      PERMANOVA_R2     = calculate_permanova_r2(pca_data = pca_data_for_metrics, variable_info = batch_info),
      Batch_Silhouette = calculate_silhouette_width(pca_data = pca_data_for_metrics, variable_info = batch_info),
      Bio_Silhouette   = calculate_silhouette_width(pca_data = pca_data_for_metrics, variable_info = biological_vars),
      kBET_Rejection   = run_kbet(pca_data = pca_data_for_metrics, batch_info = batch_info),
      PCR_R2           = calculate_pc_regression(pca_data = pca_data_for_metrics, variable_info = batch_info),
      iLISI            = run_lisi(pca_data = pca_data_for_metrics, meta_data = meta_data, var_to_check = "batch"),
      cLISI            = run_lisi(pca_data = pca_data_for_metrics, meta_data = meta_data, var_to_check = "biology"),
      Gene_R2          = if (!is_pca_method) calculate_gene_metric(method$data, batch_info, affected_genes) else NA
    )
  })
  all_metrics <- do.call(rbind, all_metrics_rows)
  R.cache::saveCache(all_metrics, key = metrics_key)
  cat("--- Metrics calculation complete. Saved to cache. ---\n")
} else {
  cat("--- Loaded metrics table from cache. ---\n")
}
# <<< CACHING END for Metrics Calculation >>>


print(knitr::kable(all_metrics, digits = 3, caption = "Comprehensive Batch Correction Method Comparison"))

# --- 7. Final Plots ---

# 7.1 PCA Grid
plots_to_show <- list(
  perform_pca_and_plot(log_counts_before, batch_info, "Pre-Correction"),
  perform_pca_and_plot(log_limma, batch_info, "Limma"),
  perform_pca_and_plot(log_combat, batch_info, "ComBat"),
  perform_pca_and_plot(log_seq, batch_info, "ComBat-Seq"),
  perform_pca_and_plot(log_ruvs, batch_info, ruvs_method_name),
  perform_pca_and_plot(log_ruvg, batch_info, ruvg_method_name),
  perform_pca_and_plot(pca_data=pcs_mnn, batch_info=batch_info, plot_title="fastMNN"),
  perform_pca_and_plot(log_pca, batch_info, "PCA Correction"),
  perform_pca_and_plot(log_sva, batch_info, "SVA"),
  perform_pca_and_plot(log_svaseq, batch_info, "SVA-Seq")
)
plot_panel <- (plots_to_show[[1]]|plots_to_show[[2]]|plots_to_show[[3]])/(plots_to_show[[4]]|plots_to_show[[5]]|plots_to_show[[6]])/(plots_to_show[[7]]|plots_to_show[[8]]|plots_to_show[[9]])/(plots_to_show[[10]]|patchwork::plot_spacer()|patchwork::plot_spacer())
print(plot_panel + plot_layout(guides='collect') & theme(aspect.ratio=1, legend.position="bottom"))

# 7.2 Metrics Barplot
metrics_melted <- tidyr::pivot_longer(all_metrics, cols = -Method, names_to = "Metric", values_to = "Value")

ordered_levels <- all_metrics$Method[!duplicated(all_metrics$Method)]
metrics_melted$Method <- factor(metrics_melted$Method, levels = ordered_levels)

metric_subtitles <- c(
  "PERMANOVA_R2"="Lower is Better", "Batch_Silhouette"="Lower is Better", "Bio_Silhouette"="Higher is Better",
  "kBET_Rejection"="Lower is Better", "PCR_R2"="Lower is Better",
  "iLISI"="Higher is Better (Mixing)", "cLISI"="Lower is Better (Purity)", "Gene_R2"="Lower is Better"
)
metrics_melted$Subtitle <- factor(metric_subtitles[metrics_melted$Metric], levels = unique(metric_subtitles))
print(ggplot(metrics_melted, aes(x=Method, y=Value, fill=Method)) +
        geom_col(show.legend=FALSE) +
        facet_grid(Subtitle ~ Metric, scales="free_y", switch="y") +
        labs(title="Batch Correction Performance Metrics", x=NULL, y=NULL) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle=45, hjust=1, size=9),
          strip.background = element_rect(fill="grey90"),
          strip.text.y.left = element_text(angle=0), panel.spacing = unit(1, "lines")
        )
)