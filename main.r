# --- 0. Setup ---
# Define high-level parameters
PACKAGES <- c(
  "SimBu", "limma",     "edgeR", "dplyr",   "ggplot2",
  "patchwork", "SummarizedExperiment", "Matrix", "tibble",
  "sva",   "RUVSeq", "EDASeq", "vegan",
  "cluster", "tidyr", "knitr", "batchelor",
  "SingleCellExperiment", "gtools", "R.cache",
  "matrixStats", "lisi", "remotes", "welch-lab/kBET"
)


NSAMPLES_PER_BATCH <- 50
N_GENES            <- 1000
N_CELLS            <- 300
# --- Adjusting gene counts to be valid ---
TOP_GENES         <- 200
N_AFFECTED_GENES  <- 100 
N_DE_GENES        <- 40 
# --- Parameter Adjustments for Realistic Batch Effects ---
NUM_SIMS           <- 4
MEAN_SHAPE   <- 1.5
MEAN_SCALE   <- 0.5
SIZE_MEANLOG <- 0
SIZE_SDLOG   <- 0.5
DISPERSION   <- 0.4
EFFECT_MULTIPLIER <- 3.5
DE_FOLD_CHANGE <- 2.5


# --- 1. Environment Setup ---
install_and_load_packages(PACKAGES)
set.seed(123) # Set seed for reproducibility

# --- 2. Data Generation and Simulation ---
# <<< CACHING START for Data Generation >>>
cat("--- Checking cache for prepared dataset... ---\n")
data_params_key <- list(
  N_GENES=N_GENES, N_CELLS=N_CELLS, TOP_GENES=TOP_GENES, NUM_SIMS=NUM_SIMS,
  NSAMPLES_PER_BATCH=NSAMPLES_PER_BATCH, EFFECT_MULTIPLIER=EFFECT_MULTIPLIER,
  N_AFFECTED_GENES=N_AFFECTED_GENES, DISPERSION=DISPERSION, 
  N_DE_GENES=N_DE_GENES, DE_FOLD_CHANGE=DE_FOLD_CHANGE,
  script_version="v8_pca_shapes"
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
    sim_list, n_genes_affected=N_AFFECTED_GENES, effect_multiplier=EFFECT_MULTIPLIER,
    n_de_genes = N_DE_GENES, de_fold_change = DE_FOLD_CHANGE
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
names(batch_info)  <- colnames(counts_with_effect)
biological_vars    <- prepared_data$biological_vars
affected_genes     <- prepared_data$affected_genes
true_de_genes      <- prepared_data$true_de_genes
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
log_combat <- run_combat(log_counts_before, batch_info, mod=mod_combat)
log_combat_ref <- run_combat(log_counts_before, batch_info, mod=mod_combat, ref_batch = "Batch1")
counts_seq <- run_combat_seq(counts_with_effect, batch_info, biological_vars = biological_vars)
log_seq    <- edgeR::cpm(counts_seq, log = TRUE, prior.count = 1)

cat("--- Preparing for RUV Methods ---\n")
true_unaffected_genes <- setdiff(rownames(counts_with_effect), c(affected_genes, true_de_genes))
if (length(true_unaffected_genes) < 50) stop("Not enough unaffected genes for ideal RUV control set.")
set.seed(456)
ideal_ctrl_genes <- sample(true_unaffected_genes, 50)
empirical_ctrl_genes <- get_empirical_controls(log_counts_before, batch_info, n_controls = 50)
cat("Using", length(empirical_ctrl_genes), "empirical control genes for RUV.\n")

members    <- split(seq_along(biological_vars), biological_vars)
max_size   <- max(lengths(members))
scIdx      <- t(sapply(members, function(idxs) c(idxs, rep(-1, max_size - length(idxs)))))
seqset_uq <- EDASeq::betweenLaneNormalization(EDASeq::newSeqExpressionSet(round(as.matrix(counts_with_effect))), which="upper")
pheno_ruvg <- data.frame(row.names=colnames(counts_with_effect), x=rep(1, ncol(counts_with_effect)))
seqset_uq_ruvg <- EDASeq::betweenLaneNormalization(EDASeq::newSeqExpressionSet(round(as.matrix(counts_with_effect)), phenoData=Biobase::AnnotatedDataFrame(pheno_ruvg)), which="upper")
k_range_ruv <- 1:5

cat("Tuning RUVs (Ideal Controls)...\n")
ruvs_ideal_log <- tune_and_run_ruv(ruv_func = run_ruvs, seqset = seqset_uq, 
                                   ctrl_genes = ideal_ctrl_genes, scIdx = scIdx, k_range = k_range_ruv, 
                                   batch_info = batch_info)
cat("Tuning RUVg (Ideal Controls)...\n")
ruvg_ideal_log <- tune_and_run_ruv(ruv_func = run_ruvg, seqset = seqset_uq_ruvg, 
                                   ctrl_genes = ideal_ctrl_genes, scIdx = NULL, k_range = k_range_ruv, 
                                   batch_info = batch_info)
cat("Tuning RUVs (Empirical Controls)...\n")
ruvs_empirical_log <- tune_and_run_ruv(ruv_func = run_ruvs, seqset = seqset_uq, 
                                       ctrl_genes = empirical_ctrl_genes, scIdx = scIdx, k_range = k_range_ruv, 
                                       batch_info = batch_info)
cat("Tuning RUVg (Empirical Controls)...\n")
ruvg_empirical_log <- tune_and_run_ruv(ruv_func = run_ruvg, seqset = seqset_uq_ruvg, 
                                       ctrl_genes = empirical_ctrl_genes, scIdx = NULL, k_range = k_range_ruv, 
                                       batch_info = batch_info)

pcs_mnn <- run_fastmnn(log_counts = log_counts_before, batch_info = batch_info)
log_pca <- run_pca_correction(log_counts_before, batch_info, sig_threshold = 1e-4)
log_sva <- run_sva(log_counts_before, mod=mod_bio, mod0=mod_null)
log_svaseq <- run_svaseq(counts_with_effect, mod=mod_bio, mod0=mod_null)

cat("--- All correction methods complete. ---\n\n")


# --- 6. Calculate Metrics and Compare All Methods ---
all_method_data <- list(
  list(name = "Pre-Correction",     type = "log", data = log_counts_before),
  list(name = "Limma",                type = "log", data = log_limma),
  list(name = "ComBat",               type = "log", data = log_combat),
  list(name = "ComBat-Ref",           type = "log", data = log_combat_ref),
  list(name = "ComBat-Seq",           type = "log", data = log_seq),
  list(name = "RUVs (Ideal)",         type = "log", data = ruvs_ideal_log$data),
  list(name = "RUVg (Ideal)",         type = "log", data = ruvg_ideal_log$data),
  list(name = "RUVs (Empirical)",     type = "log", data = ruvs_empirical_log$data),
  list(name = "RUVg (Empirical)",     type = "log", data = ruvg_empirical_log$data),
  list(name = "fastMNN",              type = "pca", data = pcs_mnn),
  list(name = "PCA Correction",       type = "log", data = log_pca),
  list(name = "SVA",                  type = "log", data = log_sva),
  list(name = "SVA-Seq",              type = "log", data = log_svaseq)
)

# <<< CACHING START for Metrics Calculation >>>
metrics_key <- list(data_key = data_params_key, corrected_data = all_method_data, metrics_version = "v2_de_metrics")
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
      return(data.frame(Method = method$name, PERMANOVA_R2 = NA, Batch_Silhouette = NA, Bio_Silhouette = NA, kBET_Rejection = NA, PCR_R2 = NA, iLISI = NA, cLISI = NA, Gene_R2 = NA, TPR=NA, FPR=NA))
    }
    
    de_metrics <- if (!is_pca_method) {
      run_de_and_evaluate(method$data, biological_vars, true_de_genes)
    } else {
      list(TPR = NA, FPR = NA)
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
      Gene_R2          = if (!is_pca_method) calculate_gene_metric(method$data, batch_info, affected_genes) else NA,
      TPR              = de_metrics$TPR,
      FPR              = de_metrics$FPR
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
pca_plots_list <- lapply(all_method_data, function(method) {
  # Pass biological_vars to the plotting function
  if (method$type == "log") {
    perform_pca_and_plot(log_counts = method$data, batch_info = batch_info, 
                         biological_vars = biological_vars, plot_title = method$name)
  } else {
    perform_pca_and_plot(pca_data = method$data, batch_info = batch_info, 
                         biological_vars = biological_vars, plot_title = method$name)
  }
})
while(length(pca_plots_list) < 15) pca_plots_list <- c(pca_plots_list, list(patchwork::plot_spacer()))
plot_panel_pca <- patchwork::wrap_plots(pca_plots_list, ncol = 5)
print(plot_panel_pca + plot_layout(guides='collect') & theme(aspect.ratio=1, legend.position="bottom"))


# 7.2 Log2CPM Boxplot Grid
boxplot_grid_plot <- plot_log2cpm_boxplot_grid(all_method_data, batch_info)
print(boxplot_grid_plot)


# 7.3 Metrics Barplot
metrics_melted <- tidyr::pivot_longer(all_metrics, cols = -Method, names_to = "Metric", values_to = "Value")
ordered_levels <- all_metrics$Method[!duplicated(all_metrics$Method)]
metrics_melted$Method <- factor(metrics_melted$Method, levels = ordered_levels)
metric_subtitles <- c(
  "PERMANOVA_R2"="Lower is Better", "Batch_Silhouette"="Lower is Better", "Bio_Silhouette"="Higher is Better",
  "kBET_Rejection"="Lower is Better", "PCR_R2"="Lower is Better",
  "iLISI"="Higher is Better (Mixing)", "cLISI"="Lower is Better (Purity)", "Gene_R2"="Lower is Better",
  "TPR"="Higher is Better", "FPR"="Lower is Better"
)
metrics_melted$Subtitle <- factor(metric_subtitles[metrics_melted$Metric], levels = unique(metric_subtitles))

metrics_melted <- metrics_melted %>% filter(!is.na(Value))

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