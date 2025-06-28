# main.R

# --- 0. Setup ---
source("functions.R")

# High-level parameters
PACKAGES <- c(
  "SimBu", "limma",   "edgeR", "dplyr",  "ggplot2",
  "patchwork", "SummarizedExperiment", "Matrix",
  "sva",  "RUVSeq", "EDASeq", "vegan",
  "cluster", "tidyr", "knitr", "batchelor",
  "SingleCellExperiment", "gtools", "R.cache",
  "matrixStats"
)
NSAMPLES_PER_BATCH <- 30
N_GENES            <- 1000
N_CELLS            <- 300
TOP_GENES          <- 200
NUM_SIMS           <- 3
# Parameters for create_base_sc_dataset (clinically tuned)
MEAN_SHAPE   <- 1.5
MEAN_SCALE   <- 0.5
SIZE_MEANLOG <- 0     # log-normal mean ~1
SIZE_SDLOG   <- 0.05  # modest library variability
DISPERSION   <- 0.2   # NB dispersion (1/size)
# Batch effect parameters
N_AFFECTED_GENES  <- 50
EFFECT_MULTIPLIER <- 7.0


# --- 1. Environment Setup ---
install_and_load_packages(PACKAGES)

# --- 2. Data Generation and Simulation ---
set.seed(123)
base_ds <- create_base_sc_dataset(
  n_genes     = N_GENES,
  n_cells     = N_CELLS,
  mean_shape  = MEAN_SHAPE,
  mean_scale  = MEAN_SCALE,
  size_shape  = SIZE_SHAPE,
  size_scale  = SIZE_SCALE,
  dispersion  = DISPERSION
)

# Filter to most specific genes
filtered_ds <- filter_genes_by_specificity(base_ds, top_n = TOP_GENES)

# Generate bulk RNA-seq simulations
sim_list <- generate_simulations(
  simbu_dataset     = filtered_ds,
  n_sims            = NUM_SIMS,
  n_samples_per_sim = NSAMPLES_PER_BATCH
)

# --- 3. Introduce Batch Effect & Prepare Data ---
batch_data       <- introduce_batch_effect(
  sim_list,
  n_genes_affected  = N_AFFECTED_GENES,
  effect_multiplier = EFFECT_MULTIPLIER
)
counts_eff       <- assay(batch_data$se_with_batch, "bulk_counts")
batch_info       <- batch_data$batch_info
affected_genes   <- batch_data$affected_genes

# Compute log2-CPM for pre-correction
df_pre           <- edgeR::cpm(
  as.matrix(counts_eff),
  log         = TRUE,
  prior.count = 1
)

# --- 4. Initialize Storage for Plots & Metrics ---
plots_list <- list(
  Pre = perform_pca_and_plot(
    log_counts = df_pre,
    batch_info = batch_info,
    plot_title = "Pre-Correction"
  )
)
metrics_list <- list(
  Pre = data.frame(
    Method           = "Pre-Correction",
    R_Squared        = calculate_permanova_r2(df_pre, batch_info),
    Silhouette_Width = calculate_silhouette_width(df_pre, batch_info)
  )
)

# --- 5. Correction Methods ---
run_method <- function(name, corrected, type = c("log", "pca")) {
  type <- match.arg(type)
  if (type == "log") {
    plot_obj <- perform_pca_and_plot(
      log_counts = corrected,
      batch_info = batch_info,
      plot_title = name
    )
    r2_val <- calculate_permanova_r2(corrected, batch_info)
    sw_val <- calculate_silhouette_width(corrected, batch_info)
  } else {
    plot_obj <- perform_pca_and_plot(
      pca_data   = corrected,
      batch_info = batch_info,
      plot_title = name
    )
    r2_val <- calculate_permanova_r2(NULL, batch_info, pca_data = corrected)
    sw_val <- calculate_silhouette_width(NULL, batch_info, pca_data = corrected)
  }
  plots_list[[name]] <<- plot_obj
  metrics_list[[name]] <<- data.frame(
    Method           = name,
    R_Squared        = r2_val,
    Silhouette_Width = sw_val
  )
}

# --- 5. Run & Evaluate Correction Methods ---

## 5.1 limma
log_limma <- limma::removeBatchEffect(log_counts_before, batch = batch_info)
plots_list$Limma <- perform_pca_and_plot(log_limma, batch_info, "Limma")
metrics_list$Limma <- data.frame(
  Method           = "Limma",
  R_Squared        = calculate_permanova_r2(log_limma, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_limma, batch_info)
)

## 5.2 ComBat
log_combat <- run_combat(log_counts_before, batch_info, ref_batch = "Batch1")
plots_list$ComBat <- perform_pca_and_plot(log_combat, batch_info, "ComBat")
metrics_list$ComBat <- data.frame(
  Method           = "ComBat",
  R_Squared        = calculate_permanova_r2(log_combat, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_combat, batch_info)
)

## 5.3 ComBat-Seq
counts_seq <- run_combat_seq(counts_with_effect, batch_info)
log_seq    <- edgeR::cpm(counts_seq, log = TRUE, prior.count = 1)
plots_list$`ComBat-Seq` <- perform_pca_and_plot(log_seq, batch_info, "ComBat-Seq")
metrics_list$`ComBat-Seq` <- data.frame(
  Method           = "ComBat-Seq",
  R_Squared        = calculate_permanova_r2(log_seq, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_seq, batch_info)
)

## 5.4 RUVr (dense + phenoData + tuning)
counts_mat <- as.matrix(counts_with_effect)
pheno_df   <- data.frame(batch = batch_info,
                         row.names = colnames(counts_mat))
phenoData  <- Biobase::AnnotatedDataFrame(pheno_df)

seqset <- EDASeq::newSeqExpressionSet(counts = counts_mat,
                                      phenoData = phenoData)

logcpm_batch <- edgeR::cpm(counts_mat, log=TRUE, prior.count=1)
gene_var     <- matrixStats::rowVars(logcpm_batch)
ctrl_cands   <- names(sort(gene_var))[1:50]
ctrl_genes   <- setdiff(ctrl_cands, affected_genes)
if (length(ctrl_genes) < 10) {
  warning("Fewer than 10 true controls; RUVr may underperform.")
}

dge    <- edgeR::DGEList(counts = counts_mat)
dge    <- edgeR::calcNormFactors(dge)
design <- model.matrix(~ batch_info)
dge    <- edgeR::estimateDisp(dge, design)
fit    <- edgeR::glmFit(dge, design)
resid_mat <- residuals(fit, type = "deviance")
# tune k = 1:3
ruvr_metrics <- lapply(1:3, function(kval) {
  seq_ruv  <- RUVSeq::RUVr(x         = seqset,
                           cIdx      = ctrl_genes,
                           k         = kval,
                           residuals = resid_mat)
  corr_mat <- Biobase::assayDataElement(seq_ruv, "normalizedCounts")
  log_ruvr <- edgeR::cpm(corr_mat, log=TRUE, prior.count=1)
  data.frame(
    k          = kval,
    R2         = calculate_permanova_r2(log_ruvr, batch_info),
    Silhouette = calculate_silhouette_width(log_ruvr, batch_info)
  )
})
ruvr_metrics_df <- do.call(rbind, ruvr_metrics)
print(knitr::kable(
  ruvr_metrics_df,
  caption = "RUVr Tuning: k vs Batch R² & Silhouette"
))
best_k   <- ruvr_metrics_df$k[which.min(ruvr_metrics_df$R2)]
seq_best <- RUVSeq::RUVr(x         = seqset,
                         cIdx      = ctrl_genes,
                         k         = best_k,
                         residuals = resid_mat)
corr_best <- Biobase::assayDataElement(seq_best, "normalizedCounts")
log_best  <- edgeR::cpm(corr_best, log=TRUE, prior.count=1)

plots_list$RUVr <- perform_pca_and_plot(
  log_best, batch_info, paste0("RUVr (k=", best_k, ")")
)
metrics_list$RUVr <- data.frame(
  Method           = paste0("RUVr (k=", best_k, ")"),
  R_Squared        = calculate_permanova_r2(log_best, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_best, batch_info)
)

## 5.5 fastMNN
pcs_mnn <- run_fastmnn(log_counts = log_counts_before, batch_info = batch_info)
plots_list$fastMNN <- perform_pca_and_plot(pca_data = pcs_mnn,
                                           batch_info = batch_info,
                                           plot_title = "fastMNN")
metrics_list$fastMNN <- data.frame(
  Method           = "fastMNN",
  R_Squared        = calculate_permanova_r2(NULL, batch_info, pca_data = pcs_mnn),
  Silhouette_Width = calculate_silhouette_width(NULL, batch_info, pca_data = pcs_mnn)
)

## 5.6 PCA-Correction
log_pca <- run_pca_correction(log_counts_before, batch_info, sig_threshold = 1e-4)
plots_list$PCA <- perform_pca_and_plot(log_pca, batch_info, "PCA Correction")
metrics_list$PCA <- data.frame(
  Method           = "PCA Correction",
  R_Squared        = calculate_permanova_r2(log_pca, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_pca, batch_info)
)

# --- 6. Compare All Methods: add kBET (k₀ = 10) ---
method_outputs <- list(
  Pre          = list(type="log", mat=log_counts_before),
  Limma        = list(type="log", mat=log_limma),
  ComBat       = list(type="log", mat=log_combat),
  `ComBat-Seq` = list(type="log", mat=log_seq),
  RUVr         = list(type="log", mat=log_best),
  fastMNN      = list(type="pca", mat=pcs_mnn),
  PCA          = list(type="log", mat=log_pca)
)

rows <- lapply(names(method_outputs), function(key) {
  out  <- method_outputs[[key]]
  base <- metrics_list[[key]]
  
  if (out$type == "log") {
    pr  <- prcomp(t(out$mat), scale.=TRUE)
    pcs <- pr$x[, 1:min(10, ncol(pr$x))]
  } else {
    pcs <- out$mat[, 1:min(10, ncol(out$mat))]
  }
  
  kbet_val <- tryCatch({
    res <- kbet::kbet(df    = pcs,
                      batch = batch_info,
                      k0    = 10,
                      plot  = FALSE)
    res$summary$kbet.observed
  }, error = function(e) NA)
  
  data.frame(
    Method           = base$Method,
    R_Squared        = base$R_Squared,
    Silhouette_Width = base$Silhouette_Width,
    kBET_Rejection   = kbet_val,
    stringsAsFactors = FALSE
  )
})

all_metrics <- do.call(rbind, rows)
print(knitr::kable(
  all_metrics,
  digits  = 3,
  caption = "All Methods: PERMANOVA R², Silhouette Width, and kBET (k₀=10)"
))

# --- 7. Final Plots ---

## 7.1 Log2CPM Boxplots
cat("\n--- Preparing data for Log2CPM distribution plot ---\n")
logcpm_list <- list(
  "Pre-Correction" = log_counts_before,
  "Limma"          = log_limma,
  "ComBat"         = log_combat,
  "ComBat-Seq"     = log_seq,
  "RUVr"           = log_best,
  "PCA Correction" = log_pca
)
sample_info <- data.frame(Sample_ID = colnames(log_counts_before),
                          Batch     = batch_info,
                          stringsAsFactors = FALSE)

tidy_df_list <- lapply(names(logcpm_list), function(meth) {
  df <- as.data.frame(logcpm_list[[meth]])
  df$Gene <- rownames(df)
  tidyr::pivot_longer(df, -Gene, names_to="Sample_ID", values_to="Log2CPM") %>%
    dplyr::mutate(Method = meth) %>%
    dplyr::left_join(sample_info, by="Sample_ID")
})
full_tidy_df <- dplyr::bind_rows(tidy_df_list)
print(generate_log2cpm_boxplot(full_tidy_df))

## 7.2 PCA Grid
plot_panel <- (plots_list[[1]] | plots_list[[2]] | plots_list[[3]]) /
  (plots_list[[4]] | plots_list[[5]] | plots_list[[6]]) /
  (plots_list[[7]] | plot_spacer()     | plot_spacer())
print(plot_panel + plot_layout(guides='collect') &
        theme(aspect.ratio=1, legend.position="bottom"))

## 7.3 Metrics Barplot
metrics_df <- do.call(rbind, metrics_list)
metrics_melted <- tidyr::pivot_longer(
  metrics_df, cols = c("R_Squared","Silhouette_Width"),
  names_to="Metric", values_to="Value"
)
metrics_melted$Method <- factor(
  metrics_melted$Method,
  levels = c(
    "Pre-Correction","Limma","ComBat",
    "ComBat-Seq",paste0("RUVr (k=",best_k,")"),
    "fastMNN","PCA Correction"
  )
)
print(ggplot2::ggplot(
  metrics_melted, aes(x=Method,y=Value,fill=Method)
) +
  ggplot2::geom_col(show.legend=FALSE) +
  ggplot2::facet_wrap(~Metric, scales="free_y") +
  ggplot2::labs(
    title="Batch Correction Performance",
    subtitle="Lower values indicate better removal of batch-driven variation",
    x=NULL,y="Metric"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.background = element_rect(fill="grey90")
  )
)
