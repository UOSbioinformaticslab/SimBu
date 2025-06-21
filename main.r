# main.R

# --- 0. Setup ---
# Source the file containing our custom functions
source("functions.R")

# Define high-level parameters for the analysis
# Add RUVSeq to the list of packages
PACKAGES <- c("SimBu", "limma", "edgeR", "dplyr", "ggplot2", "patchwork", "SummarizedExperiment", 
              "Matrix", "sva", "RUVSeq", "EDASeq", "vegan", "cluster", "tidyr", 
              "knitr", "batchelor", "SingleCellExperiment", "cachem") # Added batchelor
NSAMPLES_PER_BATCH <- 30
N_GENES <- 1000
N_CELLS <- 300


# --- 1. Environment Setup ---
# Install and load all required packages
install_and_load_packages(PACKAGES)


# --- 2. Data Generation and Simulation ---
set.seed(123)
base_dataset <- create_base_sc_dataset(n_genes = N_GENES, n_cells = N_CELLS)
simulation_list <- generate_simulations(
  simbu_dataset = base_dataset,
  n_sims = 3,
  n_samples_per_sim = NSAMPLES_PER_BATCH
)


# --- 3. Introduce Batch Effect & Prepare Data ---
# Merge simulations, add batch effect, and get list of affected genes
batch_effect_data <- introduce_batch_effect(simulation_list)
counts_with_effect <- assay(batch_effect_data$se_with_batch, "bulk_counts")
batch_info <- batch_effect_data$batch_info
affected_genes <- batch_effect_data$affected_genes # Get the affected genes

# Calculate log-counts for methods that require it
log_counts_before <- edgeR::cpm(as.matrix(counts_with_effect), log = TRUE, prior.count = 1)

# Generate PCA plot *before* batch correction
plot_before <- perform_pca_and_plot(
  log_counts = log_counts_before,
  batch_info = batch_info,
  plot_title = "Pre-Correction"
)


# --- 4. Initialize Results Storage ---
# We will store plots and metrics here
plots_list <- list()
metrics_list <- list()
plots_list[["Pre-Correction"]] <- plot_before
metrics_list[["Pre-Correction"]] <- data.frame(
  Method = "Pre-Correction",
  R_Squared = calculate_permanova_r2(log_counts_before, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_counts_before, batch_info)
)


# --- 5. Run & Evaluate Correction Methods ---

# Method 1: limma
log_counts_limma <- limma::removeBatchEffect(log_counts_before, batch = batch_info)
plots_list[["Limma"]] <- perform_pca_and_plot(log_counts_limma, batch_info, "Post-Correction (limma)")
metrics_list[["Limma"]] <- data.frame(
  Method = "Limma",
  R_Squared = calculate_permanova_r2(log_counts_limma, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_counts_limma, batch_info)
)

# Method 2: ComBat
log_counts_combat <- run_combat(log_counts_before, batch_info)
plots_list[["ComBat"]] <- perform_pca_and_plot(log_counts_combat, batch_info, "Post-Correction (ComBat)")
metrics_list[["ComBat"]] <- data.frame(
  Method = "ComBat",
  R_Squared = calculate_permanova_r2(log_counts_combat, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_counts_combat, batch_info)
)

# Method 3: ComBat-Seq
counts_combat_seq <- run_combat_seq(counts_with_effect, batch_info)
log_counts_combat_seq <- edgeR::cpm(counts_combat_seq, log = TRUE, prior.count = 1)
plots_list[["ComBat-Seq"]] <- perform_pca_and_plot(log_counts_combat_seq, batch_info, "Post-Correction (ComBat-Seq)")
metrics_list[["ComBat-Seq"]] <- data.frame(
  Method = "ComBat-Seq",
  R_Squared = calculate_permanova_r2(log_counts_combat_seq, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_counts_combat_seq, batch_info)
)

# Method 4: RUVg
unaffected_genes <- setdiff(rownames(counts_with_effect), affected_genes)
set.seed(999)
ruvg_control_genes <- sample(unaffected_genes, 500)
counts_ruvg <- run_ruvg(counts_with_effect, ruvg_control_genes, k = 1)
log_counts_ruvg <- edgeR::cpm(counts_ruvg, log = TRUE, prior.count = 1)
plots_list[["RUVg"]] <- perform_pca_and_plot(log_counts_ruvg, batch_info, "Post-Correction (RUVg)")
metrics_list[["RUVg"]] <- data.frame(
  Method = "RUVg",
  R_Squared = calculate_permanova_r2(log_counts_ruvg, batch_info),
  Silhouette_Width = calculate_silhouette_width(log_counts_ruvg, batch_info)
)

# Method 5: fastMNN
corrected_pcs_mnn <- run_fastmnn(log_counts_before, batch_info)

# Create a data frame for plotting with the correct column names
plot_data_mnn <- data.frame(
  PC1 = corrected_pcs_mnn[, 1],
  PC2 = corrected_pcs_mnn[, 2],
  batch = batch_info
)

# Pass the pre-computed PCs directly to the plotting and metrics functions
plots_list[["fastMNN"]] <- perform_pca_and_plot(
  pca_data = plot_data_mnn, 
  batch_info = batch_info, # Still need to pass this for color
  plot_title = "Post-Correction (fastMNN)"
)
metrics_list[["fastMNN"]] <- data.frame(
  Method = "fastMNN",
  R_Squared = calculate_permanova_r2(pca_data = corrected_pcs_mnn, batch_info = batch_info),
  Silhouette_Width = calculate_silhouette_width(pca_data = corrected_pcs_mnn, batch_info = batch_info)
)

# --- 6. Final Output ---

# --- 6.1. Log2CPM Distribution Boxplots ---
cat("\n--- Preparing data for Log2CPM distribution plot ---\n")

# Collect all log-CPM matrices into a list (excluding fastMNN)
logcpm_list <- list(
  "Pre-Correction" = log_counts_before,
  "Limma" = log_counts_limma,
  "ComBat" = log_counts_combat,
  "ComBat-Seq" = log_counts_combat_seq,
  "RUVg" = log_counts_ruvg
)

# Create a sample-to-batch lookup table
sample_info <- data.frame(Sample_ID = colnames(log_counts_before), Batch = batch_info)

# Process each matrix and convert to a tidy format
tidy_df_list <- lapply(names(logcpm_list), function(method_name) {
  
  # Convert matrix to a long-format data frame
  logcpm_matrix <- logcpm_list[[method_name]]
  
  logcpm_df <- as.data.frame(logcpm_matrix) %>%
    # Add gene names as a column before pivoting
    dplyr::mutate(Gene = rownames(logcpm_matrix)) %>%
    # Pivot from wide to long format
    tidyr::pivot_longer(
      cols = -Gene,
      names_to = "Sample_ID",
      values_to = "Log2CPM"
    ) %>%
    # Add the method name
    dplyr::mutate(Method = method_name) %>%
    # Add batch information by joining with the lookup table
    dplyr::left_join(sample_info, by = "Sample_ID")
  
  return(logcpm_df)
})

# Combine all tidy data frames into one
full_tidy_df <- dplyr::bind_rows(tidy_df_list)

# Generate and print the box plot
log2cpm_boxplot <- generate_log2cpm_boxplot(full_tidy_df)
print(log2cpm_boxplot)

# --- 6.2. Visual PCA Comparison ---
# Arrange the 6 plots into a 2x3 grid
plot_panel <- (plots_list[[1]] | plots_list[[2]] | plots_list[[3]]) /
  (plots_list[[4]] | plots_list[[5]] | plots_list[[6]])

final_pca_plot <- plot_panel + 
  plot_layout(guides = 'collect') & 
  theme(aspect.ratio = 1, legend.position = "bottom")

print(final_pca_plot)

# --- 6.3. Quantitative Metrics Table ---
metrics_df <- do.call(rbind, metrics_list)
rownames(metrics_df) <- NULL
cat("\n--- Quantitative Batch Correction Metrics ---\n")
print(knitr::kable(metrics_df, digits = 3, caption = "Lower R² and Silhouette Width are better."))

# --- 6.4. Quantitative Metrics Plot ---
metrics_melted <- tidyr::pivot_longer(
  metrics_df,
  cols = c("R_Squared", "Silhouette_Width"),
  names_to = "Metric",
  values_to = "Value"
)

# Update the method order for plotting
method_order <- c("Pre-Correction", "Limma", "ComBat", "ComBat-Seq", "RUVg", "fastMNN")
metrics_melted$Method <- factor(metrics_melted$Method, levels = method_order)

metrics_plot <- ggplot(metrics_melted, aes(x = Method, y = Value, fill = Method)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    title = "Quantitative Comparison of Batch Correction Methods",
    x = "Method",
    y = "Metric Value (Lower is Better)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(metrics_plot)