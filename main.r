# main.R

# --- 0. Setup ---
# Source the file containing our custom functions
source("functions.R")

# Define high-level parameters for the analysis
# Add RUVSeq to the list of packages
PACKAGES <- c("SimBu", "limma", "edgeR", "ggplot2", "patchwork", "SummarizedExperiment", "Matrix", "sva", "RUVSeq", "EDASeq")
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
log_counts_before <- cpm(counts_with_effect, log = TRUE, prior.count = 1)

# Generate PCA plot *before* batch correction
plot_before <- perform_pca_and_plot(
  log_counts = log_counts_before,
  batch_info = batch_info,
  plot_title = "Pre-Correction"
)


# --- 4. Correction Method 1: limma ---
design_limma <- model.matrix(~1, data = colData(batch_effect_data$se_with_batch))
log_counts_after_limma <- limma::removeBatchEffect(log_counts_before, batch = batch_info, design = design_limma)
plot_after_limma <- perform_pca_and_plot(
  log_counts = log_counts_after_limma, batch_info, "Post-Correction (limma)"
)


# --- 5. Correction Method 2: ComBat (Original) ---
log_counts_after_combat <- run_combat(log_counts = log_counts_before, batch_info = batch_info)
plot_after_combat <- perform_pca_and_plot(
  log_counts = log_counts_after_combat, batch_info, "Post-Correction (ComBat)"
)


# --- 6. Correction Method 3: ComBat-Seq (on raw counts) ---
counts_after_combat_seq <- run_combat_seq(raw_counts = counts_with_effect, batch_info = batch_info)
log_counts_after_combat_seq <- cpm(counts_after_combat_seq, log = TRUE, prior.count = 1)
plot_after_combat_seq <- perform_pca_and_plot(
  log_counts = log_counts_after_combat_seq, batch_info, "Post-Correction (ComBat-Seq)"
)


# --- 7. Correction Method 4: RUVg (on raw counts) ---
# Define control genes as those NOT affected by our artificial batch effect
unaffected_genes <- setdiff(rownames(counts_with_effect), affected_genes)
set.seed(999) # For reproducible sampling of control genes
ruvg_control_genes <- sample(unaffected_genes, 500) # Use a subset of 500

# Run RUVg
counts_after_ruvg <- run_ruvg(
  raw_counts = counts_with_effect,
  control_genes = ruvg_control_genes,
  k = 1 # We simulated 1 major batch effect factor
)
log_counts_after_ruvg <- cpm(counts_after_ruvg, log = TRUE, prior.count = 1)
plot_after_ruvg <- perform_pca_and_plot(
  log_counts = log_counts_after_ruvg, batch_info, "Post-Correction (RUVg)"
)


# --- 8. Final Output ---
# Combine all plots into a 2x3 panel for comparison
# plot_spacer() creates a blank plot to fill the grid
plot_panel <- plot_before + plot_after_limma + plot_after_combat +
  plot_after_combat_seq + plot_after_ruvg + plot_spacer()

final_plot <- plot_panel + 
  plot_layout(ncol = 3, guides = 'collect') & 
  theme(aspect.ratio = 1, legend.position = "bottom")

# Print the final plot
print(final_plot)

