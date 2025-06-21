# main.R

# --- 0. Setup ---
# Source the file containing our custom functions
source("functions.R")

# Define high-level parameters for the analysis
PACKAGES <- c("SimBu", "limma", "edgeR", "ggplot2", "patchwork", "SummarizedExperiment", "Matrix")
NSAMPLES_PER_BATCH <- 30
N_GENES <- 1000
N_CELLS <- 300


# --- 1. Environment Setup ---
# Install and load all required packages
install_and_load_packages(PACKAGES)


# --- 2. Data Generation and Simulation ---
set.seed(123)
# Create a base single-cell dataset
base_dataset <- create_base_sc_dataset(n_genes = N_GENES, n_cells = N_CELLS)

# Generate 3 separate simulation batches from the base dataset
simulation_list <- generate_simulations(
  simbu_dataset = base_dataset,
  n_sims = 3,
  n_samples_per_sim = NSAMPLES_PER_BATCH
)


# --- 3. Introduce and Visualize Batch Effect ---
# Merge simulations and add a strong batch effect to the count data
batch_effect_data <- introduce_batch_effect(simulation_list)
counts_with_effect <- assay(batch_effect_data$se_with_batch, "bulk_counts")
batch_info <- batch_effect_data$batch_info

# Calculate log-counts for PCA
log_counts_before <- cpm(counts_with_effect, log = TRUE, prior.count = 1)

# Generate PCA plot *before* batch correction
plot_before <- perform_pca_and_plot(
  log_counts = log_counts_before,
  batch_info = batch_info,
  plot_title = "Pre Batch Correction"
)


# --- 4. Batch Correction and Visualization ---
# Correct for the known batch effect using limma's removeBatchEffect
# Note: removeBatchEffect requires a design matrix for covariates to keep,
# here we have none, so we use a constant intercept model.
design <- model.matrix(~1, data = colData(batch_effect_data$se_with_batch))
log_counts_after <- limma::removeBatchEffect(log_counts_before, batch = batch_info, design = design)

# Generate PCA plot *after* batch correction
plot_after <- perform_pca_and_plot(
  log_counts = log_counts_after,
  batch_info = batch_info,
  plot_title = "Post removeBatchEffect Correction"
)


# --- 5. Final Output ---
# Combine the "before" and "after" plots for a side-by-side comparison
final_plot <- (plot_before + plot_after + plot_layout(guides = 'collect')) &
  theme(aspect.ratio = 1)

# Print the final plot to the display
print(final_plot)
