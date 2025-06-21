# --- 1. Installation of Required Packages ---
# This section ensures all necessary packages are installed.
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of packages to install
packages_to_install <- c("SimBu", "limma", "edgeR", "ggplot2", "patchwork")
for (pkg in packages_to_install) {
  if (!require(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- 2. Load Libraries ---
library(SimBu)
library(SummarizedExperiment) # For accessing colData
library(limma)              # For removeBatchEffect
library(edgeR)              # For cpm function (log-transformation)
library(ggplot2)            # For plotting
library(patchwork)          # For combining plots

# --- 3. Create Base Single-Cell Dataset ---
# This part is similar to your original code, setting up the reference data.
counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
colnames(counts) <- paste0("cell_", rep(1:300))
colnames(tpm) <- paste0("cell_", rep(1:300))
rownames(counts) <- paste0("gene_", rep(1:1000))
rownames(tpm) <- paste0("gene_", rep(1:1000))

annotation <- data.frame(
  "ID" = paste0("cell_", rep(1:300)),
  "cell_type" = c(
    rep("T cells CD4", 50),
    rep("T cells CD8", 50),
    rep("Macrophages", 100),
    rep("NK cells", 10),
    rep("B cells", 70),
    rep("Monocytes", 20)
  )
)

ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = tpm,
  name = "test_dataset"
)

# --- 4. Generate 3 Batches with Different Batch Effects ---
# We will create three separate simulations. The `batch_effect` parameter
# introduces a systematic bias to simulate a batch effect.

# Batch 1: No batch effect
simulation_batch1 <- SimBu::simulate_bulk(
  data = ds,
  scenario = "random",
  ncells = 100,
  nsamples = 30,
  run_parallel = FALSE # Set to FALSE for easier debugging
)
colData(simulation_batch1)$batch <- "Batch 1"

# Batch 2: With a positive batch effect
simulation_batch2 <- SimBu::simulate_bulk(
  data = ds,
  scenario = "random",
  ncells = 100,
  nsamples = 30,
  batch_effect = list(alpha = 0.2, beta = 0.3), # Introduce first batch effect
  run_parallel = FALSE
)
colData(simulation_batch2)$batch <- "Batch 2"

# Batch 3: With a negative batch effect
simulation_batch3 <- SimBu::simulate_bulk(
  data = ds,
  scenario = "random",
  ncells = 100,
  nsamples = 30,
  batch_effect = list(alpha = -0.2, beta = -0.3), # Introduce a different batch effect
  run_parallel = FALSE
)
colData(simulation_batch3)$batch <- "Batch 3"


# Merge the three simulations into a single object
merged_simulations <- SimBu::merge_simulations(list(simulation_batch1, simulation_batch2, simulation_batch3))


# --- 5. PCA Before Batch Correction ---
# Extract data for PCA
counts_before <- assay(merged_simulations)
batch_info <- merged_simulations$batch

# Normalize and log-transform the counts. This is a standard step before PCA/batch correction.
# We use log2(CPM+1) for stabilization.
log_counts_before <- cpm(counts_before, log = TRUE, prior.count = 1)

# Perform PCA
pca_before <- prcomp(t(log_counts_before))
pca_data_before <- data.frame(pca_before$x[, 1:2], batch = batch_info)

# Create the "Before" plot
plot_before <- ggplot(pca_data_before, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Before Batch Correction",
       x = paste0("PC1 (", round(summary(pca_before)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_before)$importance[2,2]*100, 1), "%)")) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  coord_fixed()


# --- 6. Perform Batch Correction ---
# Use limma's removeBatchEffect function on the log-transformed counts
log_counts_corrected <- removeBatchEffect(log_counts_before, batch = batch_info)


# --- 7. PCA After Batch Correction ---
# Perform PCA on the corrected data
pca_after <- prcomp(t(log_counts_corrected))
pca_data_after <- data.frame(pca_after$x[, 1:2], batch = batch_info)

# Create the "After" plot
plot_after <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "After Batch Correction",
       x = paste0("PC1 (", round(summary(pca_after)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_after)$importance[2,2]*100, 1), "%)")) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  coord_fixed()


# --- 8. Display Plots Side-by-Side ---
# Use patchwork to combine the plots
final_plot <- plot_before + plot_after + plot_layout(guides = 'collect')

# Print the final combined plot
print(final_plot)