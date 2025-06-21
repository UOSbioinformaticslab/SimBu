# --- 1. Installation of Required Packages ---
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_to_install <- c("SimBu", "limma", "edgeR", "ggplot2", "patchwork")
for (pkg in packages_to_install) {
  if (!require(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- 2. Load Libraries ---
library(SimBu)
library(limma)
library(edgeR)
library(ggplot2)
library(patchwork)

# --- 3. Create Base Single-Cell Dataset ---
gene_names <- paste0("gene_", 1:1000)
cell_names <- paste0("cell_", 1:300)

counts <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
rownames(counts) <- gene_names
colnames(counts) <- cell_names

tpm <- Matrix::Matrix(matrix(stats::rpois(3e5, 5), ncol = 300), sparse = TRUE)
tpm <- Matrix::t(1e6 * Matrix::t(tpm) / Matrix::colSums(tpm))
rownames(tpm) <- gene_names
colnames(tpm) <- cell_names

annotation <- data.frame(
  "ID" = cell_names,
  "cell_type" = c(
    rep("T cells CD4", 50),
    rep("T cells CD8", 50),
    rep("Macrophages", 100),
    rep("NK cells", 10),
    rep("B cells", 70),
    rep("Monocytes", 20)
  )
)

# ds is an S4 object of class "SimBu_dataset". We must use @ to access its slots.
ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = tpm,
  name = "test_dataset"
)

# --- 4. Generate 3 Simulations ---
set.seed(123) 
nsamples_per_batch <- 30

simulation_batch1 <- SimBu::simulate_bulk(
  data = ds, scenario = "random", scaling_factor = "NONE",
  ncells = 100, nsamples = nsamples_per_batch, run_parallel = FALSE
)
simulation_batch2 <- SimBu::simulate_bulk(
  data = ds, scenario = "random", scaling_factor = "NONE",
  ncells = 100, nsamples = nsamples_per_batch, run_parallel = FALSE
)
simulation_batch3 <- SimBu::simulate_bulk(
  data = ds, scenario = "random", scaling_factor = "NONE",
  ncells = 100, nsamples = nsamples_per_batch, run_parallel = FALSE
)

# --- 5. Merge Simulations and Manually Add Batch Effect ---
# merge_simulations returns a standard list, so we use $
merged_list <- SimBu::merge_simulations(list(simulation_batch1, simulation_batch2, simulation_batch3))

batch_info <- rep(c("Batch 1", "Batch 2", "Batch 3"), each = nsamples_per_batch)
counts_before <- assay(merged_list$bulk, "bulk_counts")

# *** THE DEFINITIVE FIX: Use the '@' operator to access the S4 object slot ***
rownames(counts_before) <- rownames(ds)

# Now, rownames(counts_before) is a valid character vector, and sample() will work.
set.seed(456)
n_genes_affected <- 150
affected_genes <- sample(rownames(counts_before), n_genes_affected)

batch2_cols <- which(batch_info == "Batch 2")
batch3_cols <- which(batch_info == "Batch 3")

counts_before[affected_genes, batch2_cols] <- round(counts_before[affected_genes, batch2_cols] * 2.5)
counts_before[affected_genes, batch3_cols] <- round(counts_before[affected_genes, batch3_cols] / 2.5)


# --- 6. PCA Before Batch Correction ---
log_counts_before <- cpm(counts_before, log = TRUE, prior.count = 1)
pca_before <- prcomp(t(log_counts_before))
pca_data_before <- data.frame(pca_before$x[, 1:2], batch = batch_info)

plot_before <- ggplot(pca_data_before, aes(PC1, PC2, color = batch)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(
    title = "Pre Batch Correction",
    x = paste0("PC1 (", round(summary(pca_before)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_before)$importance[2,2]*100, 1), "%)")
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  coord_fixed()

plot_after <- ggplot(pca_data_after, aes(PC1, PC2, color = batch)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(
    title = "Post removeBatchEffect Batch Correction",
    x = paste0("PC1 (", round(summary(pca_after)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_after)$importance[2,2]*100, 1), "%)")
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  coord_fixed()

final_plot <- (plot_before + plot_after + plot_layout(guides = 'collect')) &
  theme(aspect.ratio = 1)

print(final_plot)


