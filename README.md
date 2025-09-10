# SimBu - Simulations and Batch Correction Benchmarking for Bulk RNA-seq

## Overview
SimBu provides a reproducible R workflow to:
- simulate bulk RNA-seq from single-cell profiles,
- inject realistic batch and biological effects, and
- benchmark common batch-correction methods with quantitative metrics and plots.

Data type
- Outputs are bulk-level assays in a `SummarizedExperiment` (assay name `bulk_counts`).
- Single-cell data are used only as a reference to parameterize simulations; they are not analyzed directly.

This repo is organized as an analysis project (scripts + helpers). It uses functions from the SimBu toolkit plus wrappers in `functions.r` to generate data, apply corrections, compute metrics, and visualize results.

## What's Included
- `main.r`: End-to-end pipeline (data generation -> batch/biology effects -> corrections -> metrics -> plots). Uses caching for speed (`.Rcache`).
- `functions.r`: Helper functions for simulation, correction methods (limma, ComBat, ComBat-Seq, RUVg/RUVs with k-tuning, fastMNN, PCA-based, SVA/SVA-Seq), and metrics (PERMANOVA R2, silhouettes, kBET, PC regression, iLISI/cLISI, gene-level R2, DE TPR/FPR).
- `PCA_Plots.png`: Example visualization artifact.
- `.github/workflows/R-CMD-check.yaml`: CI scaffold (package setup WIP).

## Quick Start (Project Mode)
Use this mode to run the analysis directly from the scripts.

Prerequisites
- R 4.1+ (tested on recent R)
- On Windows: RTools; on macOS: Xcode CLTs; on Linux: build-essentials

Run
```r
# From the project root in an R session
source("functions.r")
source("main.r")
```
Or from a shell:
```sh
Rscript -e "source('functions.r'); source('main.r')"
```

Outputs
- Printed validation summary and method comparison table (bulk-level)
- PCA and boxplot panels displayed in the plotting device
- Cached intermediates in `.Rcache/` (speeds up re-runs)

## Configuration
Tune the core parameters at the top of `main.r`:
- `NSAMPLES_PER_BATCH`, `N_GENES`, `N_CELLS`, `TOP_GENES`
- Batch effect size: `EFFECT_MULTIPLIER`
- Differential expression: `N_DE_GENES`, `DE_FOLD_CHANGE`
- Simulation count: `NUM_SIMS`
- Reproducibility: `set.seed(123)`

Caching
- The pipeline caches generated data and computed metrics keyed by the parameter set. Change parameters to trigger regeneration; keep them the same to reuse cached results.

## Programmatic Use (Selected Building Blocks)
These helpers are defined in `functions.r`:
```r
# Create single-cell base dataset from empirical means/libsizes
base_ds <- create_base_sc_dataset(
  n_genes = 1000, n_cells = 300,
  empirical_means = rgamma(1000, 1.5, 1/0.5),
  empirical_libsizes = rlnorm(300, 0, 0.5),
  dispersion = 0.4
)

# Focus on top cell-type-specific genes
filtered_ds <- filter_genes_by_specificity(base_ds, top_n = 200)

# Simulate bulk samples across batches
sim_list <- generate_simulations(filtered_ds, n_sims = 4, n_samples_per_sim = 50)

# Inject batch + biology (e.g., Tumor vs Control) and mark ground truth
prep <- create_confound_free_dataset(
  sim_list, n_genes_affected = 100, effect_multiplier = 3.5,
  n_de_genes = 40, de_fold_change = 2.5
)

# Validate simulation realism and independence
val <- validate_simulated_data(
  counts_with_effect = assay(prep$se_with_batch, "bulk_counts"),
  batch_info = prep$batch_info,
  biological_vars = prep$biological_vars,
  affected_genes = prep$affected_genes,
  true_de_genes = prep$true_de_genes
)
```

## Methods Benchmarked
- limma `removeBatchEffect`
- ComBat / ComBat-Seq
- RUVg / RUVs with automatic k tuning (empirical and ideal controls)
- PCA-based regression of batch-associated PCs
- SVA / SVA-Seq
- fastMNN (reduced-dimension integration)

Metrics include: PERMANOVA R2 (batch), silhouette widths (batch/biology), kBET rejection rate, PC regression R2, iLISI/cLISI, gene-level R2 (affected genes), and DE TPR/FPR.

## Installation (Package Mode - optional)
This repository currently ships as a script-based project. If you intend to use it as a package, ensure a standard R package layout (`R/`, `NAMESPACE`, man docs). When the package structure is complete, you can install via:
```r
install.packages("remotes")
remotes::install_github("stef1949/SimBu")
```

## Troubleshooting
- First run installs several CRAN/Bioconductor/GitHub deps (e.g., `welch-lab/kBET`); allow time and compilers.
- If kBET warns about small batches, increase `NSAMPLES_PER_BATCH`.
- On Windows, install RTools before running; on Linux/macOS, ensure build tools are present.

## Contributing
- Open an issue to discuss ideas or bugs
- Fork -> feature branch -> PR
- Keep changes minimal and focused; add concise examples where helpful

## License
MIT - see `LICENSE`.
