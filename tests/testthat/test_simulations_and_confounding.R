test_that("generate_simulations returns list with bulk entries", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("SimBu")
  testthat::skip_if_not_installed("gtools")

  n_genes <- 50
  n_cells <- 30
  empirical_means <- rgamma(n_genes, shape = 1.5, rate = 1/0.5)
  empirical_libsizes <- rlnorm(n_cells, meanlog = 0, sdlog = 0.3)
  ds <- create_base_sc_dataset(
    n_genes = n_genes, n_cells = n_cells,
    empirical_means = empirical_means,
    empirical_libsizes = empirical_libsizes,
    dispersion = 0.3
  )

  sim_list <- generate_simulations(ds, n_sims = 2, n_samples_per_sim = 4, n_cells_per_sample = 20)
  expect_type(sim_list, "list")
  expect_equal(length(sim_list), 2)
  expect_true(all(vapply(sim_list, function(x) "bulk" %in% names(x), logical(1))))
})

test_that("create_confound_free_dataset produces bulk_counts assay and metadata", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("SimBu")

  n_genes <- 60
  n_cells <- 40
  base_ds <- create_base_sc_dataset(
    n_genes = n_genes, n_cells = n_cells,
    empirical_means = rgamma(n_genes, 1.5, rate = 2),
    empirical_libsizes = rlnorm(n_cells, 0, 0.4),
    dispersion = 0.4
  )
  sim_list <- generate_simulations(base_ds, n_sims = 3, n_samples_per_sim = 6, n_cells_per_sample = 15)

  prep <- create_confound_free_dataset(
    sim_list,
    n_genes_affected = 10,
    effect_multiplier = 2,
    n_de_genes = 6,
    de_fold_change = 2
  )

  expect_true(all(c("se_with_batch", "batch_info", "biological_vars", "affected_genes", "true_de_genes") %in% names(prep)))
  se <- prep$se_with_batch
  expect_s4_class(se, "SummarizedExperiment")
  expect_true("bulk_counts" %in% SummarizedExperiment::assayNames(se))
  cm <- SummarizedExperiment::assay(se, "bulk_counts")
  expect_true(all(is.finite(cm)))
  expect_equal(length(prep$batch_info), ncol(se))
  expect_equal(length(prep$biological_vars), ncol(se))
})

