test_that("create_base_sc_dataset returns a valid SummarizedExperiment with counts", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("Matrix")
  testthat::skip_if_not_installed("SimBu")

  n_genes <- 60
  n_cells <- 40
  empirical_means <- rgamma(n_genes, shape = 1.5, rate = 1/0.5)
  empirical_libsizes <- rlnorm(n_cells, meanlog = 0, sdlog = 0.3)

  ds <- create_base_sc_dataset(
    n_genes = n_genes, n_cells = n_cells,
    empirical_means = empirical_means,
    empirical_libsizes = empirical_libsizes,
    dispersion = 0.3
  )

  expect_s4_class(ds, "SummarizedExperiment")
  expect_true("counts" %in% SummarizedExperiment::assayNames(ds))
  # Implementations may drop all-zero genes; allow <= n_genes
  expect_lte(nrow(ds), n_genes)
  expect_gt(nrow(ds), 0)
  expect_equal(ncol(ds), n_cells)
  cm <- SummarizedExperiment::assay(ds, "counts")
  expect_true(all(cm >= 0))
})

test_that("filter_genes_by_specificity reduces to top_n and preserves metadata", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("SimBu")

  n_genes <- 80
  n_cells <- 30
  empirical_means <- rgamma(n_genes, shape = 1.2, rate = 1/0.6)
  empirical_libsizes <- rlnorm(n_cells, meanlog = 0, sdlog = 0.2)

  ds <- create_base_sc_dataset(
    n_genes = n_genes, n_cells = n_cells,
    empirical_means = empirical_means,
    empirical_libsizes = empirical_libsizes,
    dispersion = 0.4
  )

  top_n <- 25
  filt <- filter_genes_by_specificity(ds, top_n = top_n)
  expect_s4_class(filt, "SummarizedExperiment")
  expect_lte(nrow(filt), top_n)
  expect_true("cell_type" %in% colnames(SummarizedExperiment::colData(filt)))
})
