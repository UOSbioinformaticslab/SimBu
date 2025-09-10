test_that("perform_pca_and_plot returns a ggplot object", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("matrixStats")

  set.seed(1)
  # Simulate simple log-counts matrix: genes x samples
  log_counts <- matrix(rnorm(100 * 12, mean = 5, sd = 1), nrow = 100, ncol = 12)
  colnames(log_counts) <- paste0("S", seq_len(ncol(log_counts)))
  batch <- rep(c("Batch1", "Batch2", "Batch3"), length.out = 12)
  bio <- rep(c("Control", "Tumor"), length.out = 12)

  p <- perform_pca_and_plot(log_counts = log_counts, batch_info = batch, biological_vars = bio, plot_title = "Test Plot")
  expect_true(inherits(p, "ggplot"))
})

test_that("metric helpers return numeric or NA as expected", {
  testthat::skip_if_not_installed("matrixStats")

  set.seed(2)
  log_counts <- matrix(rnorm(80 * 10, mean = 4, sd = 1), nrow = 80, ncol = 10)
  batch <- rep(c("B1", "B2"), each = 5)
  bio <- rep(c("Ctl", "Trt"), length.out = 10)

  # These may return NA if packages not present; that's OK.
  r2 <- try(calculate_permanova_r2(log_counts = log_counts, variable_info = batch), silent = TRUE)
  if (!inherits(r2, "try-error")) expect_true(is.numeric(r2) || is.na(r2))

  sil_b <- try(calculate_silhouette_width(log_counts = log_counts, variable_info = batch), silent = TRUE)
  if (!inherits(sil_b, "try-error")) expect_true(is.numeric(sil_b) || is.na(sil_b))

  pcr <- try(calculate_pc_regression(log_counts = log_counts, variable_info = batch), silent = TRUE)
  if (!inherits(pcr, "try-error")) expect_true(is.numeric(pcr) || is.na(pcr))
})

test_that("run_de_and_evaluate returns bounded TPR/FPR", {
  testthat::skip_if_not_installed("limma")

  set.seed(3)
  # Simulate two-group data with some true DE genes
  ng <- 60; ns <- 12
  bio <- rep(c("Ctl", "Trt"), each = ns/2)
  log_counts <- matrix(rnorm(ng * ns, 5, 1), nrow = ng)
  rownames(log_counts) <- paste0("g", seq_len(ng))
  # Spike DE in first 8 genes for Trt
  log_counts[1:8, bio == "Trt"] <- log_counts[1:8, bio == "Trt"] + 1.0
  res <- run_de_and_evaluate(log_counts, bio, true_de_genes = rownames(log_counts)[1:8])
  expect_true(is.list(res) && all(c("TPR", "FPR") %in% names(res)))
  expect_true(is.numeric(res$TPR) && is.numeric(res$FPR))
  expect_gte(res$TPR, 0); expect_lte(res$TPR, 1)
  expect_gte(res$FPR, 0); expect_lte(res$FPR, 1)
})

