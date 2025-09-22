if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install 'testthat' to run the tests.")
}

library(testthat)

# Run all tests in tests/testthat
testthat::test_dir("tests/testthat", reporter = "summary")
