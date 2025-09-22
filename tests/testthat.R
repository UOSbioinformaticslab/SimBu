if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install 'testthat' to run the tests.")
}

# Attach testthat so test_that/expect_* are available unqualified
library(testthat)

# Use a direct runner so this works even without full package structure
testthat::test_dir("tests/testthat", reporter = "summary")
