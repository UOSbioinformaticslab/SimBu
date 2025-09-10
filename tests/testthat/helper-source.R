## Helper to source project functions when not installed as a package

possible_paths <- c(
  "functions.r",
  file.path("..", "..", "functions.r"),
  file.path("..", "functions.r")
)

src <- NULL
for (p in possible_paths) {
  if (file.exists(p)) { src <- p; break }
}

if (!is.null(src)) {
  # Load functions into the calling environment so tests can see them
  sys.source(src, envir = topenv())
} else {
  stop("Could not locate functions.r from tests. Adjust helper-source.R paths.")
}

set.seed(101)

# Avoid parallel warnings on Windows from downstream packages
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  # Force sequential execution (works cross-platform, avoids MulticoreParam warnings)
  BiocParallel::register(BiocParallel::SerialParam())
}
