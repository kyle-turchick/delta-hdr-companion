# Compares file fingerprints and CSV contents without writing when sourced.

sha256_file <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot calculate a SHA-256 fingerprint because the file is missing: ", path)
  }
  if (!requireNamespace("openssl", quietly = TRUE)) {
    stop("R package 'openssl' is required to calculate SHA-256 file fingerprints.")
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  as.character(openssl::sha256(con))
}

compare_csv_files <- function(current, baseline) {
  current_tbl <- readr::read_csv(current, show_col_types = FALSE)
  baseline_tbl <- readr::read_csv(baseline, show_col_types = FALSE)
  isTRUE(all.equal(current_tbl, baseline_tbl, check.attributes = FALSE))
}

file_mtime_or_na <- function(path) {
  if (file.exists(path)) file.info(path)$mtime else as.POSIXct(NA)
}
