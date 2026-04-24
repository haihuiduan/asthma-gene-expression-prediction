# Extract expression, metadata, and feature annotation from local GEO series matrix files.
# This script does not download data. It only reads files that already exist in data/raw/.

raw_dir <- "data/raw"
processed_dir <- "data/processed"

# Define the local datasets to extract.
datasets <- list(
  GSE40732 = file.path(raw_dir, "GSE40732_series_matrix.txt.gz"),
  GSE40888 = file.path(raw_dir, "GSE40888_series_matrix.txt.gz"),
  GSE64913 = file.path(raw_dir, "GSE64913_series_matrix.txt.gz")
)

# Load required packages.
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# Ensure the output directory exists before writing CSV files.
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Extract one local GEO series matrix file and save dataset-specific outputs.
extract_series_matrix <- function(dataset_id, raw_file) {
  dataset_prefix <- tolower(dataset_id)

  expression_out <- file.path(processed_dir, paste0(dataset_prefix, "_expression.csv"))
  metadata_out <- file.path(processed_dir, paste0(dataset_prefix, "_metadata.csv"))
  features_out <- file.path(processed_dir, paste0(dataset_prefix, "_features.csv"))

  # Skip missing local files so one absent dataset does not stop the full batch.
  if (!file.exists(raw_file)) {
    warning(
      paste0(
        "Skipping ", dataset_id, ": missing raw file at ", raw_file,
        ". Place the GEO series matrix file in data/raw/ to extract it."
      ),
      call. = FALSE
    )
    return(invisible(FALSE))
  }

  message("Loading ", dataset_id, " from: ", raw_file)

  # Load the local GEO series matrix file.
  gse <- GEOquery::getGEO(filename = raw_file)

  # Extract expression values, sample metadata, and feature annotation.
  expression_matrix <- Biobase::exprs(gse)
  metadata <- Biobase::pData(gse)
  features <- Biobase::fData(gse)

  # Save all extracted tables as CSV files for downstream Python notebooks.
  write.csv(expression_matrix, expression_out, row.names = TRUE)
  write.csv(metadata, metadata_out, row.names = TRUE)
  write.csv(features, features_out, row.names = TRUE)

  message("Saved expression matrix to: ", expression_out)
  message("Saved metadata table to: ", metadata_out)
  message("Saved feature annotation to: ", features_out)

  invisible(TRUE)
}

# Process each configured dataset.
for (dataset_id in names(datasets)) {
  extract_series_matrix(dataset_id, datasets[[dataset_id]])
}

message("GEO series matrix extraction complete.")
