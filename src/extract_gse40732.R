# Extract GSE40732 expression and metadata tables for Python notebooks.
# This script uses a local GEO series matrix file and does not download data.

raw_file <- "data/raw/GSE40732_series_matrix.txt.gz"
expression_out <- "data/processed/gse40732_expression.csv"
metadata_out <- "data/processed/gse40732_metadata.csv"

# Stop early with a clear message if the expected local file is missing.
if (!file.exists(raw_file)) {
  stop(
    paste0(
      "Missing raw GEO series matrix file: ", raw_file, "\n",
      "Place GSE40732_series_matrix.txt.gz in data/raw/ before running this script."
    ),
    call. = FALSE
  )
}

# Load required packages.
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# Load the local GEO series matrix file.
gse <- GEOquery::getGEO(filename = raw_file)

# Extract the expression matrix with genes/probes as rows and samples as columns.
expression_matrix <- Biobase::exprs(gse)

# Extract sample-level metadata from the ExpressionSet.
metadata <- Biobase::pData(gse)

# Ensure the processed data directory exists before writing outputs.
dir.create(dirname(expression_out), recursive = TRUE, showWarnings = FALSE)

# Save expression and metadata tables as CSV files for downstream Python notebooks.
write.csv(expression_matrix, expression_out, row.names = TRUE)
write.csv(metadata, metadata_out, row.names = TRUE)

message("Saved expression matrix to: ", expression_out)
message("Saved metadata table to: ", metadata_out)
