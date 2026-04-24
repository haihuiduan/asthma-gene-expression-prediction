# Map GSE40732 RefSeq/GenBank accessions to standard human gene symbols.
# This script reads the processed feature annotation table and adds a GeneSymbol column.

input_file <- "data/processed/gse40732_features.csv"
output_file <- "data/processed/gse40732_features_with_symbols.csv"

# Stop early if the expected processed feature table is missing.
if (!file.exists(input_file)) {
  stop(
    paste0(
      "Missing input file: ", input_file, "\n",
      "Run the GEO extraction script before mapping GSE40732 accessions."
    ),
    call. = FALSE
  )
}

# Load Bioconductor annotation packages.
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Read the processed feature table.
features <- read.csv(input_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Prefer GB_ACC, but keep a fallback for similarly named accession columns.
accession_candidates <- c("GB_ACC", "REFSEQ", "RefSeq", "refseq", "ID")
accession_column <- accession_candidates[accession_candidates %in% colnames(features)][1]

if (is.na(accession_column)) {
  stop(
    paste0(
      "No accession column found in ", input_file, ". Expected one of: ",
      paste(accession_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

message("Using accession column: ", accession_column)

# Clean accession IDs by removing version suffixes such as NM_001234.5.
raw_accessions <- trimws(as.character(features[[accession_column]]))
clean_accessions <- sub("\\.[0-9]+$", "", raw_accessions)
clean_accessions[clean_accessions == ""] <- NA

# Map RefSeq accessions to human gene symbols.
gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(na.omit(clean_accessions)),
  column = "SYMBOL",
  keytype = "REFSEQ",
  multiVals = "first"
)

# Add mapped symbols back to the feature table.
features$CleanRefSeq <- clean_accessions
features$GeneSymbol <- unname(gene_symbols[clean_accessions])

mapped_count <- sum(!is.na(features$GeneSymbol) & features$GeneSymbol != "")
unmapped_count <- nrow(features) - mapped_count

message("Rows mapped successfully: ", mapped_count)
message("Rows unmapped: ", unmapped_count)

# Save the updated feature table for downstream Python notebooks.
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(features, output_file, row.names = TRUE)

message("Saved updated feature table to: ", output_file)
