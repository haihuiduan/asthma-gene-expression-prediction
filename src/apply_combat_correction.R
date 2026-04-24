# Apply ComBat batch effect correction to the merged multi-dataset expression matrix.
# Input expression is samples x genes; ComBat expects genes x samples, so the matrix
# is transposed before correction and transposed back before saving.

expression_file <- "data/processed/merged_expression_common_genes.csv"
labels_file <- "data/processed/merged_labels.csv"
batch_file <- "data/processed/merged_batch_labels.csv"
output_file <- "data/processed/merged_expression_combat_corrected.csv"

# Check required input files before doing any work.
required_files <- c(expression_file, labels_file, batch_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required input file(s):\n",
      paste(missing_files, collapse = "\n"),
      "\nRun the multi-dataset gene alignment notebook before applying ComBat."
    ),
    call. = FALSE
  )
}

# Load sva if available. Provide a clear installation message if it is missing.
if (!requireNamespace("sva", quietly = TRUE)) {
  stop(
    paste0(
      "The Bioconductor package 'sva' is required for ComBat correction.\n",
      "Install it with:\n",
      "install.packages(\"BiocManager\")\n",
      "BiocManager::install(\"sva\")"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(sva)
})

# Read merged expression, labels, and batch labels.
expression <- read.csv(expression_file, row.names = 1, check.names = FALSE)
labels <- read.csv(labels_file, row.names = 1, check.names = FALSE)
batch_labels <- read.csv(batch_file, row.names = 1, check.names = FALSE)

# Validate required columns.
if (!"label" %in% colnames(labels)) {
  stop("merged_labels.csv must contain a column named 'label'.", call. = FALSE)
}

if (!"dataset" %in% colnames(batch_labels)) {
  stop("merged_batch_labels.csv must contain a column named 'dataset'.", call. = FALSE)
}

# Align all tables by shared sample ID.
common_samples <- Reduce(
  intersect,
  list(rownames(expression), rownames(labels), rownames(batch_labels))
)

if (length(common_samples) == 0) {
  stop("No overlapping sample IDs found across expression, labels, and batch labels.", call. = FALSE)
}

expression <- expression[common_samples, , drop = FALSE]
labels <- labels[common_samples, , drop = FALSE]
batch_labels <- batch_labels[common_samples, , drop = FALSE]

# Convert expression values to a numeric matrix.
expression_matrix <- as.matrix(expression)
mode(expression_matrix) <- "numeric"

if (anyNA(expression_matrix)) {
  stop("Expression matrix contains missing or non-numeric values after conversion.", call. = FALSE)
}

label <- factor(labels$label)
batch <- factor(batch_labels$dataset)

# Print dimensions and counts before ComBat correction.
cat("Aligned expression matrix dimensions before ComBat (samples x genes):\n")
print(dim(expression_matrix))

cat("\nBatch counts:\n")
print(table(batch))

cat("\nLabel counts:\n")
print(table(label))

# Preserve the biological label effect while correcting dataset/batch effects.
mod <- model.matrix(~ label)

# ComBat expects genes x samples.
combat_input <- t(expression_matrix)

cat("\nComBat input dimensions (genes x samples):\n")
print(dim(combat_input))

combat_corrected <- sva::ComBat(
  dat = combat_input,
  batch = batch,
  mod = mod,
  par.prior = TRUE
)

# Return to samples x genes for downstream Python notebooks.
corrected_expression <- t(combat_corrected)

cat("\nCorrected expression matrix dimensions (samples x genes):\n")
print(dim(corrected_expression))

# Save corrected expression with sample IDs as rows and gene symbols as columns.
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(corrected_expression, output_file, row.names = TRUE)

message("Saved ComBat-corrected expression matrix to: ", output_file)
