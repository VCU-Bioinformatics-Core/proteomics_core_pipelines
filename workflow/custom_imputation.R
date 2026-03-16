# ==============================================================================
# Custom Imputation Methods
# ==============================================================================
# Define functions here to extend the pipeline with your own imputation strategies.
#
# Convention: name your function  impute_<method_name>
# Then pass  --imputation <method_name>  on the command line.
#
# Interface:
#   Input:  mat  — named numeric matrix (rows = proteins/peptides, cols = samples)
#                  values are log2-transformed intensities; NA = missing value
#   Output: matrix of identical dimensions and dimnames, with NAs replaced
#
# Example usage on the command line:
#   Rscript de.regular.R --imputation half_min ...
#   Rscript de.phospho.R --imputation half_min ...
#
# ------------------------------------------------------------------------------
# Example: per-column half-minimum imputation
# Replace each NA with half the minimum observed value in that sample column.
# ------------------------------------------------------------------------------
# impute_half_min <- function(mat, ...) {
#   col_half_min <- apply(mat, 2, function(x) min(x, na.rm = TRUE) / 2)
#   for (j in seq_len(ncol(mat))) {
#     mat[is.na(mat[, j]), j] <- col_half_min[j]
#   }
#   mat
# }

# Add your custom imputation functions below:
