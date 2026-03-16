# ==============================================================================
# Custom Imputation Methods
# ==============================================================================
# Define functions here to extend the pipeline with your own imputation strategies.
#
# Convention: name your function  impute_<method_name>
# Then pass  --imputation <method_name>  on the command line.
#
# Interface:
#   Input:  mat    — named numeric matrix (rows = proteins/peptides, cols = samples)
#                    values are log2-transformed intensities; NA = missing value
#           groups — named character vector mapping sample name → group/condition
#                    (names match colnames of mat)
#           ...    — absorb any additional arguments for forward compatibility
#   Output: matrix of identical dimensions and dimnames, with NAs replaced
#
# Example usage on the command line:
#   Rscript de.regular.R --imputation 3by3 ...
#   Rscript de.phospho.R --imputation 3by3 ...
#
# ------------------------------------------------------------------------------
# Example: per-column half-minimum imputation (ignores groups)
# ------------------------------------------------------------------------------
# impute_half_min <- function(mat, groups, ...) {
#   col_half_min <- apply(mat, 2, function(x) min(x, na.rm = TRUE) / 2)
#   for (j in seq_len(ncol(mat))) {
#     mat[is.na(mat[, j]), j] <- col_half_min[j]
#   }
#   mat
# }


# ------------------------------------------------------------------------------
# 3by3: group-aware imputation for experiments with >= 3 samples per group
#
# Per protein, per group:
#   - n_obs >= 2 (1+ missing) → impute all missing as median of the observed values
#   - n_obs == 1 (2+ missing) → global-MinProb (left-tail Gaussian) for missing
#   - n_obs == 0              → global-MinProb for all missing positions
#
# The fallback uses a global distribution (all non-NA values across the full
# matrix) rather than per-column estimates. With only 3 samples per group a
# per-column distribution is too noisy; pooling globally gives stable estimates.
#
# Parameters:
#   q  — quantile used as the mean for the left-tail draw (default 0.01)
# ------------------------------------------------------------------------------
impute_3by3 <- function(mat, groups, q = 0.01, ...) {

  # Align groups to matrix column order by name
  if (!is.null(names(groups))) {
    groups <- groups[colnames(mat)]
  }

  # Validate: every group must have at least 3 samples
  group_sizes <- table(groups)
  small_groups <- group_sizes[group_sizes < 3]
  if (length(small_groups) > 0) {
    stop(paste0(
      "impute_3by3 requires at least 3 samples per group. ",
      "The following group(s) have fewer than 3 samples: ",
      paste(names(small_groups), "(n =", small_groups, ")", collapse = ", "),
      ". Consider using a different imputation method (e.g. --imputation DEP-MinProb)."
    ))
  }

  unique_groups <- unique(groups)

  # Build global left-tail Gaussian from ALL observed values in the matrix.
  # mu    = q-th quantile of the pooled distribution
  # sigma = MAD of the pooled distribution
  all_obs      <- mat[!is.na(mat)]
  global_mu    <- as.numeric(quantile(all_obs, q))
  global_sigma <- mad(all_obs, constant = 1)
  if (is.na(global_sigma) || global_sigma == 0) {
    global_sigma <- sd(all_obs)
  }

  for (i in seq_len(nrow(mat))) {
    for (grp in unique_groups) {
      grp_idx  <- which(groups == grp)
      vals     <- mat[i, grp_idx]
      miss_pos <- grp_idx[is.na(vals)]

      if (length(miss_pos) == 0) next

      n_obs <- sum(!is.na(vals))

      if (n_obs >= 2) {
        # Impute all missing in this group as the median of the observed values
        mat[i, miss_pos] <- median(vals[!is.na(vals)])

      } else {
        # n_obs == 1 or 0: sample from the global left-tail distribution
        mat[i, miss_pos] <- rnorm(length(miss_pos), mean = global_mu, sd = global_sigma)
      }
    }
  }

  mat
}
