# ==========================
# limma analysis functions
# ==========================

#' @title Run a limma Differential Expression Analysis
#' @param limma_params Named list containing at minimum:
#'   \describe{
#'     \item{E}{Numeric matrix of log2 intensities (features × samples).}
#'     \item{design}{Model matrix as produced by \code{model.matrix}.}
#'   }
#' @param exp Character. Column name in \code{limma_params$design} for the experimental group.
#' @param ctrl Character. Column name in \code{limma_params$design} for the control group.
#' @param min_valid_samples Integer. Minimum number of valid (non-missing) samples required
#'   per group. Default \code{2}. (Informational; filtering is expected upstream.)
#' @return Data frame of limma results from \code{topTable} sorted by p-value, with columns
#'   including \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'   \code{B}. Stops with an error if the analysis fails.
perform_limma_analysis <- function(limma_params, exp, ctrl, min_valid_samples=2) {
  tryCatch({
    flog.info("Performing limma analysis: %s vs %s", exp, ctrl)

    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(limma_params$design))

    # Fit linear model
    fit <- lmFit(limma_params$E, limma_params$design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2, robust = TRUE)

    # Extract results
    res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    flog.info("Limma analysis complete: %d results", nrow(res))
    return(as.data.frame(res))
  }, error = function(e) {
    flog.error("Limma analysis failed: %s", e$message)
    stop(e)
  })
}