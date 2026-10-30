# ==========================
# limma analysis functions
# ==========================

#' @title Run a limma Differential Expression Analysis
#' @param limma_params Named list containing at minimum:
#'   \describe{
#'     \item{E}{Numeric matrix of log2 intensities (features × samples).}
#'     \item{design}{Model matrix as produced by \code{model.matrix}.}
#'     \item{block}{Optional vector (one per sample) of block IDs, e.g. patients
#'       contributing several samples. When any block has more than one sample, the
#'       within-block correlation is estimated with \code{duplicateCorrelation} and
#'       used in \code{lmFit}. \code{NULL} fits without blocking.}
#'   }
#' @param exp Character. Column name in \code{limma_params$design} for the experimental group.
#' @param ctrl Character. Column name in \code{limma_params$design} for the control group.
#' @param min_valid_samples Integer. Minimum number of valid (non-missing) samples required
#'   per group. Default \code{2}. (Informational; filtering is expected upstream.)
#' @return Data frame of limma results from \code{topTable} sorted by p-value, with columns
#'   including \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'   \code{B}. When blocking was used, the estimated within-block correlation is
#'   attached as \code{attr(, "block_correlation")} (\code{NULL} otherwise). Stops with
#'   an error if the analysis fails.
perform_limma_analysis <- function(limma_params, exp, ctrl, min_valid_samples=2) {
  tryCatch({
    flog.info("Performing limma analysis: %s vs %s", exp, ctrl)

    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(limma_params$design))

    # Fit linear model; samples sharing a block (e.g. patient) are correlated, so the
    # within-block correlation is estimated and used in the fit
    block <- limma_params$block
    if (!is.null(block) && !any(duplicated(block))) {
      flog.info("No block has more than one sample in this comparison; fitting without blocking")
      block <- NULL
    }
    block_correlation <- NULL
    if (is.null(block)) {
      fit <- lmFit(limma_params$E, limma_params$design)
    } else {
      corfit <- duplicateCorrelation(limma_params$E, limma_params$design, block = block)
      block_correlation <- corfit$consensus.correlation
      flog.info("Within-block correlation (duplicateCorrelation): %.3f", block_correlation)
      fit <- lmFit(limma_params$E, limma_params$design,
                   block = block, correlation = block_correlation)
    }
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2, robust = TRUE)

    # Extract results
    res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    flog.info("Limma analysis complete: %d results", nrow(res))
    res <- as.data.frame(res)
    attr(res, "block_correlation") <- block_correlation
    return(res)
  }, error = function(e) {
    flog.error("Limma analysis failed: %s", e$message)
    stop(e)
  })
}