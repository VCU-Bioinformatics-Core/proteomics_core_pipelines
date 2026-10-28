# ==========================
# limma analysis functions
# ==========================

#' @title Fit the limma Linear Model Once per Run
#' @details The model is fitted once on all samples; every comparison is then a contrast
#'   on this fit (see \code{\link{perform_limma_analysis}}), so all comparisons share the
#'   same residual variance estimates.
#' @param E Numeric matrix or data frame of log2 intensities (features x samples).
#' @param design Model matrix as produced by \code{model.matrix}.
#' @return An \code{MArrayLM} object from \code{lmFit}.
fit_limma_model <- function(E, design) {
  flog.info("Fitting limma model: %d features x %d samples, coefficients: %s",
            nrow(E), ncol(E), paste(colnames(design), collapse = ", "))
  lmFit(E, design)
}

#' @title Run a limma Differential Expression Analysis for One Contrast
#' @param limma_params Named list containing at minimum:
#'   \describe{
#'     \item{fit}{\code{MArrayLM} object from \code{\link{fit_limma_model}}.}
#'     \item{design}{Model matrix used for the fit.}
#'   }
#' @param contrast Character. limma contrast in terms of the design columns,
#'   e.g. \code{"tumor - control"}.
#' @return Data frame of limma results from \code{topTable} sorted by p-value, with columns
#'   including \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val},
#'   \code{B}. Stops with an error if the analysis fails.
perform_limma_analysis <- function(limma_params, contrast) {
  tryCatch({
    flog.info("Performing limma analysis: %s", contrast)

    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = contrast,
                                     levels = colnames(limma_params$design))

    # Apply contrast to the shared fit
    fit2 <- contrasts.fit(limma_params$fit, contrast_matrix)
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
