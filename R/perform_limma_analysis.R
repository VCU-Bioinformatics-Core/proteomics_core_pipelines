# ==========================
# limma analysis functions
# ==========================

#' @title Build the limma Design Matrix
#' @details Cell-means design: one column per GroupID (\code{~0 + group}), plus
#'   treatment-coded columns for each categorical covariate. Comparisons are contrasts
#'   between the group columns, so covariates adjust every comparison without changing
#'   how contrasts are written. Stops if the design is not full rank, which happens when
#'   a covariate is confounded with the groups.
#' @param samplesheet Data frame with a \code{GroupID} column and any covariate columns.
#' @param covariates Character vector of samplesheet columns, or \code{NULL}. Each is
#'   treated as categorical (numeric codes such as batch 1/2/3 become factor levels).
#' @return Model matrix with group columns named by GroupID. The attribute
#'   \code{"formula_label"} holds a readable formula, e.g. \code{"~0 + GroupID + batch"}.
build_design <- function(samplesheet, covariates = NULL) {
  group <- factor(samplesheet$GroupID)
  if (length(covariates) == 0) {
    design <- model.matrix(~0 + group)
  } else {
    cov_df <- as.data.frame(lapply(samplesheet[covariates], function(x) factor(as.character(x))))
    design <- model.matrix(as.formula(paste("~0 + group +", paste(covariates, collapse = " + "))),
                           data = cbind(data.frame(group = group), cov_df))
  }
  colnames(design)[seq_len(nlevels(group))] <- levels(group)
  colnames(design) <- make.names(colnames(design), unique = TRUE)

  if (!is.fullrank(design))
    stop(sprintf(paste("The design matrix is not full rank; these coefficients cannot be estimated: %s.",
                       "A covariate is confounded with the groups (e.g. every batch holds a single group);",
                       "remove it from `covariates`."),
                 paste(nonEstimable(design), collapse = ", ")))

  attr(design, "formula_label") <- paste(c("~0 + GroupID", covariates), collapse = " + ")
  flog.info("limma design %s: %d samples x %d coefficients (%s)", attr(design, "formula_label"),
            nrow(design), ncol(design), paste(colnames(design), collapse = ", "))
  design
}

#' @title Fit the limma Linear Model Once per Run
#' @details The model is fitted once on all samples; every comparison is then a contrast
#'   on this fit (see \code{\link{perform_limma_analysis}}), so all comparisons share the
#'   same residual variance estimates. With \code{block} (e.g. a patient ID when patients
#'   contribute several samples), the within-block correlation is estimated with
#'   \code{duplicateCorrelation} and used in the fit.
#' @param E Numeric matrix or data frame of log2 intensities (features x samples).
#' @param design Model matrix as produced by \code{\link{build_design}}.
#' @param block Vector (one value per sample, in column order of \code{E}) or \code{NULL}.
#' @return An \code{MArrayLM} object from \code{lmFit}. With \code{block}, its
#'   \code{correlation} element holds the consensus within-block correlation.
fit_limma_model <- function(E, design, block = NULL) {
  flog.info("Fitting limma model: %d features x %d samples, coefficients: %s",
            nrow(E), ncol(E), paste(colnames(design), collapse = ", "))
  if (is.null(block)) return(lmFit(E, design))

  block_sizes <- table(block)
  if (all(block_sizes == 1))
    flog.warn("Every block level has a single sample, so blocking has no effect")
  corfit <- duplicateCorrelation(E, design, block = block)
  rho <- corfit$consensus.correlation
  if (!is.finite(rho)) stop("duplicateCorrelation could not estimate the within-block correlation")
  flog.info("duplicateCorrelation: consensus within-block correlation = %.3f (%d blocks)", rho, length(block_sizes))
  if (rho <= 0)
    flog.warn(paste("Consensus within-block correlation is %.3f (<= 0): samples from the same block are",
                    "not more alike than other samples, so blocking gives little benefit"), rho)
  lmFit(E, design, block = block, correlation = rho)
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

#' @title Contrast Matrix for All Comparisons
#' @param comparisons List of comparison definitions from \code{\link{parse_comparisons}}.
#' @param design Model matrix the contrasts refer to.
#' @return Matrix with one row per design column and one column per comparison.
build_contrast_matrix <- function(comparisons, design) {
  cm <- makeContrasts(contrasts = vapply(comparisons, `[[`, character(1), "contrast"),
                      levels = colnames(design))
  colnames(cm) <- vapply(comparisons, `[[`, character(1), "name")
  cm
}

#' @title Moderated F-test for Any Difference Between Groups
#' @details Tests all group means equal (the limma analogue of a one-way ANOVA) on the
#'   shared fit, so blocking and covariates are taken into account.
#' @param limma_params Named list with \code{fit}, \code{design} and \code{sample_groups}.
#' @return Data frame with columns \code{id} (feature), \code{P_Value} and
#'   \code{Adj_P_Value} (Benjamini-Hochberg), in feature order.
run_group_ftest <- function(limma_params) {
  groups <- unique(as.character(limma_params$sample_groups))
  contrasts <- paste(groups[-1], "-", groups[1])
  contrast_matrix <- makeContrasts(contrasts = contrasts, levels = colnames(limma_params$design))
  fit2 <- eBayes(contrasts.fit(limma_params$fit, contrast_matrix), robust = TRUE)
  res <- topTable(fit2, coef = seq_along(contrasts), number = Inf, sort.by = "none")
  data.frame(id = rownames(res), P_Value = res$P.Value, Adj_P_Value = res$adj.P.Val,
             stringsAsFactors = FALSE)
}
