# ==========================
# Report section: the limma model in equations
# ==========================
# describe_limma_model() — Markdown + LaTeX describing the model, its terms, the design
#                          matrix and the contrasts, for the HTML reports

#' @title LaTeX-safe Group or Coefficient Name
#' @param x Character vector.
#' @return Character vector of \code{\\text{...}} snippets with underscores escaped.
latex_text <- function(x) paste0("\\text{", gsub("_", "\\\\_", x), "}")

#' @title Write a limma Contrast Formula in LaTeX
#' @details Keeps the formula as the analyst wrote it: group names become
#'   \code{\\mu_{g,group}}, \code{(...)/n} becomes \code{\\tfrac{1}{n}(...)} and
#'   \code{*} becomes \code{\\cdot}.
#' @param contrast Character. Contrast formula, e.g. \code{"(V_R + L_R)/2 - (V_NR + L_NR)/2"}.
#' @param groups Character vector of group names.
#' @return Character string of LaTeX.
formula_to_latex <- function(contrast, groups) {
  ids <- gregexpr("[A-Za-z][A-Za-z0-9._]*|\\.[A-Za-z_][A-Za-z0-9._]*", contrast)
  regmatches(contrast, ids) <- lapply(regmatches(contrast, ids), function(tok) {
    ifelse(tok %in% groups, paste0("\\mu_{g,", latex_text(tok), "}"), latex_text(tok))
  })
  contrast <- gsub("\\(([^()]*)\\)\\s*/\\s*([0-9]+)", "\\\\tfrac{1}{\\2}\\\\left(\\1\\\\right)", contrast)
  gsub("\\*", " \\\\cdot ", contrast)
}

#' @title Markdown Table from a Data Frame
#' @param df Data frame.
#' @return Character string with a pipe-table.
md_table <- function(df) {
  paste(knitr::kable(df, format = "pipe", row.names = FALSE), collapse = "\n")
}

#' @title Describe the limma Model for the Report
#' @details Builds the "Statistical Model" report section: the model equation and its
#'   matrix form, a definition of every term, the sample correlation structure (when
#'   \code{block} is used), the design matrix, each comparison as a contrast, the
#'   moderated t-statistic, and the equivalent limma code. Output is Markdown with LaTeX
#'   math (rendered by MathJax in the HTML report), meant for a \code{results="asis"}
#'   chunk. Works with results saved before the design was stored, by rebuilding the
#'   cell-means design from \code{sample_info}.
#' @param analysis_params Named list saved by the pipeline (uses \code{design},
#'   \code{contrast_matrix}, \code{design_samplesheet}, \code{design_formula},
#'   \code{group_factors}, \code{covariates}, \code{block}, \code{block_correlation}).
#' @param sample_info Data frame with columns \code{sample} and \code{condition}.
#' @param comparisons List of comparison definitions.
#' @param feature Character. Singular feature name, e.g. \code{"protein"} or
#'   \code{"PTM peptide"}.
#' @param n_features Integer or \code{NULL}. Number of features tested.
#' @param imputation_method Character or \code{NULL}. Imputation method used.
#' @return Character string of Markdown.
describe_limma_model <- function(analysis_params, sample_info, comparisons, feature = "protein",
                                 n_features = NULL, imputation_method = NULL) {
  groups <- unique(as.character(sample_info$condition))
  design <- analysis_params$design
  if (is.null(design)) {
    design <- model.matrix(~0 + factor(sample_info$condition))
    colnames(design) <- levels(factor(sample_info$condition))
  }
  covariate_cols <- setdiff(colnames(design), groups)
  covariates     <- analysis_params$covariates
  block          <- analysis_params$block
  rho            <- analysis_params$block_correlation
  sheet          <- analysis_params$design_samplesheet
  n <- nrow(design); p <- ncol(design); K <- length(groups); L <- length(covariate_cols)
  formula_label  <- analysis_params$design_formula %||% "~0 + GroupID"

  cm <- analysis_params$contrast_matrix
  if (is.null(cm)) {
    cm <- makeContrasts(contrasts = vapply(comparisons, function(x) x$contrast %||% paste(x$exp, "-", x$ctrl), ""),
                        levels = colnames(design))
    colnames(cm) <- vapply(comparisons, `[[`, "", "name")
  }

  out <- character()
  add <- function(...) out <<- c(out, paste0(...))

  # ---- introduction ----
  add("Differential abundance was tested with **limma** (linear models with empirical Bayes moderation). ",
      "One linear model is fitted for every ", feature, " using **all ", n, " samples** at once; each comparison ",
      "is then a *contrast* of that single fit. Because all comparisons share one model, every ", feature,
      "'s variance is estimated from all samples, which gives more power than analyzing pairs of groups separately.\n")

  # ---- model equation ----
  add("### Model\n")
  add("For ", feature, " $g$ and sample $i$:\n")
  covariate_term <- if (L > 0) " + \\sum_{l=1}^{L} \\gamma_{gl}\\, z_{il}" else ""
  add("$$\ny_{gi} = \\mu_{g,G(i)}", covariate_term, " + \\varepsilon_{gi}\n$$\n")
  add("or, for all samples at once (matrix form):\n")
  add("$$\n\\mathbf{y}_g = \\mathbf{X}\\boldsymbol{\\beta}_g + \\boldsymbol{\\varepsilon}_g,",
      " \\qquad \\boldsymbol{\\varepsilon}_g \\sim \\mathcal{N}\\!\\left(\\mathbf{0},\\ \\sigma_g^2\\,\\mathbf{V}\\right)\n$$\n")

  # ---- term definitions ----
  group_source <- if (!is.null(analysis_params$group_factors)) {
    sprintf("built from the samplesheet columns %s (joined with `_`)",
            paste0("`", analysis_params$group_factors, "`", collapse = " × "))
  } else {
    "the samplesheet `GroupID` column"
  }
  y_desc <- paste0("log2(intensity + 1)",
                   if (!is.null(imputation_method) && imputation_method != "none")
                     paste0(", after imputation (`", imputation_method, "`)") else ", no imputation")
  terms <- data.frame(
    Term = c("$y_{gi}$", "$g$", "$i$", "$G(i)$", "$\\mu_{g,k}$"),
    Meaning = c(paste0("Abundance of ", feature, " $g$ in sample $i$"),
                paste0("Index of the ", feature),
                "Index of the sample",
                "The group of sample $i$",
                paste0("Mean log2 abundance of ", feature, " $g$ in group $k$ (estimated)")),
    `In this analysis` = c(y_desc,
                           if (!is.null(n_features)) paste0(n_features, " ", feature, "s tested") else paste0("each ", feature),
                           paste0(n, " samples (`SampleID`)"),
                           paste0(K, " groups, ", group_source, ": ", paste0("`", groups, "`", collapse = ", ")),
                           paste0("one per group: ", paste0("$\\mu_{g,", latex_text(groups), "}$", collapse = ", "))),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (L > 0) {
    cov_desc <- vapply(covariates, function(cv) {
      lv <- if (!is.null(sheet) && cv %in% colnames(sheet)) levels(factor(as.character(sheet[[cv]]))) else character()
      if (length(lv) == 0) return(paste0("`", cv, "`"))
      sprintf("`%s` (levels %s; reference level `%s` is absorbed into the group means)",
              cv, paste0("`", lv, "`", collapse = ", "), lv[1])
    }, character(1))
    terms <- rbind(terms, data.frame(
      Term = c("$z_{il}$", "$\\gamma_{gl}$"),
      Meaning = c("Covariate indicator: 1 if sample $i$ is in covariate level $l$, otherwise 0",
                  paste0("Shift in log2 abundance of ", feature, " $g$ for covariate level $l$ relative to the reference level (estimated, adjusted for)")),
      `In this analysis` = c(paste0(L, " column(s): ", paste0("`", covariate_cols, "`", collapse = ", "),
                                    ", from ", paste(cov_desc, collapse = "; ")),
                             "nuisance effects; they are not compared, only removed"),
      check.names = FALSE, stringsAsFactors = FALSE))
  }
  terms <- rbind(terms, data.frame(
    Term = c("$\\varepsilon_{gi}$", "$\\sigma_g^2$", "$\\mathbf{X}$", "$\\boldsymbol{\\beta}_g$", "$\\mathbf{V}$"),
    Meaning = c("Residual: the part of $y_{gi}$ the model does not explain",
                paste0("Residual variance of ", feature, " $g$"),
                "Design matrix: one row per sample, one column per coefficient (shown below)",
                paste0("All coefficients of ", feature, " $g$: the group means", if (L > 0) " and covariate effects" else ""),
                "Correlation between samples"),
    `In this analysis` = c("assumed normal with mean 0",
                           "moderated across all features by empirical Bayes (see *Testing* below)",
                           sprintf("%d × %d, formula `%s`", n, p, formula_label),
                           sprintf("%d coefficients: %s", p, paste0("`", colnames(design), "`", collapse = ", ")),
                           if (!is.null(block)) sprintf("samples sharing `%s` are correlated (see below)", block)
                           else "identity: samples are independent"),
    check.names = FALSE, stringsAsFactors = FALSE))
  add("**Terms**\n")
  add(md_table(terms), "\n")

  # ---- correlation structure ----
  if (!is.null(block)) {
    add("### Samples from the same `", block, "`\n")
    block_info <- ""
    if (!is.null(sheet) && block %in% colnames(sheet)) {
      sizes <- table(sheet[[block]])
      block_info <- sprintf(" There are %d `%s` levels with %s samples each.", length(sizes), block,
                            paste(sort(unique(as.integer(sizes))), collapse = " or "))
    }
    add("Several samples come from the same `", block, "`, so they are not independent.", block_info,
        " limma's `duplicateCorrelation` estimates one consensus correlation $\\rho$ shared by all ", feature, "s, ",
        "and the model uses it in $\\mathbf{V}$:\n")
    add("$$\nV_{ij} = \\begin{cases} 1 & i = j \\\\ \\rho & i \\neq j \\text{ and both samples share the same ",
        gsub("_", "\\\\_", block), "} \\\\ 0 & \\text{otherwise} \\end{cases}",
        if (!is.null(rho)) sprintf(" \\qquad \\hat{\\rho} = %.3f", rho) else "", "\n$$\n")
    add("A $\\hat\\rho$ near 0 means samples from the same `", block, "` are not more alike than other samples ",
        "(blocking changes little); a larger $\\hat\\rho$ means ignoring the blocking would have made p-values too small ",
        "for comparisons across ", block, " levels.\n")
  }

  # ---- design matrix ----
  add("### Design matrix\n")
  add("limma uses a **cell-means** design (`", formula_label, "`): there is no intercept, and each group has its own column. ",
      "A `1` in a group column means that sample's expected value is that group's mean $\\mu_{g,k}$",
      if (L > 0) "; covariate columns add the shift $\\gamma_{gl}$ for samples in that covariate level" else "",
      ". Rows are samples; columns are the coefficients $\\boldsymbol{\\beta}_g$.\n")
  design_df <- data.frame(SampleID = as.character(sample_info$sample), Group = as.character(sample_info$condition),
                          stringsAsFactors = FALSE)
  if (!is.null(block) && !is.null(sheet) && block %in% colnames(sheet)) design_df[[block]] <- as.character(sheet[[block]])
  design_df <- cbind(design_df, as.data.frame(matrix(design, nrow = n, dimnames = list(NULL, colnames(design))),
                                              check.names = FALSE))
  if (n > 24) {
    add("<details><summary>Show the design matrix (", n, " samples)</summary>\n")
    add(md_table(design_df), "\n")
    add("</details>\n")
  } else {
    add(md_table(design_df), "\n")
  }

  # ---- contrasts ----
  add("### Comparisons as contrasts\n")
  add("Each comparison is a weighted combination of the coefficients, $\\theta_g = \\mathbf{c}^\\top\\boldsymbol{\\beta}_g$. ",
      "Its estimate $\\hat\\theta_g$ is the **log2 fold change (logFC)** reported in the results: positive values mean higher ",
      "abundance on the `+` side.\n")
  for (cmp in comparisons) {
    add("- **", cmp$name, "**: $\\theta_g = ",
        formula_to_latex(cmp$contrast %||% paste(cmp$exp, "-", cmp$ctrl), groups), "$")
  }
  add("")
  cm_df <- data.frame(Coefficient = rownames(cm),
                      as.data.frame(matrix(cm, nrow = nrow(cm), dimnames = dimnames(cm)), check.names = FALSE),
                      check.names = FALSE, stringsAsFactors = FALSE)
  cm_df[-1] <- lapply(cm_df[-1], function(x) format(round(x, 3), drop0trailing = TRUE))
  add("\nContrast matrix (one column $\\mathbf{c}$ per comparison):\n")
  add(md_table(cm_df), "\n")

  # ---- testing ----
  add("### Testing\n")
  add("For each comparison, limma computes a **moderated t-statistic**:\n")
  add("$$\n\\tilde{t}_g = \\frac{\\hat{\\theta}_g}{\\tilde{s}_g\\,\\sqrt{\\mathbf{c}^\\top \\left(\\mathbf{X}^\\top \\mathbf{V}^{-1} \\mathbf{X}\\right)^{-1} \\mathbf{c}}},",
      " \\qquad \\tilde{s}_g^{2} = \\frac{d_0\\, s_0^{2} + d_g\\, s_g^{2}}{d_0 + d_g}\n$$\n")
  testing <- data.frame(
    Term = c("$\\hat{\\theta}_g$", "$s_g^2$", "$d_g$", "$s_0^2,\\ d_0$", "$\\tilde{s}_g^2$", "$\\tilde{t}_g$"),
    Meaning = c("Estimated contrast = logFC",
                paste0("Residual variance of ", feature, " $g$ from its own fit"),
                "Residual degrees of freedom",
                paste0("Prior variance and prior degrees of freedom, estimated from all ", feature, "s (robust empirical Bayes, `eBayes(robust = TRUE)`)"),
                paste0("Moderated variance: $s_g^2$ shrunk toward $s_0^2$, which stabilizes ", feature, "s with few samples or extreme variances"),
                "Moderated t-statistic (`t` column), compared with a t distribution with $d_0 + d_g$ degrees of freedom"),
    `In this analysis` = c("`logFC` column",
                           paste0("estimated separately for each ", feature),
                           sprintf("$n - p = %d - %d = %d$ (fewer for %ss with missing values)", n, p, n - p, feature),
                           paste0("estimated once, from all ", feature, "s together"),
                           "used in place of $s_g^2$",
                           "`t` and `P.Value` columns"),
    check.names = FALSE, stringsAsFactors = FALSE)
  add(md_table(testing), "\n")
  add("P-values are adjusted across ", feature, "s with the Benjamini–Hochberg method (`adj.P.Val`). ",
      "A ", feature, " is called **significant** when `adj.P.Val` < 0.05 and |logFC| ≥ 0.58 (a 1.5-fold change).",
      if (K > 2) paste0(" The **global F-test** tests $H_0: ", paste0("\\mu_{g,", latex_text(groups), "}", collapse = " = "),
                        "$ on the same fit.") else "", "\n")

  # ---- equivalent limma code ----
  add("### Equivalent limma code\n")
  contrast_args <- vapply(seq_len(ncol(cm)), function(j) {
    nm <- colnames(cm)[j]
    if (make.names(nm) != nm) nm <- paste0("`", nm, "`")
    expr <- comparisons[[j]]$contrast %||% paste(comparisons[[j]]$exp, "-", comparisons[[j]]$ctrl)
    sprintf("  %s = %s,", nm, expr)
  }, character(1))
  code <- c(
    "```r",
    "# E: log2 intensity matrix (features x samples), samples in samplesheet order",
    if (!is.null(analysis_params$group_factors))
      sprintf("samplesheet$GroupID <- paste(%s, sep = \"_\")",
              paste0("samplesheet$", analysis_params$group_factors, collapse = ", ")),
    if (length(covariates) > 0)
      sprintf("samplesheet$%s <- factor(samplesheet$%s)   # covariates are categorical", covariates, covariates),
    sprintf("design <- model.matrix(%s, data = samplesheet)", formula_label),
    "colnames(design)[1:nlevels(factor(samplesheet$GroupID))] <- levels(factor(samplesheet$GroupID))",
    if (!is.null(block)) c(
      sprintf("corfit <- duplicateCorrelation(E, design, block = samplesheet$%s)", block),
      sprintf("fit    <- lmFit(E, design, block = samplesheet$%s, correlation = corfit$consensus.correlation)", block)
    ) else "fit    <- lmFit(E, design)",
    "contrasts <- makeContrasts(",
    contrast_args,
    "  levels = design)",
    "fit2 <- eBayes(contrasts.fit(fit, contrasts), robust = TRUE)",
    sprintf("topTable(fit2, coef = \"%s\", number = Inf)   # one table per comparison", colnames(cm)[1]),
    "```"
  )
  add(paste(code, collapse = "\n"), "\n")

  paste(out, collapse = "\n")
}
