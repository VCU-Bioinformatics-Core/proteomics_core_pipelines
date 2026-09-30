# ==========================
# Single-gene visualisation functions
# ==========================

#' @title Read and Combine Per-Comparison limma Result Files
#' @details Reads each \code{<comparison>_limma.csv} file written by
#'   \code{\link{run_single_proteome_da_comparison}}, drops the row-name column
#'   added by \code{write.csv}, and adds a \code{comparison} column derived from
#'   the file name (the base name with the \code{_limma.csv} suffix removed).
#' @param limma_files Character vector of paths to limma result CSV files.
#' @return Data frame of all limma results row-bound together, with a
#'   \code{comparison} column identifying the source file.
read_limma_results <- function(limma_files) {
  missing_files <- limma_files[!file.exists(limma_files)]
  if (length(missing_files) > 0)
    stop(glue("limma file(s) not found: {paste(missing_files, collapse = ', ')}"))

  dplyr::bind_rows(lapply(limma_files, function(f) {
    res <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    # write.csv stores row names in an unnamed first column
    if (colnames(res)[1] == "") res <- res[, -1, drop = FALSE]
    res$comparison <- sub("_limma\\.csv$", "", basename(f))
    res
  }))
}

#' @title Plot Imputed Protein Levels for a Single Gene with limma Statistics
#' @details Looks up \code{gene} in \code{limma_df} by gene symbol
#'   (case-insensitive, \code{gene_name} column) or UniProt accession
#'   (\code{uniprotswissprot} column, including members of \code{;}-separated
#'   protein groups). For each matching protein, builds a strip plot of imputed
#'   log2 intensities with one strip per samplesheet group (in samplesheet order)
#'   and a mean crossbar per group, and places a table of the limma statistics
#'   from every comparison underneath. When \code{raw_matrix} is supplied,
#'   values that were missing before imputation are drawn as asterisks and a
#'   shape legend (Observed / Imputed) is added.
#' @param gene Character. Gene symbol or UniProt accession to plot.
#' @param imputed_matrix Numeric data frame or matrix (proteins x samples) of
#'   imputed log2 intensities with protein accessions as row names, e.g. the
#'   pipeline's \code{data/protein_imputed_matrix.csv}.
#' @param limma_df Data frame of limma results as returned by
#'   \code{\link{read_limma_results}}; must contain \code{uniprotswissprot},
#'   \code{comparison}, \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value},
#'   \code{adj.P.Val}, and \code{B}. \code{gene_name} is used if present.
#' @param samplesheet Data frame with \code{SampleID} and \code{GroupID} columns.
#' @param color_palette Optional character vector of colours, one per group.
#'   Default \code{NULL} uses RColorBrewer \code{"Set1"} (or ggplot's default hue
#'   palette when there are more than 9 groups).
#' @param raw_matrix Optional numeric data frame or matrix (proteins x samples)
#'   of pre-imputation log2 intensities with \code{NA} for missing values, e.g.
#'   the pipeline's \code{data/protein_raw_matrix.csv}. Used only to flag which
#'   points were imputed. Default \code{NULL} draws every point as observed.
#' @return Named list of \code{gtable} objects (plot + stats table), one per
#'   matching protein accession. Save each with \code{\link{save_plot}}.
plot_single_gene <- function(gene, imputed_matrix, limma_df, samplesheet, color_palette = NULL,
                             raw_matrix = NULL) {
  stat_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  required  <- c("uniprotswissprot", "comparison", stat_cols)
  missing_cols <- setdiff(required, colnames(limma_df))
  if (length(missing_cols) > 0)
    stop(glue("limma results are missing column(s): {paste(missing_cols, collapse = ', ')}"))

  # --- Resolve gene -> protein accession(s) ---
  gene_names <- if ("gene_name" %in% colnames(limma_df)) limma_df$gene_name else NA_character_
  by_symbol    <- !is.na(gene_names) & toupper(gene_names) == toupper(gene)
  by_accession <- vapply(strsplit(limma_df$uniprotswissprot, ";"),
                         function(ids) gene %in% ids, logical(1))
  accessions <- unique(limma_df$uniprotswissprot[by_symbol | by_accession])
  if (length(accessions) == 0)
    stop(glue("'{gene}' not found in the gene_name or uniprotswissprot columns of the limma results"))

  missing_acc <- setdiff(accessions, rownames(imputed_matrix))
  if (length(missing_acc) > 0) {
    flog.warn("Skipping accession(s) absent from the imputed matrix: %s", paste(missing_acc, collapse = ", "))
    accessions <- setdiff(accessions, missing_acc)
    if (length(accessions) == 0)
      stop(glue("No protein for '{gene}' is present in the imputed matrix"))
  }
  if (length(accessions) > 1)
    flog.info("'%s' maps to %d proteins; making one plot each: %s",
              gene, length(accessions), paste(accessions, collapse = ", "))

  # --- Sample -> group mapping in samplesheet order ---
  imputed_matrix <- align_to_samplesheet(as.data.frame(imputed_matrix), samplesheet)
  if (!is.null(raw_matrix))
    raw_matrix <- align_to_samplesheet(as.data.frame(raw_matrix), samplesheet)
  group_levels <- unique(as.character(samplesheet$GroupID))
  group_n      <- table(factor(samplesheet$GroupID, levels = group_levels))
  group_labels <- setNames(paste0(group_levels, "\n(n=", group_n, ")"), group_levels)

  if (is.null(color_palette)) {
    color_scale <- if (length(group_levels) <= 9) {
      scale_color_brewer(palette = "Set1", guide = "none")
    } else {
      scale_color_discrete(guide = "none")
    }
  } else {
    color_scale <- scale_color_manual(values = rep_len(color_palette, length(group_levels)), guide = "none")
  }

  plots <- lapply(accessions, function(acc) {
    plot_df <- data.frame(
      sample = colnames(imputed_matrix),
      group  = factor(samplesheet$GroupID, levels = group_levels),
      value  = as.numeric(unlist(imputed_matrix[acc, ]))
    )

    # Flag values that were missing before imputation
    is_imputed <- rep(FALSE, nrow(plot_df))
    if (!is.null(raw_matrix)) {
      if (acc %in% rownames(raw_matrix)) {
        is_imputed <- is.na(unlist(raw_matrix[acc, ]))
      } else {
        flog.warn("Accession %s absent from the raw matrix; treating all values as observed", acc)
      }
    }
    plot_df$status <- factor(ifelse(is_imputed, "Imputed", "Observed"), levels = c("Observed", "Imputed"))

    acc_stats <- limma_df[limma_df$uniprotswissprot == acc, , drop = FALSE]
    symbol <- if ("gene_name" %in% colnames(acc_stats)) na.omit(acc_stats$gene_name)[1] else NA
    title  <- if (!is.na(symbol)) paste0(symbol, " (", acc, ")") else acc

    p <- ggplot(plot_df, aes(x = group, y = value, color = group)) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "gray30", linewidth = 0.4) +
      geom_jitter(aes(shape = status), width = 0.15, height = 0, size = 2.5, stroke = 1, alpha = 0.8,
                  show.legend = !is.null(raw_matrix)) +  # TRUE keeps the Imputed key when unused
      color_scale +
      scale_shape_manual(values = c(Observed = 16, Imputed = 8), name = NULL, drop = FALSE,
                         guide = if (is.null(raw_matrix)) "none" else "legend") +
      scale_x_discrete(labels = group_labels) +
      labs(x = "Group", y = "Imputed log2 intensity", title = title,
           caption = if (!is.null(raw_matrix))
             paste0("Imputed values: ", sum(is_imputed), " of ", length(is_imputed)) else NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "top")

    # --- limma statistics table ---
    stats_tbl <- acc_stats[, c("comparison", stat_cols), drop = FALSE]
    stats_tbl <- stats_tbl[order(stats_tbl$adj.P.Val), , drop = FALSE]
    stats_tbl[stat_cols] <- lapply(stats_tbl[stat_cols], function(x) format(signif(x, 3)))
    colnames(stats_tbl)[1] <- "Comparison"

    tbl <- gridExtra::tableGrob(stats_tbl, rows = NULL,
                                theme = gridExtra::ttheme_minimal(base_size = 10))
    gridExtra::arrangeGrob(p, tbl, ncol = 1,
                           heights = grid::unit.c(grid::unit(1, "null"),
                                                  sum(tbl$heights) + grid::unit(0.5, "lines")))
  })
  setNames(plots, accessions)
}
