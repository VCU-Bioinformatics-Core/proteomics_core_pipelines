# ==========================
# Protein-level analysis functions
# ==========================

#' @title Generate a Protein-Level Volcano Plot
#' @details Colours points by significance and fold-change direction, labels the top 10
#'   up- and down-regulated proteins by adjusted p-value (truncating long IDs to
#'   \code{"ABC...XYZ"} format), and draws threshold lines at \code{±lfc} and at the
#'   raw p-value boundary corresponding to the significance cut-off.
#' @param data Data frame of limma results containing at minimum \code{logFC},
#'   \code{P.Value}, \code{adj.P.Val}, and the column named by \code{label_col}.
#' @param exp_name Character. Experimental group name used in axis labels.
#' @param ctrl_name Character. Control group name used in axis labels.
#' @param p_thresh Numeric. Adjusted p-value threshold for significance. Default
#'   \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change threshold. Default \code{0.58}.
#' @param sig Character. Column name to use for significance filtering. Default
#'   \code{"adj.P.Val"}.
#' @param label_col Character. Column whose values are used as point labels. Default
#'   \code{"uniprotswissprot"}.
#' @param out_dir Character. Reserved output directory parameter (currently unused;
#'   saving is the caller's responsibility). Default \code{"."}.
#' @param color1 Character. Colour for up-regulated points. Default \code{"#D55E00"}.
#' @param color2 Character. Colour for down-regulated points. Default \code{"#0072B2"}.
#' @return Named list with elements:
#'   \describe{
#'     \item{plot}{A \code{ggplot} volcano plot object.}
#'     \item{highlighted_ids}{Character vector of the labelled protein IDs (up to 20).}
#'   }
generate_volcano_protein <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                     sig = "adj.P.Val", label_col = "uniprotswissprot",
                                     out_dir = ".",
                                     color1 = "#D55E00", color2 = "#0072B2") {

  labeled_dat <- data %>%
    mutate(
      color_tag = case_when(
        eval(as.symbol(sig)) < p_thresh & logFC < -lfc ~ "Under expressed",
        eval(as.symbol(sig)) < p_thresh & logFC > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      )
    )

  top_10_genes_up <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC >= lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(.data[[label_col]])

  top_10_genes_dn <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC <= -lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(.data[[label_col]])

  labeled_dat <- labeled_dat %>%
    mutate(
      highlight = ifelse(.data[[label_col]] %in% c(top_10_genes_up, top_10_genes_dn),
                         .data[[label_col]], NA),
      label_display = ifelse(
        !is.na(highlight) & nchar(highlight) > 9,
        paste0(substr(highlight, 1, 3), "...", substr(highlight, nchar(highlight) - 2, nchar(highlight))),
        highlight
      )
    ) %>%
    filter(!is.na(logFC), !is.na(P.Value), is.finite(logFC), is.finite(P.Value))

  labeled_dat <- labeled_dat %>%
    mutate(sig_label = case_when(
      eval(as.symbol(sig)) < p_thresh & logFC >  lfc ~ "Up-Regulated",
      eval(as.symbol(sig)) < p_thresh & logFC < -lfc ~ "Down-Regulated",
      TRUE ~ "Not Significant"
    ))

  # Place the hline at the raw P.Value boundary that corresponds to adj.P.Val < p_thresh,
  # so the line aligns with the coloring. Fall back to -log10(p_thresh) if nothing is significant.
  sig_boundary <- labeled_dat %>%
    filter(sig_label != "Not Significant") %>%
    summarise(y = if (n() > 0) min(-log10(P.Value)) else -log10(p_thresh)) %>%
    pull(y)

  volcano_plot <- ggplot(labeled_dat, aes(x = logFC, y = -log10(P.Value), color = sig_label)) +
    geom_point(data = labeled_dat %>% filter(sig_label == "Not Significant"), alpha = 0.4) +
    geom_point(data = labeled_dat %>% filter(sig_label != "Not Significant"), alpha = 0.7) +
    scale_color_manual(
      values = c("Up-Regulated" = color1, "Down-Regulated" = color2, "Not Significant" = "gray70"),
      name = NULL
    ) +
    geom_label_repel(data = labeled_dat %>% filter(!is.na(label_display)),
                     aes(x = logFC, y = -log10(P.Value), label = label_display),
                     max.overlaps = Inf, show.legend = FALSE) +
    geom_hline(yintercept = sig_boundary, col = "gray40", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right",
          legend.text = element_text(size = 13),
          legend.key.size = unit(1.2, "cm"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13)) +
    labs(x = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")\n",
                    "← higher in ", ctrl_name, "          higher in ", exp_name, " →"),
         y = "-log10(P-value)",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


#' @title Generate a Protein-Level MA Plot
#' @details Plots average log2 intensity (AveExpr) on the x-axis against log2
#'   fold-change on the y-axis, coloured by significance. Proteins in
#'   \code{highlighted_ids} are labelled with the same truncation logic used in
#'   \code{\link{generate_volcano_protein}}.
#' @param data Data frame of limma results containing at minimum \code{AveExpr},
#'   \code{logFC}, \code{adj.P.Val}, and the column named by \code{label_col}.
#' @param exp_name Character. Experimental group name used in axis labels.
#' @param ctrl_name Character. Control group name used in axis labels.
#' @param highlighted_ids Character vector of protein IDs to label. Typically the
#'   \code{highlighted_ids} element returned by \code{\link{generate_volcano_protein}}.
#'   Pass \code{NULL} to skip labelling. Default \code{NULL}.
#' @param p_thresh Numeric. Adjusted p-value threshold for colouring. Default \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change for colouring. Default
#'   \code{0.58}.
#' @param label_col Character. Column used for point labels. Default
#'   \code{"uniprotswissprot"}.
#' @param color1 Character. Colour for up-regulated points. Default \code{"#D55E00"}.
#' @param color2 Character. Colour for down-regulated points. Default \code{"#0072B2"}.
#' @return A \code{ggplot} MA plot object. Not saved to disk — pass to
#'   \code{\link{save_plot}}.
generate_ma_plot_protein <- function(data, exp_name, ctrl_name, highlighted_ids = NULL,
                                      p_thresh = 0.05, lfc = 0.58,
                                      label_col = "uniprotswissprot",
                                      color1 = "#D55E00", color2 = "#0072B2") {
  ma_df <- data %>%
    filter(!is.na(AveExpr), !is.na(logFC), is.finite(AveExpr), is.finite(logFC)) %>%
    mutate(sig = case_when(
      adj.P.Val < p_thresh & logFC >  lfc ~ "Up",
      adj.P.Val < p_thresh & logFC < -lfc ~ "Down",
      TRUE ~ "NS"
    ))

  ma_df <- ma_df %>%
    mutate(sig = case_when(
      sig == "Up"   ~ "Up-Regulated",
      sig == "Down" ~ "Down-Regulated",
      TRUE          ~ "Not Significant"
    ))

  ma_labeled <- if (!is.null(highlighted_ids)) {
    ma_df %>%
      filter(.data[[label_col]] %in% highlighted_ids) %>%
      mutate(label_display = ifelse(
        nchar(.data[[label_col]]) > 9,
        paste0(substr(.data[[label_col]], 1, 3), "...",
               substr(.data[[label_col]], nchar(.data[[label_col]]) - 2, nchar(.data[[label_col]]))),
        .data[[label_col]]
      ))
  } else NULL

  p <- ggplot(ma_df, aes(x = AveExpr, y = logFC, color = sig)) +
    geom_point(data = ma_df %>% filter(sig == "Not Significant"), alpha = 0.3, size = 1) +
    geom_point(data = ma_df %>% filter(sig != "Not Significant"), alpha = 0.6, size = 1) +
    scale_color_manual(
      values = c("Up-Regulated" = color1, "Down-Regulated" = color2, "Not Significant" = "gray60"),
      name = NULL
    ) +
    geom_hline(yintercept = 0,    linetype = "dashed", color = "black") +
    geom_hline(yintercept =  lfc, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = -lfc, linetype = "dotted", color = "gray40") +
    labs(x = "Average log2 intensity (AveExpr)",
         y = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")\n",
                    "\n← higher in ", ctrl_name, "          higher in ", exp_name, " →"),
         title = create_comparison_name(exp_name, ctrl_name, "MA Plot - ")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right",
          legend.text = element_text(size = 13),
          legend.key.size = unit(1.2, "cm"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13))

  if (!is.null(ma_labeled) && nrow(ma_labeled) > 0) {
    p <- p + geom_label_repel(data = ma_labeled, aes(x = AveExpr, y = logFC, label = label_display, color = sig),
                               max.overlaps = Inf, show.legend = FALSE, size = 3)
  }
  p
}


# ==========================
# Protein-level orchestration
# ==========================

#' @title Run Full Protein-Level Differential Analysis for a Single Comparison
#' @details Orchestrates the end-to-end protein analysis pipeline for one contrast:
#'   \enumerate{
#'     \item Runs limma differential analysis via \code{\link{perform_limma_analysis}}.
#'     \item Optionally joins protein metadata (e.g. gene names).
#'     \item Optionally annotates results with imputation categories based on pre-imputation
#'       NA counts.
#'     \item Saves differential results to \code{out_dirs$de_data}.
#'     \item Generates volcano, MA, and heatmap figures.
#'     \item Optionally runs GSEA via \code{\link{process_gsea}} and saves results.
#'   }
#' @param comparison Named list with elements \code{name}, \code{exp}, and \code{ctrl}.
#' @param limma_params Named list with elements \code{E} (intensity matrix) and
#'   \code{design} (model matrix).
#' @param normalized_counts Numeric matrix or data frame (proteins x samples) of
#'   imputed log2 intensities passed to \code{\link{generate_heatmap}}.
#' @param out_dirs Named list of output directories as returned by
#'   \code{\link{setup_directories}}.
#' @param intensity_matrix_raw Numeric matrix or data frame of pre-imputation intensities
#'   (with \code{NA}s) used to compute imputation categories. Pass \code{NULL} to skip.
#' @param ont_option Character. GO ontology for GSEA: \code{"BP"} (default), \code{"MF"},
#'   or \code{"CC"}.
#' @param skip_gsea Logical. If \code{TRUE}, GSEA is skipped entirely. Default
#'   \code{FALSE}.
#' @param protein_metadata Data frame with a \code{uniprotswissprot} column for joining
#'   additional annotation (e.g. \code{PG.Genes}). Pass \code{NULL} to skip.
#' @param heatmap_norm Character. \code{"zscore"} (default) or \code{"intensity"} for
#'   the heatmap colour scale.
#' @param color1 Character. Primary plot colour (up-regulated / imputed).
#'   Default \code{"#D55E00"}.
#' @param color2 Character. Secondary plot colour (down-regulated / observed).
#'   Default \code{"#0072B2"}.
#' @return Named list with elements:
#'   \describe{
#'     \item{limma}{Data frame of annotated differential analysis results.}
#'     \item{gsea}{\code{gseaResult} object or \code{NULL} if GSEA was skipped or
#'       failed.}
#'     \item{highlighted_ids}{Character vector of protein IDs labelled on the volcano
#'       plot.}
#'   }
#'   Returns \code{NULL} if the analysis fails entirely.
run_analysis <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, ont_option = "BP", skip_gsea = FALSE, protein_metadata = NULL, heatmap_norm = "zscore", color1 = "#D55E00", color2 = "#0072B2") {
  tryCatch({
    flog.info("=== Starting analysis for comparison: %s ===", comparison$name)

    limma_results <- perform_limma_analysis(limma_params, comparison$exp, comparison$ctrl)
    limma_results <- limma_results %>% rownames_to_column(var = "uniprotswissprot")

    # Join gene names from input data (PG.Genes column)
    if (!is.null(protein_metadata)) {
      annotated_results <- limma_results %>%
        left_join(protein_metadata, by = "uniprotswissprot")
    } else {
      annotated_results <- limma_results
    }

    # Compute per-protein imputation category using pre-imputation NA counts
    if (!is.null(intensity_matrix_raw)) {
      all_samples  <- colnames(limma_params$E)
      exp_samples  <- all_samples[limma_params$design[, comparison$exp]  == 1]
      ctrl_samples <- all_samples[limma_params$design[, comparison$ctrl] == 1]
      n_exp  <- length(exp_samples)
      n_ctrl <- length(ctrl_samples)

      na_exp   <- rowSums(is.na(intensity_matrix_raw[annotated_results$uniprotswissprot, exp_samples,  drop = FALSE]))
      na_ctrl  <- rowSums(is.na(intensity_matrix_raw[annotated_results$uniprotswissprot, ctrl_samples, drop = FALSE]))
      na_total <- na_exp + na_ctrl

      annotated_results <- annotated_results %>%
        mutate(
          imputation_category = case_when(
            adj.P.Val >= 0.05                                                      ~ "not-significant",
            na_total == 0                                                           ~ "complete-data",
            (na_exp == n_exp & na_ctrl == 0) | (na_ctrl == n_ctrl & na_exp == 0)  ~ "on-off",
            na_total == 1                                                           ~ "imputation-low",
            na_total == 2                                                           ~ "imputation-medium",
            na_total >= 3                                                           ~ "imputation-high",
            TRUE                                                                    ~ "other"
          )
        )
      flog.info("Imputation category distribution for %s:\n%s",
                comparison$name,
                paste(capture.output(
                  print(table(annotated_results$imputation_category))
                ), collapse = "\n"))
    }

    output_file <- create_file_path(out_dirs$de_data, comparison$name, "_limma")
    write.csv(annotated_results, output_file)

    flog.info("Generating volcano plot: %s", comparison$name)
    volcano_result <- generate_volcano_protein(annotated_results, comparison$exp, comparison$ctrl, color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    flog.info("Generating MA plot: %s", comparison$name)
    ma_plot <- generate_ma_plot_protein(annotated_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids,
                                        color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    flog.info("Generating heatmap: %s", comparison$name)
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, normalized_counts,
                     exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                     fig_dir = out_dirs$heatmap, design = limma_params$design,
                     heatmap_norm = heatmap_norm, color1 = color1, color2 = color2)
    draw(ht, column_title = paste0(comparison$exp, " vs ", comparison$ctrl, " — Differentially Abundant Proteins"),
         column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
    dev.off()

    if (skip_gsea) {
      flog.info("Skipping GSEA for %s (--skip-gsea flag set)", comparison$name)
      gse <- NULL
    } else {
      flog.info("Running GSEA for %s", comparison$name)
      gse <- process_gsea(annotated_results, ont_option = ont_option)

      if(!is.null(gse)) {
        flog.info("GSEA returned results for %s", comparison$name)
        write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, "go_analysis_", comparison$name))

        gsea_plot <- create_barplot(gse, create_comparison_name(comparison$exp, comparison$ctrl, "GSEA "), color1 = color1, color2 = color2)

        save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                  width = 10, height = 12)
      }
    }

    return(list(limma = annotated_results, gsea = gse, highlighted_ids = volcano_result$highlighted_ids))


  }, error = function(e) {
    flog.error("run_analysis failed for comparison '%s': %s", comparison$name, e$message)
    return(NULL)
  })
}
