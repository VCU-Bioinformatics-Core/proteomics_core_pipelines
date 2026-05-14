# ==========================
# PTM-level analysis functions
# ==========================
library(tidyr)
library(dplyr)
library(progress)

#' @title Generate a PTM-Level Volcano Plot
#' @details Colours points by significance and fold-change direction, and labels the top
#'   10 up- and down-regulated PTM peptides by adjusted p-value. Labels split the
#'   \code{peptide_id} on \code{" -- "} to show a protein prefix plus a truncated
#'   precursor suffix (e.g. \code{"P12345 -- AB...YZ"}).
#' @param data Data frame of PTM differential results containing at minimum \code{logFC},
#'   \code{P.Value}, \code{adj.P.Val}, and \code{peptide_id}.
#' @param exp_name Character. Experimental group name used in axis labels.
#' @param ctrl_name Character. Control group name used in axis labels.
#' @param p_thresh Numeric. Adjusted p-value threshold for significance. Default
#'   \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change threshold. Default \code{0.58}.
#' @param sig Character. Column name to use for significance filtering. Default
#'   \code{"adj.P.Val"}.
#' @param logFC_colname Character. Column name for log2 fold-change values. Default
#'   \code{"logFC"}.
#' @param out_dir Character. Reserved output directory parameter (saving is the caller's
#'   responsibility). Default \code{"."}.
#' @param color1 Character. Colour for up-regulated points. Default \code{"#D55E00"}.
#' @param color2 Character. Colour for down-regulated points. Default \code{"#0072B2"}.
#' @return Named list with elements:
#'   \describe{
#'     \item{plot}{A \code{ggplot} volcano plot object.}
#'     \item{highlighted_ids}{Character vector of the labelled peptide IDs (up to 20).}
#'   }
generate_volcano_ptm <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                    sig = "adj.P.Val", logFC_colname = "logFC", out_dir = ".",
                                    color1 = "#D55E00", color2 = "#0072B2") {

  labeled_dat <- data %>%
    mutate(
      color_tag = case_when(
        .data[[sig]] < p_thresh & .data[[logFC_colname]] < -lfc ~ "Under expressed",
        .data[[sig]] < p_thresh & .data[[logFC_colname]] > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      )
    )
  
  top_10_genes_up <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC >= lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(peptide_id)

  top_10_genes_dn <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC <= -lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(peptide_id)

  labeled_dat <- labeled_dat %>%
    mutate(
      highlight = ifelse(peptide_id %in% c(top_10_genes_up, top_10_genes_dn), peptide_id, NA),
      label_display = ifelse(
        !is.na(highlight),
        {
          parts <- strsplit(highlight, " -- ", fixed = TRUE)
          protein_label <- sapply(parts, `[`, 1)
          precursor_part <- sapply(parts, function(x) paste(x[-1], collapse = " -- "))
          ifelse(
            nchar(precursor_part) > 10,
            paste0(protein_label, " -- ", substr(precursor_part, 1, 2), "...",
                   substr(precursor_part, nchar(precursor_part) - 1, nchar(precursor_part))),
            paste0(protein_label, " -- ", precursor_part)
          )
        },
        NA_character_
      )
    ) %>%
    filter(!is.na(logFC), !is.na(P.Value), is.finite(logFC), is.finite(P.Value))

  labeled_dat <- labeled_dat %>%
    mutate(sig_label = case_when(
      .data[[sig]] < p_thresh & .data[[logFC_colname]] >  lfc ~ "Up-Regulated",
      .data[[sig]] < p_thresh & .data[[logFC_colname]] < -lfc ~ "Down-Regulated",
      TRUE ~ "Not Significant"
    ))

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


#' @title Generate a PTM-Level MA Plot
#' @details Plots average log2 intensity (AveExpr) against log2 fold-change, coloured by
#'   significance. Peptides in \code{highlighted_ids} are labelled using the same
#'   \code{"protein -- precursor"} truncation logic as
#'   \code{\link{generate_volcano_ptm}}.
#' @param data Data frame of PTM differential results containing at minimum \code{AveExpr},
#'   \code{logFC}, \code{adj.P.Val}, and \code{peptide_id}.
#' @param exp_name Character. Experimental group name used in axis labels.
#' @param ctrl_name Character. Control group name used in axis labels.
#' @param highlighted_ids Character vector of peptide IDs to label. Typically the
#'   \code{highlighted_ids} element returned by \code{\link{generate_volcano_ptm}}.
#'   Pass \code{NULL} to skip labelling. Default \code{NULL}.
#' @param p_thresh Numeric. Adjusted p-value threshold for colouring. Default \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change for colouring. Default
#'   \code{0.58}.
#' @param color1 Character. Colour for up-regulated points. Default \code{"#D55E00"}.
#' @param color2 Character. Colour for down-regulated points. Default \code{"#0072B2"}.
#' @return A \code{ggplot} MA plot object. Not saved to disk — pass to
#'   \code{\link{save_plot}}.
generate_ma_plot_ptm <- function(data, exp_name, ctrl_name, highlighted_ids = NULL,
                                      p_thresh = 0.05, lfc = 0.58,
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
      filter(peptide_id %in% highlighted_ids) %>%
      mutate(label_display = {
        parts <- strsplit(peptide_id, " -- ", fixed = TRUE)
        protein_label <- sapply(parts, `[`, 1)
        precursor_part <- sapply(parts, function(x) paste(x[-1], collapse = " -- "))
        ifelse(
          nchar(precursor_part) > 10,
          paste0(protein_label, " -- ", substr(precursor_part, 1, 2), "...",
                 substr(precursor_part, nchar(precursor_part) - 1, nchar(precursor_part))),
          paste0(protein_label, " -- ", precursor_part)
        )
      })
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
                    "↑ higher in ", exp_name, "\n↓ higher in ", ctrl_name),
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


#' @title Aggregate PTM Peptide Results to Protein Level for GSEA
#' @details Filters to significant PTM peptides (\code{adj.P.Val < p} and
#'   \code{|logFC| >= lfc}), then collapses to one representative logFC per UniProt ID
#'   by selecting the peptide with the smallest adjusted p-value. Maps UniProt IDs to
#'   Ensembl gene IDs via BioMart and drops entries without a mapping.
#'
#'   Requires the \code{ensembl} mart object to be available in the calling environment
#'   (not passed as a parameter — relies on lexical scoping from the pipeline runner).
#' @param limma_results Data frame of PTM limma results. Must contain \code{adj.P.Val},
#'   \code{logFC}, and \code{peptide_id}. If \code{uniprot_id} is absent it is derived
#'   by joining \code{peptide_metadata} and splitting \code{PG.UniProtIds} on \code{";"}.
#' @param peptide_metadata Data frame with columns \code{peptide_id} and
#'   \code{PG.UniProtIds} used as a fallback join when \code{uniprot_id} is missing from
#'   \code{limma_results}.
#' @param p Numeric. Adjusted p-value threshold. Default \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change threshold. Default \code{0.58}.
#' @return Named list with elements:
#'   \describe{
#'     \item{data}{Data frame with columns \code{uniprot_id}, \code{logFC},
#'       \code{ensembl_gene_id}, and \code{hgnc_symbol} ready for
#'       \code{\link{process_gsea}}.}
#'     \item{counts}{Named list of diagnostic counts: \code{n_sig_peptides},
#'       \code{n_sig_uniprots}, \code{n_mapped}.}
#'   }
#'   Returns \code{NULL} if no significant peptides are found, BioMart fails, or no
#'   Ensembl IDs are mapped.
aggregate_ptm_for_gsea <- function(limma_results, peptide_metadata, p = 0.05, lfc = 0.58) {
  # uniprot_id is already present on limma_results (joined upstream in
  # run_analysis_ptm); if for any reason it is absent, fall back to
  # joining peptide_metadata and deriving it from PG.UniProtIds
  if (!"uniprot_id" %in% colnames(limma_results)) {
    limma_results <- limma_results %>%
      left_join(peptide_metadata, by = "peptide_id") %>%
      mutate(uniprot_id = sapply(strsplit(PG.UniProtIds, ";"), `[`, 1))
  }

  sig_peptides <- limma_results %>%
    filter(adj.P.Val < p & abs(logFC) >= lfc)
  n_sig_peptides <- nrow(sig_peptides)

  sig_uniprots <- sig_peptides %>%
    pull(uniprot_id) %>%
    unique()
  n_sig_uniprots <- length(sig_uniprots)

  if (n_sig_uniprots == 0) {
    flog.warn("No significant PTM peptides found for GSEA aggregation")
    return(NULL)
  }

  # For each qualifying UniProt ID, take the logFC of the most significant
  # PTM peptide (smallest adj.P.Val) as the representative value
  aggregated <- limma_results %>%
    filter(uniprot_id %in% sig_uniprots) %>%
    group_by(uniprot_id) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(uniprot_id, logFC)

  # Map UniProt IDs to Ensembl IDs via biomart (same approach as regular pipeline)
  if (is.null(ensembl)) {
    flog.warn("BioMart mart object is NULL — skipping Ensembl mapping for PTM GSEA")
    return(NULL)
  }
  mapping <- tryCatch(
    getBM(
      attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol"),
      filters    = "uniprotswissprot",
      values     = unique(aggregated$uniprot_id),
      mart       = ensembl
    ) %>% distinct(uniprotswissprot, .keep_all = TRUE),
    error = function(e) {
      flog.error("BioMart query failed in PTM GSEA aggregation: %s",
                 e$message)
      NULL
    }
  )
  if (is.null(mapping)) return(NULL)

  aggregated <- aggregated %>%
    left_join(mapping, by = c("uniprot_id" = "uniprotswissprot")) %>%
    filter(!is.na(ensembl_gene_id))
  n_mapped <- nrow(aggregated)

  if (n_mapped == 0) {
    flog.warn("No Ensembl IDs found for aggregated PTM UniProt IDs — GSEA skipped")
    return(NULL)
  }

  counts <- list(n_sig_peptides = n_sig_peptides, n_sig_uniprots = n_sig_uniprots, n_mapped = n_mapped)
  return(list(data = aggregated, counts = counts))
}


# ==========================
# PTM-level orchestration for PTM only datasets
# ==========================

#' @title Run Full PTM-Only Differential Analysis for a Single Comparison
#' @details Orchestrates the end-to-end PTM analysis pipeline for one contrast using
#'   standard limma (no protein-level correction — use
#'   \code{run_single_ptm_ida_comparison} for MSstatsPTM-style correction):
#'   \enumerate{
#'     \item Runs limma differential analysis via \code{\link{perform_limma_analysis}}.
#'     \item Optionally joins peptide metadata and BioMart gene symbols.
#'     \item Optionally annotates results with imputation categories.
#'     \item Saves differential results to \code{out_dirs$de_data}.
#'     \item Generates volcano, MA, and heatmap figures.
#'     \item Optionally aggregates to protein level and runs GSEA.
#'   }
#' @param comparison Named list with elements \code{name}, \code{exp}, and \code{ctrl}.
#' @param limma_params Named list with elements \code{E} (intensity matrix) and
#'   \code{design} (model matrix).
#' @param normalized_counts Numeric matrix or data frame (PTM peptides x samples) of
#'   imputed log2 intensities passed to \code{\link{generate_heatmap}}.
#' @param out_dirs Named list of output directories as returned by
#'   \code{\link{setup_directories}}.
#' @param intensity_matrix_raw Numeric matrix or data frame of pre-imputation intensities
#'   (with \code{NA}s) used to compute imputation categories. Pass \code{NULL} to skip.
#' @param peptide_metadata Data frame with columns \code{peptide_id} and
#'   \code{PG.UniProtIds} for annotation joins. Pass \code{NULL} to skip.
#' @param ont_option Character. GO ontology for GSEA: \code{"BP"} (default), \code{"MF"},
#'   or \code{"CC"}.
#' @param ensembl A \code{Mart} object from \code{biomaRt::useMart} for gene symbol
#'   annotation. Pass \code{NULL} to skip BioMart lookup.
#' @param skip_gsea Logical. If \code{TRUE}, GSEA is skipped entirely. Default
#'   \code{FALSE}.
#' @param heatmap_norm Character. \code{"zscore"} (default) or \code{"intensity"} for
#'   the heatmap colour scale.
#' @param color1 Character. Primary plot colour. Default \code{"#D55E00"}.
#' @param color2 Character. Secondary plot colour. Default \code{"#0072B2"}.
#' @return Named list with elements:
#'   \describe{
#'     \item{limma}{Data frame of annotated differential analysis results.}
#'     \item{gsea}{\code{gseaResult} object or \code{NULL}.}
#'     \item{protein_counts}{Named list of GSEA aggregation counts, or \code{NULL}.}
#'     \item{highlighted_ids}{Character vector of peptide IDs labelled on the volcano
#'       plot.}
#'   }
#'   Returns \code{NULL} if the analysis fails entirely.
run_analysis_ptm <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, peptide_metadata = NULL, ont_option = "BP", ensembl=NULL, skip_gsea = FALSE, heatmap_norm = "zscore", color1 = "#D55E00", color2 = "#0072B2") {
  tryCatch({
    flog.info("=== Starting PTM analysis for comparison: %s ===", comparison$name)

    limma_results <- perform_limma_analysis(limma_params, comparison$exp, comparison$ctrl)
    limma_results <- limma_results %>% rownames_to_column(var = "peptide_id")

    # Join peptide metadata to get UniProt IDs, then annotate with gene symbols
    if (!is.null(peptide_metadata)) {
      limma_results <- limma_results %>% left_join(peptide_metadata, by = "peptide_id") %>%
        mutate(uniprot_id = sapply(strsplit(PG.UniProtIds, ";"), `[`, 1))
      annotation <- tryCatch({
        m <- getBM(
          attributes = c("uniprotswissprot", "hgnc_symbol"),
          filters    = "uniprotswissprot",
          values     = unique(na.omit(limma_results$uniprot_id)),
          mart       = ensembl
        )
        distinct(m, uniprotswissprot, .keep_all = TRUE)
      }, error = function(e) {
        flog.error("BioMart annotation failed for PTM comparison '%s': %s",
                   comparison$name, e$message)
        data.frame(uniprotswissprot = character(), hgnc_symbol = character())
      })
      limma_results <- limma_results %>%
        left_join(annotation, by = c("uniprot_id" = "uniprotswissprot"))
    }

    # Compute per-protein imputation category using pre-imputation NA counts
    if (!is.null(intensity_matrix_raw)) {
      all_samples  <- colnames(limma_params$E)
      exp_samples  <- all_samples[limma_params$design[, comparison$exp]  == 1]
      ctrl_samples <- all_samples[limma_params$design[, comparison$ctrl] == 1]
      n_exp  <- length(exp_samples)
      n_ctrl <- length(ctrl_samples)

      na_exp   <- rowSums(is.na(intensity_matrix_raw[limma_results$peptide_id, exp_samples,  drop = FALSE]))
      na_ctrl  <- rowSums(is.na(intensity_matrix_raw[limma_results$peptide_id, ctrl_samples, drop = FALSE]))
      na_total <- na_exp + na_ctrl

      limma_results <- limma_results %>%
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
                  print(table(limma_results$imputation_category))
                ), collapse = "\n"))
    }

    output_file <- create_file_path(out_dirs$de_data, comparison$name, "_limma")
    write.csv(limma_results, output_file)

    flog.info("Generating PTM volcano plot: %s", comparison$name)
    volcano_result <- generate_volcano_ptm(limma_results, comparison$exp, comparison$ctrl, color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    flog.info("Generating PTM MA plot: %s", comparison$name)
    ma_plot <- generate_ma_plot_ptm(limma_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids,
                                        color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    flog.info("Generating PTM heatmap: %s", comparison$name)
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, normalized_counts,
                           exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                           fig_dir = out_dirs$heatmap, design = limma_params$design,
                           row_id_col = "peptide_id",
                           heatmap_norm = heatmap_norm, color1 = color1, color2 = color2)
    draw(ht, column_title = paste0(comparison$exp, " vs ", comparison$ctrl, " — Differentially Abundant PTM Peptides"),
         column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
    dev.off()

    if (skip_gsea) {
      flog.info("Skipping GSEA for %s (--skip-gsea flag set)", comparison$name)
      gse <- NULL
      protein_counts <- NULL
    } else {
      flog.info("Running PTM GSEA for %s", comparison$name)
      agg_result <- tryCatch(
        aggregate_ptm_for_gsea(limma_results, peptide_metadata),
        error = function(e) {
          flog.error("PTM GSEA aggregation failed for '%s': %s",
                     comparison$name, e$message)
          NULL
        }
      )
      protein_counts <- if (!is.null(agg_result)) agg_result$counts else NULL
      gse <- tryCatch({
        if (!is.null(agg_result)) process_gsea(agg_result$data, ont_option = ont_option) else NULL
      }, error = function(e) {
        flog.error("PTM GSEA failed for '%s': %s", comparison$name, e$message)
        NULL
      })
      if (!is.null(gse)) {
        flog.info("PTM GSEA returned results for %s", comparison$name)
        tryCatch({
          write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, comparison$name, "_go_analysis.csv"))
          gsea_plot <- create_barplot(gse, create_comparison_name(comparison$ctrl, comparison$exp, "GSEA "), color1 = color1, color2 = color2)
          save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                    width = 10, height = 12)
        }, error = function(e) {
          flog.error("Failed to save PTM GSEA outputs for '%s': %s",
                     comparison$name, e$message)
        })
      }
    }

    return(list(limma = limma_results, gsea = gse, protein_counts = protein_counts, highlighted_ids = volcano_result$highlighted_ids))

  }, error = function(e) {
    flog.error("run_analysis_ptm failed for comparison '%s': %s",
               comparison$name, e$message)
    return(NULL)
  })
}