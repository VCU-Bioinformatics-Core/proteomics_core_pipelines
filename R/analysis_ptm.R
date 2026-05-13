# ==========================
# PTM-level analysis functions
# ==========================
library(tidyr)
library(dplyr)
library(progress)

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







# ==============================================================
# PTM-level orchestration for PTM and protein level integrations
# ==============================================================

calc_lm_stats <- function(exp_vals, ctrl_vals, contrast_name = "exp_vs_ctrl") {
  # combine into a data frame for lm
  df <- data.frame(
    value = c(exp_vals, ctrl_vals),
    group = factor(c(rep("exp", length(exp_vals)), rep("ctrl", length(ctrl_vals))),
                   levels = c("ctrl", "exp"))  # ctrl is reference
  )
  
  # drop NAs (from imputation failures etc.)
  df <- df[is.finite(df$value), ]
  
  # need at least 2 groups with observations
  if (length(unique(df$group)) < 2 || sum(df$group == "exp") < 1 || sum(df$group == "ctrl") < 1) {
    return(list(logFC = NA, SE = NA, DF = NA, t = NA, P.Value = NA))
  }
  
  fit <- lm(value ~ group, data = df)
  coef_summary <- summary(fit)$coefficients
  
  # "groupexp" row is the exp vs ctrl contrast
  logFC  <- coef_summary["groupexp", "Estimate"]
  SE     <- coef_summary["groupexp", "Std. Error"]
  DF     <- fit$df.residual
  t      <- coef_summary["groupexp", "t value"]
  pval   <- coef_summary["groupexp", "Pr(>|t|)"]

  list(logFC = logFC, SE = SE, DF = DF, t = t, P.Value = pval)
}

perform_ptm_protein_diff_analysis <- function(comparison, exp_samples, ctrl_samples,
                                              ptm_imp_matrix, protein_imp_matrix){
  results <- list()
  uniq_uniprot_ids <- unique(protein_imp_matrix$PG.ProteinAccessions)
  uniq_uniprot_ids <- head(uniq_uniprot_ids, 50)  # DEBUG: limit to 500
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta",
    total = length(uniq_uniprot_ids), clear = FALSE
  )
  for (uniprot_id in uniq_uniprot_ids){
    pb$tick()

    # extract the protein subset
    curr_prots <- protein_imp_matrix %>% filter(PG.ProteinAccessions == uniprot_id)
    curr_prots.exp <- curr_prots %>% filter(sample %in% exp_samples)
    curr_prots.ctrl <- curr_prots %>% filter(sample %in% ctrl_samples)

    # calculating the protein-level stats
    curr_prot_stats <- calc_lm_stats(curr_prots.exp$intensity, curr_prots.ctrl$intensity)

    # get the ptms of the current protein accession
    curr_ptms <- ptm_imp_matrix %>% filter(PG.UniProtIds == uniprot_id)
    uniq_ptms <- unique(curr_ptms$peptide_id)

    # cycle through the unique ptms
    for (ptm in uniq_ptms){

      #cat("\t", glue("ptm: {ptm}"), "\n")

      # extract the ptm subset
      curr_ptms.exp <- curr_ptms %>% filter(sample %in% exp_samples, peptide_id == ptm)
      curr_ptms.ctrl <- curr_ptms %>% filter(sample %in% ctrl_samples, peptide_id == ptm)

      # calculating the PTM-level stats
      curr_ptm_stats <- calc_lm_stats(curr_ptms.exp$intensity, curr_ptms.ctrl$intensity)

      # perform MSstatsPTM-style adjusted test
      logfc <- curr_ptm_stats$logFC - curr_prot_stats$logFC

      s2_ptm  <- curr_ptm_stats$SE ^ 2
      s2_prot <- curr_prot_stats$SE ^ 2

      stderr <- sqrt(s2_ptm + s2_prot)
      numer  <- (s2_ptm + s2_prot) ^ 2
      denom  <- (s2_ptm ^ 2 / curr_ptm_stats$DF) + (s2_prot ^ 2 / curr_prot_stats$DF)
      df     <- numer / denom

      t      <- logfc / stderr
      pval   <- 2 * pt(abs(t), df, lower.tail = FALSE)

      ave_expr <- mean(c(curr_ptms.exp$intensity, curr_ptms.ctrl$intensity), na.rm = TRUE)

      results[[glue("{uniprot_id}__{ptm}")]] <- data.frame(
        uniprot_id      = uniprot_id,
        peptide_id      = ptm,
        logFC_ptm       = curr_ptm_stats$logFC,
        logFC_prot      = curr_prot_stats$logFC,
        logFC           = logfc,
        AveExpr         = ave_expr,
        SE_ptm          = curr_ptm_stats$SE,
        SE_prot         = curr_prot_stats$SE,
        DF_ptm          = curr_ptm_stats$DF,
        DF_prot         = curr_prot_stats$DF,
        t               = t,
        P.Value         = pval
      )
    }
    #cat("\n")
  }
  results_df <- bind_rows(results)
  results_df <- bind_rows(results) %>% mutate(adj.P.Val = p.adjust(P.Value, method = "BH"))
  return(results_df)
}

run_ptm_integration_analysis <- function(comparison, limma_params, ptm_imp_matrix,
                                          protein_imp_matrix, out_dirs, intensity_matrix_raw = NULL,
                                          peptide_metadata = NULL, ont_option = "BP", ensembl = NULL,
                                          skip_gsea = FALSE, heatmap_norm = "zscore",
                                          color1 = "#D55E00", color2 = "#0072B2") {




  tryCatch({

    flog.info("=== Starting PTM analysis for comparison: %s ===", comparison$name)

    # define the list of samples for ctrl and exp
    all_samples  <- colnames(limma_params$E)
    exp_samples  <- all_samples[limma_params$design[, comparison$exp]  == 1]
    ctrl_samples <- all_samples[limma_params$design[, comparison$ctrl] == 1]
    n_exp  <- length(exp_samples)
    n_ctrl <- length(ctrl_samples)

    # remove rows with semicolons in the protein ID column ──────────────────
    ptm_imp_matrix_clean <- ptm_imp_matrix[!grepl(";", ptm_imp_matrix$PG.UniProtIds), ]
    protein_imp_matrix_clean <- protein_imp_matrix[!grepl(";", protein_imp_matrix$PG.ProteinAccessions), ]

    # find shared IDs and keep only overlapping rows in each ────────────────
    shared_ids <- intersect(ptm_imp_matrix_clean$PG.UniProtIds,
                            protein_imp_matrix_clean$PG.ProteinAccessions)
    ptm_imp_matrix_filtered    <- ptm_imp_matrix_clean[ptm_imp_matrix_clean$PG.UniProtIds %in% shared_ids, ]
    protein_imp_matrix_filtered <- protein_imp_matrix_clean[protein_imp_matrix_clean$PG.ProteinAccessions %in% shared_ids, ]

    # pivot the ptm data
    ptm_imp_matrix_filtered_long <- ptm_imp_matrix_filtered %>%
      pivot_longer(
        cols      = -c(peptide_id, PG.UniProtIds),  # pivot everything except these two
        names_to  = "sample",
        values_to = "intensity"
      )

    # pivot the protein data
    protein_imp_matrix_filtered_long <- protein_imp_matrix_filtered %>%
      pivot_longer(
        cols      = -c(PG.ProteinAccessions),  # pivot everything except these two
        names_to  = "sample",
        values_to = "intensity"
      )

    diff_results <- perform_ptm_protein_diff_analysis(comparison, exp_samples, ctrl_samples, ptm_imp_matrix_filtered_long, protein_imp_matrix_filtered_long)
    
    # Join peptide metadata to get UniProt IDs, then annotate with gene symbols
    if (!is.null(peptide_metadata)) {

      diff_results <- diff_results %>% left_join(peptide_metadata, by="peptide_id") %>%
        mutate(uniprot_id = sapply(strsplit(PG.UniProtIds, ";"), `[`, 1))
      
      annotation <- tryCatch({
        m <- getBM(
          attributes = c("uniprotswissprot", "hgnc_symbol"),
          filters    = "uniprotswissprot",
          values     = unique(na.omit(diff_results$uniprot_id)),
          mart       = ensembl
        )
        distinct(m, uniprotswissprot, .keep_all = TRUE)
      }, error = function(e) {
        flog.error("BioMart annotation failed for PTM comparison '%s': %s",
                   comparison$name, e$message)
        data.frame(uniprotswissprot = character(), hgnc_symbol = character())
      })
      diff_results <- diff_results %>%
        left_join(annotation, by = c("uniprot_id" = "uniprotswissprot"))
    }

    # Compute per-protein imputation category using pre-imputation NA counts
    if (!is.null(intensity_matrix_raw)) {      
      na_exp   <- rowSums(is.na(intensity_matrix_raw[diff_results$peptide_id, exp_samples,  drop = FALSE]))
      na_ctrl  <- rowSums(is.na(intensity_matrix_raw[diff_results$peptide_id, ctrl_samples, drop = FALSE]))
      na_total <- na_exp + na_ctrl
      
      diff_results <- diff_results %>%
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
                  print(table(diff_results$imputation_category))
                ), collapse = "\n"))
    }
    
    output_file <- create_file_path(out_dirs$de_data, comparison$name, "_limma")
    write.csv(diff_results, output_file)
    
    flog.info("Generating PTM volcano plot: %s", comparison$name)
    volcano_result <- generate_volcano_ptm(diff_results, comparison$exp, comparison$ctrl, 
                                          logFC_colname="logFC", color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)
    
    flog.info("Generating PTM MA plot: %s", comparison$name)
    ma_plot <- generate_ma_plot_ptm(diff_results, comparison$exp, comparison$ctrl,
                                    highlighted_ids = volcano_result$highlighted_ids,
                                    color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)
    
    flog.info("Generating PTM heatmap: %s", comparison$name)
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)

    heatmap_df <- ptm_imp_matrix_filtered %>% select(-c(peptide_id, PG.UniProtIds))
    ht <- generate_heatmap(diff_results, heatmap_df,
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
        aggregate_ptm_for_gsea(diff_results, peptide_metadata),
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
    
    return(list(limma = diff_results, gsea = gse, protein_counts = protein_counts, highlighted_ids = volcano_result$highlighted_ids))

  }, error = function(e) {
    flog.error("run_analysis_ptm_integrated failed for comparison '%s': %s\n%s",
               comparison$name, e$message, paste(deparse(e$call), collapse = " "))
    return(NULL)
  })



  
}







