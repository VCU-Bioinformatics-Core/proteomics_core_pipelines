# ==========================
# Helper Functions
# ==========================
create_file_path <- function(base_dir, prefix, name, extension = ".csv") {
  file.path(base_dir, paste0(prefix, name, extension))
}

create_comparison_name <- function(exp, ctrl, prefix = "") {
  paste0(prefix, exp, " vs. ", ctrl)
}

save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  ggsave(plot, filename = filename, width = width, height = height, dpi = dpi, bg = "white")
}

export_plotly_to_html <- function(plotly_obj, file_path) {
  tryCatch({
    if (!inherits(plotly_obj, "plotly")) stop("Invalid plotly object")
    htmlwidgets::saveWidget(plotly_obj, file_path, selfcontained = TRUE)
  }, error = function(e) {
    cat(paste0("Error: ", e$message, "\n"))
  })
}

setup_directories <- function(base_dir) {
  dirs <- list(
    data = file.path(base_dir, "data"),
    figures = file.path(base_dir, "figures"),
    de_data = file.path(base_dir, "data/de_data"),
    gsea_data = file.path(base_dir, "data/gsea_data"),
    volcano = file.path(base_dir, "figures/volcano"),
    ma = file.path(base_dir, "figures/ma"),
    heatmap = file.path(base_dir, "figures/heatmap"),
    gsea = file.path(base_dir, "figures/gsea"),
    pca = file.path(base_dir, "figures/pca")
  )
  
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# ==========================
# limma analysis functions
# ==========================
perform_limma_analysis <- function(limma_params, exp, ctrl, min_valid_samples=2) {
  tryCatch({
    print(paste("Performing limma-voom analysis for", exp, "vs", ctrl))
    
    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(limma_params$design))
    
    # Fit linear model
    fit <- lmFit(limma_params$E, limma_params$design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2, robust = TRUE)
    
    # Extract results
    res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    print(paste("Analysis complete. Number of results:", nrow(res)))
    return(as.data.frame(res))
  }, error = function(e) {
    print(paste("Error in limma analysis:", e$message))
    stop(e)
  })
}

generate_volcano_protein <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                     sig = "adj.P.Val", label_col = "uniprotswissprot",
                                     out_dir = ".") {

  print('Labeling the data')
  labeled_dat <- data %>%
    mutate(
      color_tag = case_when(
        eval(as.symbol(sig)) < p_thresh & logFC < -lfc ~ "Under expressed",
        eval(as.symbol(sig)) < p_thresh & logFC > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      )
    )

  print('Extracting the top 10 genes')
  top_10_genes_up <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC >= lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(.data[[label_col]])

  print('Extracting the bottom 10 genes')
  top_10_genes_dn <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC <= -lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(.data[[label_col]])

  print('Plotting')
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

  volcano_plot <- labeled_dat %>%
    ggplot(aes(x = logFC, y = -log10(P.Value), color = color_tag)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_label_repel(aes(label = label_display), max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = c("firebrick", "steelblue")) +
    geom_hline(yintercept = -log10(p_thresh), col = "red", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme(legend.title = element_blank()) +
    labs(x = "Log2 Fold-Change (FC)",
         y = "-log10( P-value )",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


generate_volcano_phospho <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                     sig = "adj.P.Val", out_dir = ".") {

  print('Labeling the data')
  labeled_dat <- data %>%
    mutate(
      color_tag = case_when(
        eval(as.symbol(sig)) < p_thresh & logFC < -lfc ~ "Under expressed",
        eval(as.symbol(sig)) < p_thresh & logFC > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      )
    )

  print('Extracting the top 10 genes')
  top_10_genes_up <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC >= lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(peptide_id)

  print('Extracting the bottom 10 genes')
  top_10_genes_dn <- labeled_dat %>%
    filter(P.Value <= p_thresh & logFC <= -lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 10) %>%
    pull(peptide_id)

  print('Plotting')
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

  volcano_plot <- labeled_dat %>%
    ggplot(aes(x = logFC, y = -log10(P.Value), color = color_tag)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_label_repel(aes(label = label_display), max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = c("firebrick", "steelblue")) +
    geom_hline(yintercept = -log10(p_thresh), col = "red", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme(legend.title = element_blank()) +
    labs(x = "Log2 Fold-Change (FC)",
         y = "-log10( P-value )",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


generate_ma_plot_protein <- function(data, exp_name, ctrl_name, highlighted_ids = NULL,
                                      p_thresh = 0.05, lfc = 0.58,
                                      label_col = "uniprotswissprot") {
  ma_df <- data %>%
    filter(!is.na(AveExpr), !is.na(logFC), is.finite(AveExpr), is.finite(logFC)) %>%
    mutate(sig = case_when(
      adj.P.Val < p_thresh & logFC >  lfc ~ "Up",
      adj.P.Val < p_thresh & logFC < -lfc ~ "Down",
      TRUE ~ "NS"
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
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "gray60"),
                       name = NULL) +
    geom_hline(yintercept = 0,    linetype = "dashed", color = "black") +
    geom_hline(yintercept =  lfc, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = -lfc, linetype = "dotted", color = "gray40") +
    labs(x = "Average log2 intensity (AveExpr)",
         y = "log2 fold-change",
         title = create_comparison_name(exp_name, ctrl_name, "MA Plot - ")) +
    theme_bw() +
    theme(legend.position = "top")

  if (!is.null(ma_labeled) && nrow(ma_labeled) > 0) {
    p <- p + geom_label_repel(data = ma_labeled, aes(label = label_display),
                               max.overlaps = Inf, show.legend = FALSE, size = 3)
  }
  p
}


generate_ma_plot_phospho <- function(data, exp_name, ctrl_name, highlighted_ids = NULL,
                                      p_thresh = 0.05, lfc = 0.58) {
  ma_df <- data %>%
    filter(!is.na(AveExpr), !is.na(logFC), is.finite(AveExpr), is.finite(logFC)) %>%
    mutate(sig = case_when(
      adj.P.Val < p_thresh & logFC >  lfc ~ "Up",
      adj.P.Val < p_thresh & logFC < -lfc ~ "Down",
      TRUE ~ "NS"
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
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "gray60"),
                       name = NULL) +
    geom_hline(yintercept = 0,    linetype = "dashed", color = "black") +
    geom_hline(yintercept =  lfc, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = -lfc, linetype = "dotted", color = "gray40") +
    labs(x = "Average log2 intensity (AveExpr)",
         y = "log2 fold-change",
         title = create_comparison_name(exp_name, ctrl_name, "MA Plot - ")) +
    theme_bw() +
    theme(legend.position = "top")

  if (!is.null(ma_labeled) && nrow(ma_labeled) > 0) {
    p <- p + geom_label_repel(data = ma_labeled, aes(label = label_display),
                               max.overlaps = Inf, show.legend = FALSE, size = 3)
  }
  p
}


generate_global_heatmap <- function(intensity_matrix, out_dirs, top_n = 1000, molecule_label = "Proteins") {
  mat <- as.matrix(intensity_matrix)

  # Drop rows with any non-finite values before clustering
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]

  # Filter to top N most variable molecules by coefficient of variation (CV = sd / |mean|)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds   <- apply(mat, 1, sd, na.rm = TRUE)
  cv        <- row_sds / abs(row_means)
  cv[is.nan(cv) | is.infinite(cv)] <- 0

  n_select <- min(top_n, nrow(mat))
  mat <- mat[order(cv, decreasing = TRUE)[1:n_select], , drop = FALSE]

  zscores <- t(scale(t(mat)))
  zscores[is.nan(zscores) | is.infinite(zscores)] <- 0

  out_path <- file.path(out_dirs$heatmap, "global_heatmap.png")
  png(out_path, width = 2400, height = 3200, res = 300)
  ht <- Heatmap(
    zscores,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "firebrick")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_title = paste0(molecule_label, " (top ", n_select, " by CV)"),
    column_title = "Samples",
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 12)
  )
  draw(ht)
  dev.off()

  return(out_path)
}

generate_heatmap <- function(results_df, normalized_counts, p = 0.05, lfc = 0.58,
                             exp_name, ctrl_name, fig_dir, design,
                             row_id_col = "uniprotswissprot") {

  # filter the data for certain pvalues and lfc
  filtered_data <- results_df %>%
    dplyr::filter((adj.P.Val < p & abs(logFC) >= lfc))

  # restrict columns to only the two groups being compared
  all_samples  <- colnames(normalized_counts)
  exp_samples  <- all_samples[design[, exp_name]  == 1]
  ctrl_samples <- all_samples[design[, ctrl_name] == 1]
  comparison_samples <- c(exp_samples, ctrl_samples)

  # add jitter
  values <- normalized_counts[filtered_data[[row_id_col]], comparison_samples] %>%
    as.matrix() %>%
    jitter(factor = 1, amount = 0.00001)
  
  # normalize the values
  zscores <- t(scale(t(values)))
  
  # Replace NaN with row means
  print("WARNING JR: REMOVE THIS REPLACING OF NAN WITH SOMETHING MORE ACCURATE WHEN READY")
  zscores[is.na(zscores) | is.nan(zscores) | is.infinite(zscores)] <- 0
  
  # plot the heatmap
  ht <- Heatmap(
    zscores,
    name = "Z-score",  # nombre de la leyenda
    col = colorRamp2(
      c(min(zscores), 0, max(zscores)), 
      c("blue", "white", "firebrick")
    ),
    cluster_rows = TRUE,     # dendrograma de filas
    cluster_columns = TRUE,  # dendrograma de columnas
    show_row_names = FALSE,  # equivalente a labRow = NA
    show_column_names = TRUE,
    row_title = "Features",
    column_title = "Samples",
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 12)
  )
  
  return(ht)
  
}

# ==========================
# gsea analysis functions
# ==========================

aggregate_phospho_for_gsea <- function(limma_results, peptide_metadata, p = 0.05, lfc = 0.58) {
  # Join DE results with metadata to get PG.UniProtIds per peptide
  limma_results <- limma_results %>%
    left_join(peptide_metadata, by = "peptide_id")

  # Identify peptides with at least one significant result, then get their
  # parent UniProt IDs (PG.UniProtIds may be semicolon-separated; take first)
  limma_results <- limma_results %>%
    mutate(uniprot_id = sapply(strsplit(PG.UniProtIds, ";"), `[`, 1))

  sig_peptides <- limma_results %>%
    filter(adj.P.Val < p & abs(logFC) >= lfc)
  n_sig_peptides <- nrow(sig_peptides)

  sig_uniprots <- sig_peptides %>%
    pull(uniprot_id) %>%
    unique()
  n_sig_uniprots <- length(sig_uniprots)

  if (n_sig_uniprots == 0) {
    message("No significant phosphopeptides found for GSEA aggregation")
    return(NULL)
  }

  # For each qualifying UniProt ID, take the logFC of the most significant
  # phosphopeptide (smallest adj.P.Val) as the representative value
  aggregated <- limma_results %>%
    filter(uniprot_id %in% sig_uniprots) %>%
    group_by(uniprot_id) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(uniprot_id, logFC)

  # Map UniProt IDs to Ensembl IDs via biomart (same approach as regular pipeline)
  mapping <- getBM(
    attributes = c("uniprotswissprot", "ensembl_gene_id"),
    filters    = "uniprotswissprot",
    values     = unique(aggregated$uniprot_id),
    mart       = ensembl
  ) %>%
    distinct(uniprotswissprot, .keep_all = TRUE)

  aggregated <- aggregated %>%
    left_join(mapping, by = c("uniprot_id" = "uniprotswissprot")) %>%
    filter(!is.na(ensembl_gene_id))
  n_mapped <- nrow(aggregated)

  if (n_mapped == 0) {
    message("No Ensembl IDs found for aggregated phospho UniProt IDs")
    return(NULL)
  }

  counts <- list(n_sig_peptides = n_sig_peptides, n_sig_uniprots = n_sig_uniprots, n_mapped = n_mapped)
  return(list(data = aggregated, counts = counts))
}

process_gsea <- function(result, p = 0.05, lfc = 0.58, ont_option = "BP") {
  tryCatch({
    
    sig_genes <- result %>% arrange(desc(logFC))
    sig_genes <- sig_genes %>% filter(!is.na(logFC))
    
    if(nrow(sig_genes) < 2) {
      message("Not enough significant genes for GSEA analysis")
      return(NULL)
    }
    
    # create a gene list where the names are ensembl genes
    # and the values are logFC
    gene_list <- sig_genes$logFC
    names(gene_list) <- sig_genes$ensembl_gene_id
    gene_list <- gene_list[!is.na(names(gene_list))]
    gene_list <- gene_list[!duplicated(names(gene_list))]
    
    valid_genes <- bitr(names(gene_list),
                        fromType = "ENSEMBL",
                        toType = "ENTREZID",
                        OrgDb = annotation_db)
    gene_list <- gene_list[names(gene_list) %in% valid_genes$ENSEMBL]
    
    gse_result <- gseGO(geneList = gene_list,
                        ont = ont_option,
                        minGSSize = 1,
                        maxGSSize = 900,
                        keyType = "ENSEMBL",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        OrgDb = annotation_db)
    
    if(nrow(gse_result@result) == 0) {
      message("No enriched terms found in GSEA analysis")
      return(NULL)
    }
    
    setReadable(gse_result, OrgDb = annotation_db, keyType = "ENSEMBL")
  }, error = function(e) {
    message("Error in GSEA processing: ", e$message)
    return(NULL)
  })
}

create_dotplot <- function(gse, title) {
  dotplot(gse,
          showCategory = 15,
          title = title,
          split = ".sign",
          orderBy = "p.adjust",
          label_format = 31,
          font.size = 9) +
    facet_grid(.~.sign)
}

# ==========================
# centralizing function
# ==========================

run_analysis <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, ont_option = "BP") {
  tryCatch({
    print(paste("\nStarting analysis for comparison:", comparison$name))

    print('Running limma')
    limma_results <- perform_limma_analysis(limma_params, comparison$exp, comparison$ctrl)
    limma_results <- limma_results %>% rownames_to_column(var = "uniprotswissprot")

    print('Annotating the results')
    mapping <- getBM(
      attributes = c("uniprotswissprot", "ensembl_gene_id"), # "external_gene_name"
      filters = "uniprotswissprot",
      values = uniprot_clean,
      mart = ensembl
    )
    mapping <- mapping %>% distinct(uniprotswissprot, .keep_all = TRUE)
    annotated_results <- limma_results %>% left_join(mapping, by=join_by(uniprotswissprot))

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
      print("Imputation category distribution:")
      print(table(annotated_results$imputation_category))
    }

    output_file <- create_file_path(out_dirs$de_data, "limma_", comparison$name)
    write.csv(annotated_results, output_file)
    
    print('Generating the volcano plot')
    volcano_result <- generate_volcano_protein(annotated_results, comparison$exp, comparison$ctrl)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    print('Generating the MA plot')
    ma_plot <- generate_ma_plot_protein(annotated_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    print('Generating the heatmap')
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, intensity_matrix,
                     exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                     fig_dir = out_dirs$heatmap, design = limma_params$design)
    draw(ht)
    dev.off()
    
    print('Running GSEA')
    gse <- NULL # process_gsea(annotated_results, ont_option = ont_option)
    #gse <- process_gsea(annotated_results, ont_option = ont_option)
    
    if(!is.null(gse)) {
      print('GSE has data')
      write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, "GO_Analysis_", comparison$name))
      
      print('create dotplot')
      gsea_plot <- create_dotplot(gse, create_comparison_name(comparison$exp, comparison$ctrl, "GSEA "))
      
      print('save plot')
      save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_GSEA.png"),
                width = 10, height = 12)
    }

    return(list(limma = annotated_results, gsea = gse, highlighted_ids = volcano_result$highlighted_ids))
    
    
  }, error = function(e) {
    message("Error in run_analysis: ", e$message)
    return(NULL)
  })
}

run_analysis_phospho <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, peptide_metadata = NULL, ont_option = "BP") {
  tryCatch({
    print(paste("Starting analysis for comparison:", comparison$name))

    print('Running limma')
    limma_results <- perform_limma_analysis(limma_params, comparison$exp, comparison$ctrl)
    limma_results <- limma_results %>% rownames_to_column(var = "peptide_id")

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
      print("Imputation category distribution:")
      print(table(limma_results$imputation_category))
    }

    output_file <- create_file_path(out_dirs$de_data, "limma_", comparison$name)
    write.csv(limma_results, output_file)

    print('Generating the volcano plot')
    volcano_result <- generate_volcano_phospho(limma_results, comparison$exp, comparison$ctrl)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    print('Generating the MA plot')
    ma_plot <- generate_ma_plot_phospho(limma_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    print('Generating the heatmap')
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, normalized_counts,
                           exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                           fig_dir = out_dirs$heatmap, design = limma_params$design,
                           row_id_col = "peptide_id")
    draw(ht)
    dev.off()
    
    print('Running GSEA')
    agg_result <- tryCatch(
      aggregate_phospho_for_gsea(limma_results, peptide_metadata),
      error = function(e) { message("Error in phospho aggregation: ", e$message); NULL }
    )
    protein_counts <- if (!is.null(agg_result)) agg_result$counts else NULL
    gse <- tryCatch({
      if (!is.null(agg_result)) process_gsea(agg_result$data, ont_option = ont_option) else NULL
    }, error = function(e) {
      message("Error in phospho GSEA: ", e$message)
      NULL
    })
    if(!is.null(gse)) {
      print('GSE has data')
      tryCatch({
        write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, "GO_Analysis_", comparison$name))

        print('create dotplot')
        gsea_plot <- create_dotplot(gse, create_comparison_name(comparison$ctrl, comparison$exp, "GSEA "))

        print('save plot')
        save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_GSEA.png"),
                  width = 10, height = 12)
      }, error = function(e) {
        message("Error saving phospho GSEA outputs: ", e$message)
      })
    }

    return(list(limma = limma_results, gsea = gse, protein_counts = protein_counts, highlighted_ids = volcano_result$highlighted_ids))
    
    
  }, error = function(e) {
    message("Error in run_analysis: ", e$message)
    return(NULL)
  })
}
