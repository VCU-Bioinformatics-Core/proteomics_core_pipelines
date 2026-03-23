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
    anova = file.path(base_dir, "data/anova"),
    gsea_data = file.path(base_dir, "data/gsea_data"),
    volcano = file.path(base_dir, "figures/volcano"),
    ma = file.path(base_dir, "figures/ma"),
    heatmap = file.path(base_dir, "figures/heatmap"),
    imputation = file.path(base_dir, "figures/imputation"),
    gsea = file.path(base_dir, "figures/gsea"),
    pca = file.path(base_dir, "figures/pca")
  )
  
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# Build a tree string of actual output files under out_dirs$data and out_dirs$figures.
# Returns a single character string suitable for printing inside a code block.
build_output_tree <- function(out_dirs) {
  B <- "\u251C\u2500\u2500 "  # branch:  ├──
  L <- "\u2514\u2500\u2500 "  # last:    └──
  P <- "\u2502   "            # pipe:    │
  S <- "    "                 # space

  # Helper: render a named list of subdirs as tree lines
  render_dir <- function(path, prefix = "") {
    files <- sort(list.files(path, full.names = FALSE))
    if (length(files) == 0) return(character(0))
    lines <- character(0)
    for (i in seq_along(files)) {
      connector <- if (i < length(files)) B else L
      lines <- c(lines, paste0(prefix, connector, files[i]))
    }
    lines
  }

  # Top-level dirs to render, in display order
  top_dirs <- list(
    data    = out_dirs$data,
    figures = out_dirs$figures
  )

  # Second-level subdirs under data/ and figures/
  data_subdirs    <- list(anova = out_dirs$anova, de_data = out_dirs$de_data, gsea_data = out_dirs$gsea_data)
  figures_subdirs <- list(gsea       = out_dirs$gsea,
                          heatmap    = out_dirs$heatmap,
                          imputation = out_dirs$imputation,
                          ma         = out_dirs$ma,
                          pca        = out_dirs$pca,
                          volcano    = out_dirs$volcano)

  render_section <- function(subdirs, outer_prefix, last_top) {
    lines <- character(0)
    subnames <- names(subdirs)
    # filter to only dirs that exist and have files
    subnames <- subnames[sapply(subdirs[subnames], function(d)
      dir.exists(d) && length(list.files(d)) > 0)]
    for (j in seq_along(subnames)) {
      nm   <- subnames[j]
      path <- subdirs[[nm]]
      is_last_sub <- j == length(subnames)
      sub_conn  <- if (is_last_sub) L else B
      sub_pipe  <- if (last_top) S else P
      file_pipe <- if (is_last_sub) paste0(outer_prefix, sub_pipe, S) else paste0(outer_prefix, sub_pipe, P)
      lines <- c(lines, paste0(outer_prefix, sub_pipe, sub_conn, nm))
      lines <- c(lines, render_dir(path, prefix = file_pipe))
    }
    lines
  }

  all_lines <- "."
  top_names <- names(top_dirs)
  for (i in seq_along(top_names)) {
    nm      <- top_names[i]
    is_last <- i == length(top_names)
    conn    <- if (is_last) L else B
    all_lines <- c(all_lines, paste0(conn, nm))
    subdirs <- if (nm == "data") data_subdirs else figures_subdirs
    all_lines <- c(all_lines, render_section(subdirs, outer_prefix = "", last_top = is_last))
  }

  paste(all_lines, collapse = "\n")
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
                                     out_dir = ".",
                                     color1 = "#D55E00", color2 = "#0072B2") {

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

  ns_dat  <- labeled_dat %>% filter(is.na(color_tag))
  sig_dat <- labeled_dat %>% filter(!is.na(color_tag))
  max_lfc <- max(abs(labeled_dat$logFC), na.rm = TRUE)

  volcano_plot <- ggplot() +
    geom_point(data = ns_dat,  aes(x = logFC, y = -log10(P.Value)),
               color = "gray70", alpha = 0.4) +
    geom_point(data = sig_dat, aes(x = logFC, y = -log10(P.Value), color = logFC),
               alpha = 0.7) +
    scale_color_gradientn(
      colours = c(color2, "white", color1),
      limits  = c(-max_lfc, max_lfc),
      name    = paste0("higher in\n", exp_name, "\n\u2191\n\n\u2193\nhigher in\n", ctrl_name),
      guide   = guide_colorbar(direction = "vertical",
                               title.position = "right",
                               title.hjust = 0.5,
                               title.theme = element_text(angle = 90, hjust = 0.5, size = 9))
    ) +
    geom_label_repel(data = labeled_dat %>% filter(!is.na(label_display)),
                     aes(x = logFC, y = -log10(P.Value), label = label_display),
                     max.overlaps = Inf, show.legend = FALSE) +
    geom_hline(yintercept = -log10(p_thresh), col = "gray40", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")"),
         y = "-log10(P-value)",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


generate_volcano_phospho <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                     sig = "adj.P.Val", out_dir = ".",
                                     color1 = "#D55E00", color2 = "#0072B2") {

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

  ns_dat  <- labeled_dat %>% filter(is.na(color_tag))
  sig_dat <- labeled_dat %>% filter(!is.na(color_tag))
  max_lfc <- max(abs(labeled_dat$logFC), na.rm = TRUE)

  volcano_plot <- ggplot() +
    geom_point(data = ns_dat,  aes(x = logFC, y = -log10(P.Value)),
               color = "gray70", alpha = 0.4) +
    geom_point(data = sig_dat, aes(x = logFC, y = -log10(P.Value), color = logFC),
               alpha = 0.7) +
    scale_color_gradientn(
      colours = c(color2, "white", color1),
      limits  = c(-max_lfc, max_lfc),
      name    = paste0("higher in\n", exp_name, "\n\u2191\n\n\u2193\nhigher in\n", ctrl_name),
      guide   = guide_colorbar(direction = "vertical",
                               title.position = "right",
                               title.hjust = 0.5,
                               title.theme = element_text(angle = 90, hjust = 0.5, size = 9))
    ) +
    geom_label_repel(data = labeled_dat %>% filter(!is.na(label_display)),
                     aes(x = logFC, y = -log10(P.Value), label = label_display),
                     max.overlaps = Inf, show.legend = FALSE) +
    geom_hline(yintercept = -log10(p_thresh), col = "gray40", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")"),
         y = "-log10(P-value)",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


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

  ns_df  <- ma_df %>% filter(sig == "NS")
  sig_df <- ma_df %>% filter(sig != "NS")
  max_lfc <- max(abs(ma_df$logFC), na.rm = TRUE)

  p <- ggplot() +
    geom_point(data = ns_df,  aes(x = AveExpr, y = logFC),
               color = "gray60", alpha = 0.3, size = 1) +
    geom_point(data = sig_df, aes(x = AveExpr, y = logFC, color = logFC),
               alpha = 0.6, size = 1) +
    scale_color_gradientn(
      colours = c(color2, "white", color1),
      limits  = c(-max_lfc, max_lfc),
      name    = paste0("higher in\n", exp_name, "\n\u2191\n\n\u2193\nhigher in\n", ctrl_name),
      guide   = guide_colorbar(direction = "vertical",
                               title.position = "right",
                               title.hjust = 0.5,
                               title.theme = element_text(angle = 90, hjust = 0.5, size = 9))
    ) +
    geom_hline(yintercept = 0,    linetype = "dashed", color = "black") +
    geom_hline(yintercept =  lfc, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = -lfc, linetype = "dotted", color = "gray40") +
    labs(x = "Average log2 intensity (AveExpr)",
         y = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")"),
         title = create_comparison_name(exp_name, ctrl_name, "MA Plot - ")) +
    theme_bw() +
    theme(legend.position = "right")

  if (!is.null(ma_labeled) && nrow(ma_labeled) > 0) {
    p <- p + geom_label_repel(data = ma_labeled, aes(label = label_display),
                               max.overlaps = Inf, show.legend = FALSE, size = 3)
  }
  p
}


generate_ma_plot_phospho <- function(data, exp_name, ctrl_name, highlighted_ids = NULL,
                                      p_thresh = 0.05, lfc = 0.58,
                                      color1 = "#D55E00", color2 = "#0072B2") {
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

  ns_df  <- ma_df %>% filter(sig == "NS")
  sig_df <- ma_df %>% filter(sig != "NS")
  max_lfc <- max(abs(ma_df$logFC), na.rm = TRUE)

  p <- ggplot() +
    geom_point(data = ns_df,  aes(x = AveExpr, y = logFC),
               color = "gray60", alpha = 0.3, size = 1) +
    geom_point(data = sig_df, aes(x = AveExpr, y = logFC, color = logFC),
               alpha = 0.6, size = 1) +
    scale_color_gradientn(
      colours = c(color2, "white", color1),
      limits  = c(-max_lfc, max_lfc),
      name    = paste0("higher in\n", exp_name, "\n\u2191\n\n\u2193\nhigher in\n", ctrl_name),
      guide   = guide_colorbar(direction = "vertical",
                               title.position = "right",
                               title.hjust = 0.5,
                               title.theme = element_text(angle = 90, hjust = 0.5, size = 9))
    ) +
    geom_hline(yintercept = 0,    linetype = "dashed", color = "black") +
    geom_hline(yintercept =  lfc, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = -lfc, linetype = "dotted", color = "gray40") +
    labs(x = "Average log2 intensity (AveExpr)",
         y = paste0("log2 fold-change (", exp_name, " / ", ctrl_name, ")"),
         title = create_comparison_name(exp_name, ctrl_name, "MA Plot - ")) +
    theme_bw() +
    theme(legend.position = "right")

  if (!is.null(ma_labeled) && nrow(ma_labeled) > 0) {
    p <- p + geom_label_repel(data = ma_labeled, aes(label = label_display),
                               max.overlaps = Inf, show.legend = FALSE, size = 3)
  }
  p
}


generate_global_heatmap <- function(intensity_matrix, out_dirs, top_n = 1000, molecule_label = "Proteins",
                                    color1 = "#D55E00", color2 = "#0072B2") {
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
    col = colorRamp2(c(-2, 0, 2), c(color2, "white", color1)),
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
  draw(ht, column_title = paste0("Global ", molecule_label, " Heatmap"),
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()

  return(out_path)
}

# Generate four imputation QC figures and save them to figures/imputation/.
# molecule_label: "protein" or "peptide" — used in axis labels and filenames.
generate_imputation_figures <- function(intensity_raw, intensity_imputed, out_dirs,
                                        molecule_label = "protein",
                                        color1 = "#D55E00", color2 = "#0072B2") {
  fig_dir <- out_dirs$imputation

  raw_mat     <- as.matrix(intensity_raw)
  imputed_mat <- as.matrix(intensity_imputed)
  # Align raw_mat to post-filter rows (not-imputable rows may have been dropped)
  raw_mat  <- raw_mat[rownames(imputed_mat), , drop = FALSE]
  obs_mask <- !is.na(raw_mat)

  # --- 1. Observed vs Imputed histogram ---
  raw_vals <- data.frame(value = raw_mat[obs_mask],  type = "Observed")
  imp_vals <- data.frame(value = imputed_mat[!obs_mask], type = "Imputed")
  all_vals      <- rbind(raw_vals, imp_vals)
  all_vals$type <- factor(all_vals$type, levels = c("Observed", "Imputed"))

  p1 <- ggplot(all_vals, aes(x = value, fill = type)) +
    geom_histogram(alpha = 0.6, bins = 80, position = "identity") +
    scale_fill_manual(values = c("Observed" = color2, "Imputed" = color1)) +
    labs(x = "log2 Intensity", y = "Count", fill = NULL,
         title = "Observed vs. Imputed Intensity Values") +
    theme_bw() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, size = 11)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
  save_plot(p1, file.path(fig_dir, "global_imputation_histogram.png"), width = 9, height = 5)

  # --- 2. Missing values per sample ---
  miss_per_sample <- colSums(is.na(raw_mat))
  miss_df <- data.frame(
    sample  = factor(names(miss_per_sample), levels = names(miss_per_sample)),
    missing = as.integer(miss_per_sample)
  )
  p2 <- ggplot(miss_df, aes(x = sample, y = missing)) +
    geom_bar(stat = "identity", fill = color2, color = "white") +
    geom_text(aes(label = missing), vjust = -0.4, size = 3) +
    labs(x = "Sample", y = "Number of missing values",
         title = "Number of Missing Values per Sample") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 11))
  fig_width_samples <- max(6, ncol(obs_mask) * 0.6)
  save_plot(p2, file.path(fig_dir, "global_imputation_missing_per_sample.png"),
            width = fig_width_samples, height = 5)

  # --- 3. Total imputed values per molecule ---
  n_obs  <- sum(obs_mask)
  n_tot  <- length(obs_mask)
  imp_per_mol <- rowSums(!obs_mask)
  max_imp     <- ncol(obs_mask)
  imp_counts  <- table(factor(imp_per_mol, levels = 0:max_imp))
  imp_df <- data.frame(n_imputed  = as.integer(names(imp_counts)),
                       n_molecules = as.integer(imp_counts))
  imp_df$total_val <- imp_df$n_imputed * imp_df$n_molecules
  imp_df$total_val[imp_df$n_imputed == 0] <- n_obs
  imp_df$bar_type  <- ifelse(imp_df$n_imputed == 0, "Observed", "Imputed")
  imp_df$pct_total <- round(100 * imp_df$total_val / n_tot, 1)
  x_labels <- setNames(c("no-imputation", as.character(seq_len(max_imp))), 0:max_imp)

  p3 <- ggplot(imp_df, aes(x = factor(n_imputed), y = total_val, fill = bar_type)) +
    geom_bar(stat = "identity", color = "white") +
    geom_text(aes(label = paste0(pct_total, "%")), vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("Observed" = color2, "Imputed" = color1)) +
    scale_x_discrete(labels = x_labels) +
    labs(x = paste0("Number of imputed values per ", molecule_label),
         y = "Total values",
         title = paste0("Total imputed values by imputation count per ", molecule_label),
         fill = NULL) +
    theme_bw() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p3, file.path(fig_dir, "global_imputation_total_counts.png"), width = 7, height = 5)

  # --- 4. Distribution of imputed values per molecule ---
  imp_df$pct <- round(100 * imp_df$n_molecules / nrow(obs_mask), 1)

  p4 <- ggplot(imp_df, aes(x = factor(n_imputed), y = n_molecules)) +
    geom_bar(stat = "identity", fill = color2, color = "white") +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.4, size = 3) +
    scale_x_discrete(labels = x_labels) +
    labs(x = paste0("Number of imputed values per ", molecule_label),
         y = paste0("Number of ", molecule_label, "s"),
         title = paste0("Distribution of imputed values per ", molecule_label)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p4, file.path(fig_dir, "global_imputation_distribution.png"), width = 7, height = 5)

  invisible(fig_dir)
}

generate_heatmap <- function(results_df, normalized_counts, p = 0.05, lfc = 0.58,
                             exp_name, ctrl_name, fig_dir, design,
                             row_id_col = "uniprotswissprot",
                             color1 = "#D55E00", color2 = "#0072B2") {

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
  #print("WARNING JR: REMOVE THIS REPLACING OF NAN WITH SOMETHING MORE ACCURATE WHEN READY")
  #zscores[is.na(zscores) | is.nan(zscores) | is.infinite(zscores)] <- 0
  
  # plot the heatmap
  ht <- Heatmap(
    zscores,
    name = "Z-score",  # nombre de la leyenda
    col = colorRamp2(
      c(min(zscores), 0, max(zscores)), 
      c(color2, "white", color1)
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
    attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol"),
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

create_barplot <- function(gse, title, color1 = "#D55E00", color2 = "#0072B2") {
  res <- as.data.frame(gse)

  top_pos <- res %>% filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = 15)
  top_neg <- res %>% filter(NES < 0) %>% arrange(NES)        %>% slice_head(n = 15)
  plot_df <- bind_rows(top_pos, top_neg) %>%
    mutate(Description = stringr::str_wrap(Description, width = 40))

  ggplot(plot_df, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, color = "gray30", linetype = "dashed") +
    scale_fill_gradient(low = color1, high = color2, name = "adj. p-value") +
    theme_classic() +
    theme(
      axis.title.y       = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.4, linetype = 3),
      plot.title         = element_text(hjust = 0.5, size = 11)
    ) +
    labs(x = "Normalized Enrichment Score (NES)", title = title)
}

# ==========================
# centralizing function
# ==========================

run_analysis <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, ont_option = "BP", skip_gsea = FALSE, protein_metadata = NULL, color1 = "#D55E00", color2 = "#0072B2") {
  tryCatch({
    print(paste("\nStarting analysis for comparison:", comparison$name))

    print('Running limma')
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
      print("Imputation category distribution:")
      print(table(annotated_results$imputation_category))
    }

    output_file <- create_file_path(out_dirs$de_data, comparison$name, "_limma")
    write.csv(annotated_results, output_file)

    print('Generating the volcano plot')
    volcano_result <- generate_volcano_protein(annotated_results, comparison$exp, comparison$ctrl, color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    print('Generating the MA plot')
    ma_plot <- generate_ma_plot_protein(annotated_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids,
                                        color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    print('Generating the heatmap')
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, intensity_matrix,
                     exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                     fig_dir = out_dirs$heatmap, design = limma_params$design,
                     color1 = color1, color2 = color2)
    draw(ht, column_title = paste0(comparison$exp, " vs ", comparison$ctrl, " — Differentially Abundant Proteins"),
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()

    if (skip_gsea) {
      message("Skipping GSEA (--skip-gsea)")
      gse <- NULL
    } else {
      print('Running GSEA')
      gse <- NULL # process_gsea(annotated_results, ont_option = ont_option)
      #gse <- process_gsea(annotated_results, ont_option = ont_option)

      if(!is.null(gse)) {
        print('GSE has data')
        write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, "go_analysis_", comparison$name))

        print('create dotplot')
        gsea_plot <- create_barplot(gse, create_comparison_name(comparison$exp, comparison$ctrl, "GSEA "), color1 = color1, color2 = color2)

        print('save plot')
        save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                  width = 10, height = 12)
      }
    }

    return(list(limma = annotated_results, gsea = gse, highlighted_ids = volcano_result$highlighted_ids))
    
    
  }, error = function(e) {
    message("Error in run_analysis: ", e$message)
    return(NULL)
  })
}

run_analysis_phospho <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, peptide_metadata = NULL, ont_option = "BP", skip_gsea = FALSE, color1 = "#D55E00", color2 = "#0072B2") {
  tryCatch({
    print(paste("Starting analysis for comparison:", comparison$name))

    print('Running limma')
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
        message("BioMart annotation failed for phospho: ", e$message)
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
      print("Imputation category distribution:")
      print(table(limma_results$imputation_category))
    }

    output_file <- create_file_path(out_dirs$de_data, comparison$name, "_limma")
    write.csv(limma_results, output_file)

    print('Generating the volcano plot')
    volcano_result <- generate_volcano_phospho(limma_results, comparison$exp, comparison$ctrl, color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)

    print('Generating the MA plot')
    ma_plot <- generate_ma_plot_phospho(limma_results, comparison$exp, comparison$ctrl,
                                        highlighted_ids = volcano_result$highlighted_ids,
                                        color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)

    print('Generating the heatmap')
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 2400, height = 3200, res = 300)
    ht <- generate_heatmap(limma_results, normalized_counts,
                           exp_name = comparison$exp, ctrl_name = comparison$ctrl,
                           fig_dir = out_dirs$heatmap, design = limma_params$design,
                           row_id_col = "peptide_id",
                           color1 = color1, color2 = color2)
    draw(ht, column_title = paste0(comparison$exp, " vs ", comparison$ctrl, " — Differentially Abundant Phosphopeptides"),
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()

    if (skip_gsea) {
      message("Skipping GSEA (--skip-gsea)")
      gse <- NULL
      protein_counts <- NULL
    } else {
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
          write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, comparison$name, "_go_analysis.csv"))

          print('create dotplot')
          gsea_plot <- create_barplot(gse, create_comparison_name(comparison$ctrl, comparison$exp, "GSEA "), color1 = color1, color2 = color2)

          print('save plot')
          save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                    width = 10, height = 12)
        }, error = function(e) {
          message("Error saving phospho GSEA outputs: ", e$message)
        })
      }
    }

    return(list(limma = limma_results, gsea = gse, protein_counts = protein_counts, highlighted_ids = volcano_result$highlighted_ids))
    
    
  }, error = function(e) {
    message("Error in run_analysis: ", e$message)
    return(NULL)
  })
}
