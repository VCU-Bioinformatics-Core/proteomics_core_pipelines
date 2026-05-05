# ==========================
# Logging setup
# ==========================
setup_logging <- function(outdir, run_id = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  log_name <- if (!is.null(run_id) && nzchar(run_id)) {
    paste0(run_id, "_pipeline.log")
  } else {
    "pipeline.log"
  }
  log_path <- file.path(outdir, log_name)
  flog.appender(appender.tee(log_path))
  flog.threshold(INFO)
  invisible(log_path)
}

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
    flog.warn("Failed to export plotly to HTML '%s': %s", file_path, e$message)
  })
}

setup_directories <- function(base_dir) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
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
    flog.info("Performing limma analysis: %s vs %s", exp, ctrl)
    
    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(limma_params$design))
    
    # Fit linear model
    fit <- lmFit(limma_params$E, limma_params$design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
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
                    "\u2190 higher in ", ctrl_name, "          higher in ", exp_name, " \u2192"),
         y = "-log10(P-value)",
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(list(plot = volcano_plot, highlighted_ids = c(top_10_genes_up, top_10_genes_dn)))
}


generate_volcano_ptm <- function(data, exp_name, ctrl_name, p_thresh = 0.05, lfc = 0.58,
                                     sig = "adj.P.Val", out_dir = ".",
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
      eval(as.symbol(sig)) < p_thresh & logFC >  lfc ~ "Up-Regulated",
      eval(as.symbol(sig)) < p_thresh & logFC < -lfc ~ "Down-Regulated",
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
                    "\u2190 higher in ", ctrl_name, "          higher in ", exp_name, " \u2192"),
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
                    "\n\u2190 higher in ", ctrl_name, "          higher in ", exp_name, " \u2192"),
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
                    "\u2191 higher in ", exp_name, "\n\u2193 higher in ", ctrl_name),
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


generate_global_heatmap <- function(intensity_matrix, out_dirs, top_n = 1000, molecule_label = "Proteins",
                                    heatmap_norm = "zscore",
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

  if (heatmap_norm == "zscore") {
    plot_mat <- t(scale(t(mat)))
    plot_mat[is.nan(plot_mat) | is.infinite(plot_mat)] <- 0
    legend_name <- "Z-score"
    col_scale <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  } else {
    plot_mat <- mat
    legend_name <- "Intensity"
    mid <- median(plot_mat, na.rm = TRUE)
    col_scale <- colorRamp2(c(min(plot_mat, na.rm = TRUE), mid, max(plot_mat, na.rm = TRUE)),
                            c("#2166AC", "white", "#B2182B"))
  }

  out_path <- file.path(out_dirs$heatmap, "global_heatmap.png")
  png(out_path, width = 2400, height = 3200, res = 300)
  ht <- Heatmap(
    plot_mat,
    name = legend_name,
    col = col_scale,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_title = paste0(molecule_label, " (top ", n_select, " by CV)"),
    column_title = "Samples",
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_rot = 45,
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_gp = grid::gpar(fontsize = 12)
  )
  draw(ht, column_title = paste0("Global ", molecule_label, " Heatmap"),
       column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
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
  imp_vals <- data.frame(value = imputed_mat[!obs_mask], type = rep("Imputed", sum(!obs_mask)))
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
  imp_df$bar_type  <- factor(ifelse(imp_df$n_imputed == 0, "Observed", "Imputed"),
                             levels = c("Observed", "Imputed"))
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
                             heatmap_norm = "zscore",
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

  if (heatmap_norm == "zscore") {
    plot_mat <- t(scale(t(values)))
    legend_name <- "Z-score"
    col_scale <- colorRamp2(c(min(plot_mat, na.rm = TRUE), 0, max(plot_mat, na.rm = TRUE)),
                            c("#2166AC", "white", "#B2182B"))
  } else {
    plot_mat <- values
    legend_name <- "Intensity"
    mid <- median(plot_mat, na.rm = TRUE)
    col_scale <- colorRamp2(c(min(plot_mat, na.rm = TRUE), mid, max(plot_mat, na.rm = TRUE)),
                            c("#2166AC", "white", "#B2182B"))
  }

  # plot the heatmap
  ht <- Heatmap(
    plot_mat,
    name = legend_name,
    col = col_scale,
    cluster_rows = TRUE,     # dendrograma de filas
    cluster_columns = TRUE,  # dendrograma de columnas
    show_row_names = FALSE,  # equivalente a labRow = NA
    show_column_names = TRUE,
    row_title = "Features",
    column_title = "Samples",
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_rot = 45,
    column_title_gp = grid::gpar(fontsize = 12),
    row_title_gp = grid::gpar(fontsize = 12)
  )
  
  return(ht)
  
}

# ==========================
# gsea analysis functions
# ==========================

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
    flog.warn("No Ensembl IDs found for aggregated PTM UniProt IDs — GSEA skipped"
    )
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
      flog.warn("Not enough significant genes for GSEA (n < 2)")
      return(NULL)
    }
    
    # create a gene list where the names are ensembl genes
    # and the values are logFC
    gene_list <- sig_genes$logFC
    names(gene_list) <- sig_genes$ensembl_gene_id
    gene_list <- gene_list[!is.na(names(gene_list))]
    gene_list <- gene_list[!duplicated(names(gene_list))]
    gene_list <- gene_list[abs(gene_list) >= 0.1]  # drop near-zero logFC genes

    if (length(gene_list) < 2) {
      flog.warn("Not enough genes with Ensembl IDs for GSEA (n < 2)")
      return(NULL)
    }

    gse_result <- gseGO(geneList = gene_list,
                        ont = ont_option,
                        minGSSize = 5,
                        maxGSSize = 500,
                        keyType = "ENSEMBL",
                        scoreType = "std",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        OrgDb = annotation_db)

    if (is.null(gse_result) || nrow(gse_result@result) == 0) {
      flog.warn("No enriched GO terms found in GSEA result")
      return(NULL)
    }

    setReadable(gse_result, OrgDb = annotation_db, keyType = "ENSEMBL")
  }, error = function(e) {
    flog.error("GSEA processing failed: %s", e$message)
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

run_analysis_ptm <- function(comparison, limma_params, normalized_counts, out_dirs, intensity_matrix_raw = NULL, peptide_metadata = NULL, ont_option = "BP", skip_gsea = FALSE, heatmap_norm = "zscore", color1 = "#D55E00", color2 = "#0072B2") {
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
