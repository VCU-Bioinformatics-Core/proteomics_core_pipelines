# ==========================
# Logging setup
# ==========================

#' @title Set Up Pipeline Logging
#' @param outdir Character. Path to the output directory where the log file will be written.
#'   Created recursively if it does not exist.
#' @param run_id Character or NULL. Optional run identifier prepended to the log filename
#'   (e.g. \code{"run01"} produces \code{run01_pipeline.log}). If NULL or empty, defaults
#'   to \code{pipeline.log}.
#' @return Invisibly returns the full path to the log file as a character string.
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

#' @title Build a Standardised Output File Path
#' @param base_dir Character. Base directory for the file.
#' @param prefix Character. String prepended to \code{name} in the filename.
#' @param name Character. Core part of the filename.
#' @param extension Character. File extension including the leading dot. Default \code{".csv"}.
#' @return Character string with the full file path.
create_file_path <- function(base_dir, prefix, name, extension = ".csv") {
  file.path(base_dir, paste0(prefix, name, extension))
}

#' @title Build a Comparison Label
#' @param exp Character. Name of the experimental group.
#' @param ctrl Character. Name of the control group.
#' @param prefix Character. Optional string prepended to the label. Default \code{""}.
#' @return Character string of the form \code{"<prefix><exp> vs. <ctrl>"}.
create_comparison_name <- function(exp, ctrl, prefix = "") {
  paste0(prefix, exp, " vs. ", ctrl)
}

#' @title Save a ggplot to Disk
#' @param plot A \code{ggplot} object to save.
#' @param filename Character. Full output file path including extension.
#' @param width Numeric. Plot width in inches. Default \code{10}.
#' @param height Numeric. Plot height in inches. Default \code{8}.
#' @param dpi Numeric. Resolution in dots per inch. Default \code{300}.
#' @return Invisibly returns \code{filename} (via \code{ggsave}).
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  ggsave(plot, filename = filename, width = width, height = height, dpi = dpi, bg = "white")
}

#' @title Export a Plotly Object to a Self-Contained HTML File
#' @param plotly_obj A \code{plotly} object to export.
#' @param file_path Character. Destination HTML file path.
#' @return Invisibly returns \code{NULL}. Logs a warning if the export fails.
export_plotly_to_html <- function(plotly_obj, file_path) {
  tryCatch({
    if (!inherits(plotly_obj, "plotly")) stop("Invalid plotly object")
    htmlwidgets::saveWidget(plotly_obj, file_path, selfcontained = TRUE)
  }, error = function(e) {
    flog.warn("Failed to export plotly to HTML '%s': %s", file_path, e$message)
  })
}

#' @title Create Standard Pipeline Output Directory Structure
#' @param base_dir Character. Root output directory. Will be normalised via
#'   \code{normalizePath(..., mustWork = FALSE)}.
#' @return Named list of character paths for each subdirectory:
#'   \code{data}, \code{figures}, \code{de_data}, \code{anova}, \code{gsea_data},
#'   \code{volcano}, \code{ma}, \code{heatmap}, \code{imputation}, \code{gsea}, \code{pca}.
#'   All directories are created (recursively) as a side effect.
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

#' @title Build a Tree-Style String of Pipeline Output Files
#' @details Renders the actual files present under \code{out_dirs$data} and
#'   \code{out_dirs$figures} (and their standard subdirectories) as a Unicode
#'   tree string suitable for embedding in a report code block. Subdirectories
#'   with no files are silently omitted.
#' @param out_dirs Named list of directory paths as returned by
#'   \code{\link{setup_directories}}.
#' @return Single character string with newline-separated tree lines rooted at \code{"."}.
build_output_tree <- function(out_dirs) {
  B <- "├── "  # branch:  ├──
  L <- "└── "  # last:    └──
  P <- "│   "            # pipe:    │
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
# Shared visualisation functions
# ==========================

#' @title Generate a Global Heatmap of All Detected Molecules
#' @details Filters to the top \code{top_n} most variable rows by coefficient of variation,
#'   applies z-score or raw intensity normalisation, clusters rows and columns, and writes
#'   a PNG to \code{out_dirs$heatmap/global_heatmap.png}.
#' @param intensity_matrix Numeric matrix or data frame (features x samples) of log2
#'   intensities. Rows with non-finite values are dropped before clustering.
#' @param out_dirs Named list of output directories as returned by
#'   \code{\link{setup_directories}}.
#' @param top_n Integer. Maximum number of rows (by CV) to include. Default \code{1000}.
#' @param molecule_label Character. Label used in plot titles and axis (e.g. \code{"Proteins"},
#'   \code{"Peptides"}). Default \code{"Proteins"}.
#' @param heatmap_norm Character. Normalisation method: \code{"zscore"} (default) scales each
#'   row to mean 0 / sd 1; any other value uses raw intensities.
#' @param color1 Character. Reserved colour for consistency with other plot helpers.
#'   Default \code{"#D55E00"}.
#' @param color2 Character. Reserved colour for consistency with other plot helpers.
#'   Default \code{"#0072B2"}.
#' @return Character string giving the full path to the saved PNG file.
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

#' @title Generate Imputation QC Figures
#' @details Produces four PNG figures saved to \code{out_dirs$imputation}:
#'   \enumerate{
#'     \item Observed vs. imputed intensity histogram.
#'     \item Number of missing values per sample (bar chart).
#'     \item Total observed/imputed values by imputation count per molecule (stacked bar).
#'     \item Distribution of imputation counts per molecule (bar chart with percentages).
#'   }
#' @param intensity_raw Numeric matrix or data frame (features x samples) of raw log2
#'   intensities, with \code{NA} where values are missing.
#' @param intensity_imputed Numeric matrix or data frame (features x samples) of imputed
#'   log2 intensities. Row names must be a subset of \code{intensity_raw} row names.
#' @param out_dirs Named list of output directories as returned by
#'   \code{\link{setup_directories}}.
#' @param molecule_label Character. Singular label for the molecule type used in axis and
#'   title text (e.g. \code{"protein"}, \code{"peptide"}). Default \code{"protein"}.
#' @param color1 Character. Hex colour for "Imputed" bars/fills. Default \code{"#D55E00"}.
#' @param color2 Character. Hex colour for "Observed" bars/fills. Default \code{"#0072B2"}.
#' @return Invisibly returns the path to the imputation figure directory
#'   (\code{out_dirs$imputation}).
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

#' @title Generate a Per-Comparison Differential Abundance Heatmap
#' @details Filters \code{diff_df} to significant features (\code{adj.P.Val < p} and
#'   \code{|logFC| >= lfc}), subsets \code{intensity_df} to the relevant samples, adds
#'   a small jitter to break ties, and returns a \code{ComplexHeatmap::Heatmap} object
#'   (not yet drawn or saved).
#' @param diff_df Data frame of limma results (output of \code{\link{perform_limma_analysis}})
#'   containing at minimum \code{adj.P.Val}, \code{logFC}, and the column named by
#'   \code{row_id_col}.
#' @param intensity_df Numeric data frame or matrix (features x samples) of log2 intensities.
#'   Row names must match values in \code{diff_df[[row_id_col]]}.
#' @param p Numeric. Adjusted p-value threshold for significance. Default \code{0.05}.
#' @param lfc Numeric. Minimum absolute log2 fold-change threshold. Default \code{0.58}
#'   (~1.5-fold).
#' @param exp_name Character. Column name in \code{design} identifying experimental samples.
#' @param ctrl_name Character. Column name in \code{design} identifying control samples.
#' @param fig_dir Character. Directory for figure output (saving is the caller's
#'   responsibility).
#' @param design Model matrix (samples x groups) used to identify which samples belong to
#'   each group.
#' @param row_id_col Character. Column in \code{diff_df} whose values index rows of
#'   \code{intensity_df}. Default \code{"uniprotswissprot"}.
#' @param heatmap_norm Character. \code{"zscore"} (default) or any other value for raw
#'   intensity.
#' @param color1 Character. Reserved colour parameter for consistency. Default
#'   \code{"#D55E00"}.
#' @param color2 Character. Reserved colour parameter for consistency. Default
#'   \code{"#0072B2"}.
#' @return A \code{ComplexHeatmap::Heatmap} object ready to be drawn with \code{draw()}.
generate_heatmap <- function(diff_df, intensity_df, p = 0.05, lfc = 0.58,
                             exp_name, ctrl_name, fig_dir, design,
                             row_id_col = "uniprotswissprot",
                             heatmap_norm = "zscore",
                             color1 = "#D55E00", color2 = "#0072B2") {

  # filter the data for certain pvalues and lfc
  filtered_data <- diff_df %>%
    dplyr::filter((adj.P.Val < p & abs(logFC) >= lfc))
  flog.info("generate_heatmap: %d rows after significance filter (p=%.2f, lfc=%.2f)", nrow(filtered_data), p, lfc)
  flog.info("generate_heatmap: row_id_col='%s', values: %s", row_id_col, paste(head(filtered_data[[row_id_col]], 3), collapse=", "))

  # restrict columns to only the two groups being compared
  all_samples  <- colnames(intensity_df)
  exp_samples  <- all_samples[design[, exp_name]  == 1]
  ctrl_samples <- all_samples[design[, ctrl_name] == 1]
  comparison_samples <- c(exp_samples, ctrl_samples)
  flog.info("generate_heatmap: exp_samples=%s", paste(exp_samples, collapse=", "))
  flog.info("generate_heatmap: ctrl_samples=%s", paste(ctrl_samples, collapse=", "))

  # add jitter
  values <- intensity_df[filtered_data[[row_id_col]], comparison_samples] %>% as.matrix()
  values <- values %>% jitter(factor = 1, amount = 0.00001)

  flog.info("generate_heatmap: values matrix dim=%s, NA count=%d, Inf count=%d",
            paste(dim(values), collapse="x"),
            sum(is.na(values)),
            sum(is.infinite(values)))

  if (heatmap_norm == "zscore") {
    plot_mat <- t(scale(t(values)))
    flog.info("generate_heatmap: plot_mat after zscore dim=%s, NA count=%d, Inf count=%d",
              paste(dim(plot_mat), collapse="x"),
              sum(is.na(plot_mat)),
              sum(is.infinite(plot_mat)))
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

  ht <- Heatmap(
    plot_mat,
    name = legend_name,
    col = col_scale,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
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

#' @title Create a GSEA NES Bar Plot
#' @details Selects up to 15 top positively enriched and 15 top negatively enriched GO
#'   terms by NES, wraps long descriptions, and plots a horizontal bar chart coloured by
#'   adjusted p-value.
#' @param gse A \code{gseaResult} object (output of \code{\link{process_gsea}}).
#' @param title Character. Plot title.
#' @param color1 Character. Colour for the low end of the p-value gradient (most
#'   significant). Default \code{"#D55E00"}.
#' @param color2 Character. Colour for the high end (least significant). Default
#'   \code{"#0072B2"}.
#' @return A \code{ggplot} object. Not saved to disk — pass to \code{\link{save_plot}}.
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
