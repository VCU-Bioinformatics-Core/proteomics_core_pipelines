# ==========================
# Main PTM pipeline function
# ==========================
# run_ptm_pipeline()      — PTM peptide-level differential abundance
#
# These functions contain the full analysis logic previously in
# workflow/de.regular.R and workflow/de.ptm.R. The thin CLI wrappers
# in inst/scripts/ parse command-line arguments and call these functions.

#' Run the PTM (post-translational modification) differential abundance pipeline
#'
#' Reads PTM-level and protein-level intensity data, performs imputation,
#' limma-based pairwise differential analysis, PCA, optional ANOVA, GSEA,
#' and writes an HTML report plus all intermediate RDS/figure outputs.
#'
#' @param run_id Character. Unique identifier for this analysis run;
#'   used in output file names and the HTML report title.
#' @param ptm_input_file Character. Path to the PTM-level intensity matrix
#'   (tab-separated, rows = PTM sites, cols = samples).
#' @param samplesheet_file Character. Path to the sample sheet CSV describing
#'   sample groups and pairwise comparisons.
#' @param protein_input_file Character. Path to the protein-level intensity
#'   matrix used for protein-level normalisation context.
#' @param out_dir Character. Directory where all outputs will be written.
#'   Created if it does not already exist.
#' @param genome Character. Organism genome to use for annotation
#'   (\code{"mouse"} or \code{"human"}). Default \code{"mouse"}.
#' @param imputation_method Character. Imputation strategy. One of
#'   \code{"DEP-MinProb"}, \code{"DEP-BPCA"}, \code{"none"}, or a custom
#'   method name. Default \code{"DEP-MinProb"}.
#' @param imputation_q Numeric. Quantile threshold passed to DEP imputation
#'   methods. Default \code{0.01}.
#' @param imputation_seed Integer. Random seed for reproducible imputation.
#'   Default \code{42}.
#' @param heatmap_top_n Integer. Number of top variable PTM sites to include
#'   in the global heatmap. Default \code{1000}.
#' @param heatmap_norm Character. Normalisation applied before heatmap
#'   (\code{"zscore"} or \code{"none"}). Default \code{"zscore"}.
#' @param gsea_ont Character. Gene Ontology domain for GSEA
#'   (\code{"BP"}, \code{"MF"}, or \code{"CC"}). Default \code{"BP"}.
#' @param skip_gsea Logical. If \code{TRUE}, skip GSEA entirely.
#'   Default \code{FALSE}.
#' @param skip_anova Logical. If \code{TRUE}, skip one-way ANOVA across all
#'   groups. Default \code{FALSE}.
#' @param group_color1 Character. Hex colour for the first group in plots.
#'   Default \code{"#D55E00"}.
#' @param group_color2 Character. Hex colour for the second group in plots.
#'   Default \code{"#0072B2"}.
#'
#' @return Invisibly returns the path to the generated HTML report.
#'
#' @examples
#' \dontrun{
#' run_ptm_pipeline(
#'   run_id             = "exp_001",
#'   ptm_input_file     = "data/ptm_intensities.tsv",
#'   samplesheet_file   = "data/samplesheet.csv",
#'   protein_input_file = "data/protein_intensities.tsv",
#'   out_dir            = "results/exp_001"
#' )
#' }
run_ptm_pipeline <- function(
  run_id,
  ptm_input_file,
  samplesheet_file,
  protein_input_file,
  out_dir,
  genome            = "mouse",
  imputation_method = "DEP-MinProb",
  imputation_q      = 0.01,
  imputation_seed   = 42,
  heatmap_top_n     = 1000,
  heatmap_norm      = "zscore",
  gsea_ont          = "BP",
  skip_gsea         = FALSE,
  skip_anova        = FALSE,
  group_color1      = "#D55E00",
  group_color2      = "#0072B2"
) {
  setup_logging(out_dir, run_id = run_id)
  flog.info("PTM pipeline started: runID=%s, genome=%s, outdir=%s", run_id, genome, out_dir)

  # ==========================
  # Load annotation DB
  # ==========================
  if (genome == "human") {
    if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
    annotation_db <- org.Hs.eg.db
    ensembl <- tryCatch(
      useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl"),
      error = function(e) {
        flog.warn("BioMart connection failed for human: %s — GSEA will be skipped", e$message)
        NULL
      }
    )
  } else if (genome == "mouse") {
    if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
    annotation_db <- org.Mm.eg.db
    ensembl <- tryCatch(
      useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl"),
      error = function(e) {
        flog.warn("BioMart connection failed for mouse: %s — GSEA will be skipped", e$message)
        NULL
      }
    )
  } else {
    flog.fatal("Invalid genome '%s'. Use 'mouse' or 'human'", genome)
    stop("Invalid genome specified. Use 'mouse' or 'human'")
  }

  # ==========================
  # Read and prepare PTM intensity data
  # ==========================

  out_dirs <- setup_directories(out_dir)
  full_peptide_levels <- data.frame(
    read_csv(ptm_input_file, col_names = TRUE, na = c("", "NA", "Filtered"))
  )

  log_peptide_count <- function(df, prefix) flog.info("%s: %d", prefix, nrow(df))

  n_peptides_total <- nrow(full_peptide_levels)
  log_peptide_count(full_peptide_levels, "Peptides - Total Amount")

  full_peptide_levels <- full_peptide_levels %>%
    filter(!grepl("cRAP[0-9]+", PG.ProteinAccessions))
  n_peptides_no_crap <- nrow(full_peptide_levels)
  log_peptide_count(full_peptide_levels, "Peptides - without cRAP")

  full_peptide_levels <- full_peptide_levels %>%
    filter(grepl("Phospho \\(STY\\)", EG.PrecursorId))
  n_peptides_ptm <- nrow(full_peptide_levels)
  log_peptide_count(full_peptide_levels, "Peptides - only PTM")

  rownames(full_peptide_levels) <- paste0(
    full_peptide_levels$PG.Genes, " -- ", full_peptide_levels$EG.PrecursorId
  )

  non_numeric_cols <- c("PG.ProteinAccessions", "PG.Genes", "PG.UniProtIds",
                        "PG.FASTAName", "EG.PrecursorId")
  counts <- full_peptide_levels %>% dplyr::select(-any_of(non_numeric_cols))

  intensity_matrix     <- log2(counts + 1)
  intensity_matrix_raw <- intensity_matrix

  peptide_metadata <- full_peptide_levels %>%
    transmute(peptide_id = rownames(full_peptide_levels), PG.UniProtIds)

  # ==========================
  # Read and prepare the comparisons
  # ==========================
  comparisons_raw  <- read.csv(samplesheet_file)
  comparisons_cols <- colnames(comparisons_raw)[3:ncol(comparisons_raw)]
  comparisons <- list()
  for (comparisons_name in comparisons_cols) {
    comparisons_vector <- comparisons_raw[[comparisons_name]]
    exp_group  <- unique(comparisons_raw$GroupID[comparisons_vector == 1 & !is.na(comparisons_vector)])
    ctrl_group <- unique(comparisons_raw$GroupID[comparisons_vector == 0 & !is.na(comparisons_vector)])
    if (length(exp_group) > 0 && length(ctrl_group) > 0) {
      comparisons[[length(comparisons) + 1]] <- list(
        name = comparisons_name, exp = exp_group, ctrl = ctrl_group
      )
    }
  }

  # ==========================
  # Imputation
  # ==========================
  DEP_METHODS <- c("MinProb", "knn", "bpca", "QRILC", "man")
  n_peptides_not_imputable <- 0L

  if (imputation_method == "none") {
    flog.info("Imputation skipped (--imputation none)")

  } else if (startsWith(imputation_method, "DEP-")) {
    dep_fun <- sub("^DEP-", "", imputation_method)
    if (!dep_fun %in% DEP_METHODS) {
      stop(glue("Unknown DEP method '{dep_fun}'. Valid DEP methods: {paste(DEP_METHODS, collapse=', ')}"))
    }
    set.seed(imputation_seed)
    se <- SummarizedExperiment::SummarizedExperiment(
      assays  = list(intensity = as.matrix(intensity_matrix)),
      colData = S4Vectors::DataFrame(label = comparisons_raw$SampleID, condition = comparisons_raw$GroupID),
      rowData = S4Vectors::DataFrame(name  = rownames(intensity_matrix), ID = rownames(intensity_matrix))
    )
    se_imputed <- if (dep_fun %in% c("MinProb", "QRILC")) {
      DEP::impute(se, fun = dep_fun, q = imputation_q)
    } else {
      DEP::impute(se, fun = dep_fun)
    }
    intensity_matrix <- as.data.frame(SummarizedExperiment::assay(se_imputed))
    flog.info("Imputation applied: %s (q=%s)", imputation_method, imputation_q)

  } else {
    fn_name <- paste0("impute_", imputation_method)
    if (!exists(fn_name, mode = "function")) {
      flog.fatal(
        "Unknown imputation method '%s'. For DEP methods prefix with 'DEP-'. For custom methods define %s() in R/custom_imputation.R.",
        imputation_method, fn_name
      )
      stop(glue(
        "Unknown imputation method '{imputation_method}'. ",
        "For DEP methods prefix with 'DEP-' (e.g. 'DEP-MinProb'). ",
        "For custom methods, define `{fn_name}()` in R/custom_imputation.R."
      ))
    }
    set.seed(imputation_seed)
    fn <- get(fn_name)
    groups_vec <- setNames(comparisons_raw$GroupID, colnames(intensity_matrix))
    imputed_mat <- fn(as.matrix(intensity_matrix), groups = groups_vec)
    n_peptides_not_imputable <- attr(imputed_mat, "n_not_imputable") %||% 0L
    intensity_matrix <- as.data.frame(imputed_mat)
    flog.info("Custom imputation applied: %s", imputation_method)
    if (n_peptides_not_imputable > 0)
      flog.warn("Discarded %d not-imputable peptide(s): all groups had <=1 valid value",
                n_peptides_not_imputable)
  }




  # ==========================
  # Read and prepare protein data
  # ==========================
  protein_intensities <- data.frame(
    read_csv(protein_input_file, col_names = TRUE, na = c("", "NA", "Filtered"))
  )










    








  # ==========================
  # Differential abundance analysis loop for Limma design matrix + analysis loop
  # ==========================
  # Currently relies on limma but I will introduce a msstatsPTM correction to run_analysis_ptm
  # DEVELOP BOOKMARK

  # Limma design matrix + analysis loop
  sample_info <- data.frame(sample = comparisons_raw$SampleID, condition = comparisons_raw$GroupID)
  design      <- model.matrix(~0 + condition, data = sample_info)
  colnames(design) <- levels(factor(sample_info$condition))
  limma_params <- list(E = intensity_matrix, design = design)

  results <- vector("list", length(comparisons))
  for (i in seq_along(comparisons)) {
    flog.info("=== Analysis loop iteration %d of %d ===", i, length(comparisons))
    curr_result <- run_analysis_ptm(
      comparisons[[i]], limma_params, intensity_matrix, out_dirs, intensity_matrix_raw,
      peptide_metadata, ont_option = gsea_ont, skip_gsea = skip_gsea,
      heatmap_norm = heatmap_norm, color1 = group_color1, color2 = group_color2
    )
    if (!is.null(curr_result)) results[[i]] <- curr_result
  }








  # ==========================
  # PCA
  # ==========================
  flog.info("Running Principal Component Analysis")
  input_pca_matrix <- as.matrix(intensity_matrix)
  finite_rows      <- apply(input_pca_matrix, 1, function(x) all(is.finite(x)))
  input_pca_matrix <- t(input_pca_matrix[finite_rows, ])
  pca_results  <- prcomp(input_pca_matrix, rank. = 3)
  pc_scores    <- pca_results$x
  pca_var_pct  <- round((pca_results$sdev^2) / sum(pca_results$sdev^2) * 100, 1)
  two_comps_pca_df <- data.frame(
    Sample = rownames(pc_scores), X = pc_scores[, 1],
    Y = pc_scores[, 2], Group = sample_info$condition
  )
  pca_plot <- ggplot(data = two_comps_pca_df, aes(x = X, y = Y, label = Sample, color = Group)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = as.character(Sample)), show.legend = FALSE, size = 2.5) +
    xlab(paste0("PC1 - ", pca_var_pct[1], "%")) +
    ylab(paste0("PC2 - ", pca_var_pct[2], "%")) +
    labs(title = "Principal Component Analysis") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
  ggsave(pca_plot, filename = str_c(out_dirs$pca, "/global_pca_plot.png"),
         width = 7, height = 5, dpi = 300)

  # ==========================
  # One-Way ANOVA
  # ==========================
  flog.info("Running one-way ANOVA")
  n_anova_groups <- length(unique(sample_info$condition))
  if (skip_anova) {
    flog.info("Skipping ANOVA (--skip-anova flag set)")
    n_anova_sig <- n_anova_total <- NULL
  } else if (n_anova_groups > 2) {
    group       <- factor(sample_info$condition)
    mat         <- as.matrix(intensity_matrix)
    anova_pvals <- apply(mat, 1, function(x) {
      tryCatch(summary(aov(x ~ group))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    })
    anova_padj    <- p.adjust(anova_pvals, method = "BH")
    n_anova_sig   <- sum(anova_padj < 0.05, na.rm = TRUE)
    n_anova_total <- nrow(mat)
    anova_df <- data.frame(
      PTM_peptide = rownames(mat),
      P_Value        = signif(anova_pvals, 3),
      Adj_P_Value    = signif(anova_padj, 3),
      stringsAsFactors = FALSE
    ) %>% dplyr::arrange(Adj_P_Value)
    write.csv(anova_df, file.path(out_dirs$anova, "global_anova.csv"), row.names = FALSE)
  } else {
    n_anova_sig <- n_anova_total <- NULL
  }
  anova_summary <- list(n_groups = n_anova_groups, n_sig = n_anova_sig,
                        n_total = n_anova_total, skipped = skip_anova)

  # ==========================
  # Global heatmap + imputation figures
  # ==========================
  flog.info("Generating global heatmap")
  generate_global_heatmap(intensity_matrix, out_dirs, top_n = heatmap_top_n,
                          molecule_label = "PTM Peptides",
                          heatmap_norm = heatmap_norm, color1 = group_color1, color2 = group_color2)

  flog.info("Generating imputation figures")
  generate_imputation_figures(intensity_matrix_raw, intensity_matrix, out_dirs,
                              molecule_label = "peptide",
                              color1 = group_color1, color2 = group_color2)

  # ==========================
  # Master query table
  # ==========================
  flog.info("Building master query table")
  all_rows <- dplyr::bind_rows(lapply(seq_along(comparisons), function(i) {
    res <- results[[i]]$limma
    if (is.null(res)) return(NULL)
    res$comparison_name <- comparisons[[i]]$name
    res$is_sig <- !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.58
    res[, intersect(c("peptide_id", "hgnc_symbol", "uniprot_id", "comparison_name", "is_sig"), colnames(res))]
  }))
  query_df <- all_rows %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(c("peptide_id", "hgnc_symbol", "uniprot_id")))) %>%
    dplyr::summarise(
      Significant_Comparisons = {
        sig_comps <- comparison_name[is_sig]
        if (length(sig_comps) == 0) "Not Significant" else paste(sig_comps, collapse = "; ")
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(Significant_Comparisons == "Not Significant", peptide_id)
  write.csv(query_df, file.path(out_dirs$de_data, "master_query_table.csv"), row.names = FALSE)

  # ==========================
  # Save RDS + generate report
  # ==========================
  imputation_params <- list(method = imputation_method, q = imputation_q)
  peptide_counts    <- list(total = n_peptides_total, no_crap = n_peptides_no_crap,
                            ptm = n_peptides_ptm, not_imputable = n_peptides_not_imputable)
  analysis_params   <- list(genome = genome, gsea_ont = gsea_ont, skip_gsea = skip_gsea,
                            heatmap_top_n = heatmap_top_n, heatmap_norm = heatmap_norm,
                            color1 = group_color1, color2 = group_color2)
  rds      <- list(results, comparisons, out_dirs, pca_plot, intensity_matrix_raw,
                   intensity_matrix, imputation_params, sample_info, peptide_counts,
                   analysis_params, anova_summary)
  rds_path <- file.path(out_dir, "analysis_results.rds")
  flog.info("Saving analysis RDS to %s", rds_path)
  saveRDS(rds, rds_path)

  generate_report_ptm(rds_path, output_dir = out_dir)
  flog.info("Pipeline complete: runID=%s", run_id)
  invisible(rds_path)
}
