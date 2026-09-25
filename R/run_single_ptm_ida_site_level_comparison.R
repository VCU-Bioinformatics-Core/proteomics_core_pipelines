# ==============================================================
# PTM-level orchestration for PTM and protein level integrations
# ==============================================================

#' @title Perform MSstatsPTM-Style PTM Differential Analysis with Protein-Level Correction
#' @details For each UniProt protein, fits separate linear models at the protein level
#'   and at each PTM peptide level, then computes a protein-corrected PTM log fold-change
#'   (\code{logFC_ptm - logFC_prot}) with Welch-Satterthwaite degrees of freedom.
#'   Adjusted p-values are computed via Benjamini-Hochberg across all PTM-protein pairs.
#'
#'   Note: currently iterates over the first 50 unique protein IDs (debug limit on line 43).
#' @param comparison Named list with elements \code{name}, \code{exp}, and \code{ctrl}
#'   describing the contrast (same structure used elsewhere in the pipeline).
#' @param exp_samples Character vector. Sample names belonging to the experimental group.
#' @param ctrl_samples Character vector. Sample names belonging to the control group.
#' @param ptm_imp_matrix Data frame in long format with columns \code{peptide_id},
#'   \code{PG.UniProtIds}, \code{sample}, and \code{intensity}.
#' @param protein_imp_matrix Data frame in long format with columns
#'   \code{PG.ProteinAccessions}, \code{sample}, and \code{intensity}.
#' @return Data frame with one row per (protein, PTM peptide) pair and columns:
#'   \code{uniprot_id}, \code{peptide_id}, \code{logFC_ptm}, \code{logFC_prot},
#'   \code{logFC} (corrected), \code{AveExpr}, \code{SE_ptm}, \code{SE_prot},
#'   \code{DF_ptm}, \code{DF_prot}, \code{t}, \code{P.Value}, \code{adj.P.Val}.
perform_multitest_ptm_ida_analysis2 <- function(comparison, exp_samples, ctrl_samples,
                                              ptm_imp_matrix, protein_imp_matrix){
  results <- list()
  uniq_uniprot_ids <- unique(protein_imp_matrix$PG.ProteinAccessions)
  # uniq_uniprot_ids <- head(uniq_uniprot_ids, 500)  # DEBUG: limit to 500
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
    curr_prot_stats <- calculate_da_ttest_fast(curr_prots.exp$intensity, curr_prots.ctrl$intensity)

    # get the ptms of the current protein accession and pre-split by peptide_id
    curr_ptms <- ptm_imp_matrix %>% filter(PG.UniProtIds == uniprot_id)
    ptms_by_id <- split(curr_ptms, curr_ptms$peptide_id)

    # cycle through the unique ptms
    for (ptm in names(ptms_by_id)){
      ptm_rows <- ptms_by_id[[ptm]]

      # extract the ptm subset
      curr_ptms.exp  <- ptm_rows[ptm_rows$sample %in% exp_samples,]
      curr_ptms.ctrl <- ptm_rows[ptm_rows$sample %in% ctrl_samples,]

      # calculating the PTM-level stats
      curr_ptm_stats <- calculate_da_ttest_fast(curr_ptms.exp$intensity, curr_ptms.ctrl$intensity)

      # perform MSstatsPTM-style adjusted test
      logfc <- curr_ptm_stats$logFC - curr_prot_stats$logFC

      s2_ptm  <- curr_ptm_stats$SE ^ 2
      s2_prot <- curr_prot_stats$SE ^ 2

      stderr <- sqrt(s2_ptm + s2_prot)
      numer <- (s2_ptm + s2_prot) ^ 2
      denom <- (s2_ptm ^ 2 / curr_ptm_stats$DF) + (s2_prot ^ 2 / curr_prot_stats$DF)
      df <- numer / denom

      t <- logfc / stderr
      pval <- 2 * pt(abs(t), df, lower.tail = FALSE)
      ave_expr <- mean(c(curr_ptms.exp$intensity, curr_ptms.ctrl$intensity), na.rm = TRUE)

      results[[glue("{uniprot_id}__{ptm}")]] <- data.frame(
        uniprot_id = uniprot_id,
        peptide_id = ptm,
        logFC_ptm = curr_ptm_stats$logFC,
        logFC_prot = curr_prot_stats$logFC,
        logFC = logfc,
        AveExpr = ave_expr,
        SE_ptm = curr_ptm_stats$SE,
        SE_prot = curr_prot_stats$SE,
        DF_ptm = curr_ptm_stats$DF,
        DF_prot = curr_prot_stats$DF,
        t = t,
        P.Value = pval
      )
    }
  }

  results_df <- bind_rows(results) %>% mutate(adj.P.Val = p.adjust(P.Value, method = "BH"))
  return(results_df)
}

#' @title Run Full PTM Integration Analysis for a Single Comparison
#' @details Orchestrates the end-to-end PTM analysis pipeline for one contrast:
#'   \enumerate{
#'     \item Filters PTM and protein matrices to shared, unambiguous UniProt IDs.
#'     \item Pivots both matrices to long format.
#'     \item Calls \code{\link{perform_multitest_ptm_ida_analysis}} for protein-corrected
#'       differential analysis.
#'     \item Optionally joins peptide metadata and BioMart gene symbols.
#'     \item Optionally annotates results with imputation categories.
#'     \item Saves differential results to \code{out_dirs$de_data}.
#'     \item Generates volcano, MA, and heatmap figures.
#'     \item Optionally runs GSEA via \code{\link{process_gsea}} and saves results.
#'   }
#' @param comparison Named list with elements \code{name}, \code{exp}, and \code{ctrl}.
#' @param limma_params Named list with elements \code{E} (intensity matrix) and
#'   \code{design} (model matrix), used to identify sample group membership.
#' @param ptm_imp_matrix Numeric matrix or data frame (PTM peptides x samples) of
#'   imputed log2 intensities, with columns \code{peptide_id} and \code{PG.UniProtIds}.
#' @param protein_imp_matrix Numeric matrix or data frame (proteins x samples) of
#'   imputed log2 intensities, with column \code{PG.ProteinAccessions}.
#' @param out_dirs Named list of output directories as returned by
#'   \code{\link{setup_directories}}.
#' @param intensity_matrix_raw Numeric matrix or data frame of pre-imputation intensities
#'   (with \code{NA}s) used to compute imputation categories. Pass \code{NULL} to skip.
#' @param peptide_metadata Data frame with at minimum \code{peptide_id} and
#'   \code{PG.UniProtIds} columns for joining annotation. Pass \code{NULL} to skip.
#' @param ont_option Character. GO ontology for GSEA: \code{"BP"} (default), \code{"MF"},
#'   or \code{"CC"}.
#' @param ensembl A \code{Mart} object from \code{biomaRt::useMart} used for gene symbol
#'   annotation. Pass \code{NULL} to skip BioMart lookup.
#' @param skip_gsea Logical. If \code{TRUE}, GSEA is skipped entirely. Default
#'   \code{FALSE}.
#' @param gsea_method Character. Which GSEA method to run: \code{"go"} (default)
#'   or \code{"msigdb"}. Only one method runs per comparison.
#' @param msigdb_species Character. Species name as expected by \code{msigdbr}
#'   (e.g. \code{"Homo sapiens"} or \code{"Mus musculus"}), used only when
#'   \code{gsea_method = "msigdb"}.
#' @param msigdb_collection Character. MSigDB collection code (e.g. \code{"H"}),
#'   used only when \code{gsea_method = "msigdb"}. Default \code{"H"}.
#' @param msigdb_subcollection Character. MSigDB subcollection (e.g. \code{"CP:KEGG"}),
#'   used only when \code{gsea_method = "msigdb"}. Default \code{NULL}.
#' @param heatmap_norm Character. \code{"zscore"} (default) or \code{"intensity"} for
#'   the heatmap colour scale.
#' @param color1 Character. Primary plot colour (e.g. down-regulated / imputed).
#'   Default \code{"#D55E00"}.
#' @param color2 Character. Secondary plot colour (e.g. up-regulated / observed).
#'   Default \code{"#0072B2"}.
#' @return Named list with elements:
#'   \describe{
#'     \item{limma}{Data frame of differential analysis results.}
#'     \item{gsea}{\code{gseaResult} object or \code{NULL} if GSEA was skipped or failed.}
#'     \item{protein_counts}{Data frame of per-protein PTM counts used for GSEA, or
#'       \code{NULL}.}
#'     \item{highlighted_ids}{Character vector of peptide IDs highlighted on the volcano
#'       plot.}
#'   }
#'   Returns \code{NULL} if the analysis fails entirely.
run_single_ptm_ida_site_level_comparison <- function(comparison, limma_params, ptm_imp_matrix,
                                          protein_imp_matrix, out_dirs, intensity_matrix_raw = NULL,
                                          peptide_metadata = NULL, ont_option = "BP", ensembl = NULL, org_db = NULL,
                                          skip_gsea = FALSE, gsea_method = "go", msigdb_species = "Homo sapiens",
                                          msigdb_collection = "H", msigdb_subcollection = NULL, heatmap_norm = "zscore",
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

    diff_results <- perform_multitest_ptm_ida_analysis2(comparison, exp_samples, ctrl_samples, ptm_imp_matrix_filtered_long, protein_imp_matrix_filtered_long)
    
    # Join peptide metadata to get UniProt IDs, then annotate with gene symbols
    if (!is.null(peptide_metadata)) {

      diff_results <- diff_results %>% left_join(peptide_metadata, by="peptide_id") %>%
        mutate(uniprot_id = sapply(strsplit(PG.UniProtIds, ";"), `[`, 1))
      
      annotation <- tryCatch({
        m <- getBM(
          attributes = c("uniprotswissprot", "hgnc_symbol"),
          filters = "uniprotswissprot",
          values = unique(na.omit(diff_results$uniprot_id)),
          mart = ensembl
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
      na_exp <- rowSums(is.na(intensity_matrix_raw[diff_results$peptide_id, exp_samples,  drop = FALSE]))
      na_ctrl <- rowSums(is.na(intensity_matrix_raw[diff_results$peptide_id, ctrl_samples, drop = FALSE]))
      na_total <- na_exp + na_ctrl
      
      diff_results <- diff_results %>%
        mutate(
          imputation_category = case_when(
            adj.P.Val >= 0.05 ~ "not-significant",
            na_total == 0 ~ "complete-data",
            (na_exp == n_exp & na_ctrl == 0) | (na_ctrl == n_ctrl & na_exp == 0) ~ "on-off",
            na_total == 1 ~ "imputation-low",
            na_total == 2 ~ "imputation-medium",
            na_total >= 3 ~ "imputation-high",
            TRUE ~ "other"
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
    
    # generating the PTM volcano plot
    flog.info("Generating PTM volcano plot: %s", comparison$name)
    volcano_result <- generate_volcano_ptm(diff_results, comparison$exp, comparison$ctrl, 
                                          logFC_colname="logFC", color1 = color1, color2 = color2)
    save_plot(volcano_result$plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"),
              width = 10, height = 8)
    
    # generating the PMT MA plot
    flog.info("Generating PTM MA plot: %s", comparison$name)
    ma_plot <- generate_ma_plot_ptm(diff_results, comparison$exp, comparison$ctrl,
                                    highlighted_ids = volcano_result$highlighted_ids,
                                    color1 = color1, color2 = color2)
    save_plot(ma_plot, create_file_path(out_dirs$ma, "", comparison$name, "_ma.png"),
              width = 10, height = 8)
    
    # Generating the PTM heatmap
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

    # Build KSEA input: required columns are Protein, Gene, Peptide, Residue.Both, p, FC.
    # Residue.Both is the modification site code (e.g. "C55") parsed from the second
    # underscore-delimited field of peptide_id (format: <Gene>_<Residue>_..._<Sequence>).
    # Peptide is set to "NULL" per the KSEA spec when the sequence is unavailable.
    # FC is back-transformed from log2 space (2^logFC); p uses the BH-adjusted p-value.
    flog.info("Saving KSEA input data: %s", comparison$name)
    ksea_input <- diff_results %>%
      mutate(
        Residue.Both = sub("^[^_]+_([A-Za-z][0-9]+)_.*$", "\\1", peptide_id),
        FC = 2^logFC
      ) %>%
      select(
        Protein      = uniprot_id,
        Gene         = hgnc_symbol,
        Peptide      = peptide_id,
        Residue.Both,
        p            = adj.P.Val,
        FC
      ) %>%
      mutate(Peptide = "NULL")
    write.csv(ksea_input,
              create_file_path(out_dirs$ksea, comparison$name, "_ksea_input_data"),
              row.names = FALSE)

    # # Run KSEAapp: all functions, outputs saved to out_dirs$ksea
    # flog.info("Running KSEAapp for %s", comparison$name)
    # tryCatch({
    #   # data("PSP&NetworKIN_Kinase_Substrate_Dataset_July2016", package = "KSEAapp",
    #   #      envir = environment())
    #   # KSData <- get("PSP&NetworKIN_Kinase_Substrate_Dataset_July2016",
    #   #               envir = environment())

    #   old_wd <- setwd(out_dirs$ksea)
    #   on.exit(setwd(old_wd), add = TRUE)

    #   # generate ks 
    #   ksea_ks_table <- KSEAapp::KSEA.KS_table(KSData, ksea_input, NetworKIN = FALSE)
    #   write.csv(ksea_ks_table,
    #             create_file_path(out_dirs$ksea, comparison$name, "_ksea_ks_table"),
    #             row.names = FALSE)

    #   # generate ks scores
    #   ksea_scores <- KSEAapp::KSEA.Scores(KSData, ksea_input, NetworKIN = FALSE)
    #   write.csv(ksea_scores,
    #             create_file_path(out_dirs$ksea, comparison$name, "_ksea_scores"),
    #             row.names = FALSE)

    #   # generate barplot
    #   KSEAapp::KSEA.Barplot(KSData, ksea_input,
    #                         NetworKIN = FALSE,
    #                         m.cutoff = 5, p.cutoff = 0.05,
    #                         export = TRUE)
    # }, error = function(e) {
    #   flog.error("KSEAapp failed for '%s': %s", comparison$name, e$message)
    # })

    # conduct GSEA
    if (skip_gsea == TRUE) {
      flog.info("Skipping GSEA for %s (--skip-gsea flag set)", comparison$name)
      gse <- NULL
      protein_counts <- NULL
    } else {
      agg_result <- tryCatch(
        aggregate_ptm_for_gsea(diff_results, peptide_metadata, ensembl),
        error = function(e) {
          flog.error("PTM GSEA aggregation failed for '%s': %s",
                     comparison$name, e$message)
          NULL
        }
      )
      protein_counts <- if (!is.null(agg_result)) agg_result$counts else NULL

      if (gsea_method == "msigdb") {
        flog.info("Running PTM MSigDB GSEA for %s", comparison$name)
        gse <- tryCatch({
          if (!is.null(agg_result)) process_gsea_msigdb(agg_result$data, species = msigdb_species,
                                                          collection = msigdb_collection, subcollection = msigdb_subcollection) else NULL
        }, error = function(e) {
          flog.error("PTM MSigDB GSEA failed for '%s': %s", comparison$name, e$message)
          NULL
        })
        if (!is.null(gse)) {
          flog.info("PTM MSigDB GSEA returned results for %s", comparison$name)
          tryCatch({
            write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, comparison$name, "_msigdb_analysis.csv"), row.names = FALSE, quote = FALSE)
            gsea_plot <- create_barplot(gse, create_comparison_name(comparison$ctrl, comparison$exp, "MSigDB GSEA "), color1 = color1, color2 = color2)
            save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                      width = 10, height = 12)
          }, error = function(e) {
            flog.error("Failed to save PTM MSigDB GSEA outputs for '%s': %s",
                       comparison$name, e$message)
          })
        }
      } else {
        flog.info("Running PTM GSEA for %s", comparison$name)
        gse <- tryCatch({
          if (!is.null(agg_result)) process_gsea(agg_result$data, org_db, ont_option = ont_option) else NULL
        }, error = function(e) {
          flog.error("PTM GSEA failed for '%s': %s", comparison$name, e$message)
          NULL
        })
        if (!is.null(gse)) {
          flog.info("PTM GSEA returned results for %s", comparison$name)
          tryCatch({
            write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, comparison$name, "_go_analysis.csv"), row.names = FALSE, quote = FALSE)
            gsea_plot <- create_barplot(gse, create_comparison_name(comparison$ctrl, comparison$exp, "GSEA "), color1 = color1, color2 = color2)
            save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_gsea.png"),
                      width = 10, height = 12)
          }, error = function(e) {
            flog.error("Failed to save PTM GSEA outputs for '%s': %s",
                       comparison$name, e$message)
          })
        }
      }
    }

    return(list(limma = diff_results, gsea = gse, protein_counts = protein_counts, highlighted_ids = volcano_result$highlighted_ids))

  }, error = function(e) {
    flog.error("run_analysis_ptm_integrated failed for comparison '%s': %s\n%s",
               comparison$name, e$message, paste(deparse(e$call), collapse = " "))
    return(NULL)
  })
}





