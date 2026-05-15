#' @title Run Gene Set Enrichment Analysis (GSEA) via gseGO
#' @details Builds a ranked gene list from \code{result$logFC} (keyed by
#'   \code{result$ensembl_gene_id}), drops genes with near-zero fold-change
#'   (\code{|logFC| < 0.1}), and runs \code{clusterProfiler::gseGO}. Returns
#'   \code{NULL} (with a warning) if there are fewer than two genes or no
#'   enriched terms are found.
#' @param result Data frame of limma results containing at minimum \code{logFC} and
#'   \code{ensembl_gene_id} columns.
#' @param p Numeric. Adjusted p-value threshold (filtering of \code{result} is
#'   expected upstream). Default \code{0.05}.
#' @param lfc Numeric. Log2 fold-change threshold (passed through). Default \code{0.58}.
#' @param ont_option Character. GO ontology to test: \code{"BP"} (default), \code{"MF"},
#'   or \code{"CC"}.
#' @return A \code{gseaResult} object with gene symbols resolved via
#'   \code{setReadable}, or \code{NULL} if GSEA cannot be performed or yields no
#'   significant terms.
process_gsea <- function(result, org_db, p = 0.05, lfc = 0.58, ont_option = "BP") {
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
                        OrgDb = org_db)

    if (is.null(gse_result) || nrow(gse_result@result) == 0) {
      flog.warn("No enriched GO terms found in GSEA result")
      return(NULL)
    }

    setReadable(gse_result, OrgDb = org_db, keyType = "ENSEMBL")
  }, error = function(e) {
    flog.error("GSEA processing failed: %s", e$message)
    return(NULL)
  })
}