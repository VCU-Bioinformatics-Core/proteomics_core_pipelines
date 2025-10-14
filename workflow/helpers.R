# ==========================
# Helper Functions
# ==========================
create_file_path <- function(base_dir, prefix, name, extension = ".csv") {
  file.path(base_dir, paste0(prefix, name, extension))
}

create_comparison_name <- function(exp, ctrl, prefix = "") {
  paste0(prefix, exp, " vs. ", ctrl)
}

save_plot <- function(plot, filename, width = 15, height = 17) {
  ggsave(plot, filename = filename, width = width, height = height, bg = "white")
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
perform_limma_analysis <- function(limma_params, exp, ctrl) {
  tryCatch({
    print(paste("Performing limma-voom analysis for", exp, "vs", ctrl))
    
    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(limma_params$design))
    
    # Fit linear model
    fit <- lmFit(limma_params$E, limma_params$design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Extract results
    res <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    print(paste("Analysis complete. Number of results:", nrow(res)))
    return(as.data.frame(res))
  }, error = function(e) {
    print(paste("Error in limma analysis:", e$message))
    stop(e)
  })
}

generate_volcano <- function(data, exp_name, ctrl_name, p = 0.05, lfc = 0.58, 
                             sig = "adj.P.Val", out_dir = ".") {
  
  print('Labeling the data')
  labeled_dat <- data %>% 
    mutate(
      color_tag = case_when(
        eval(as.symbol(sig)) < p & logFC < -lfc ~ "Under expressed",
        eval(as.symbol(sig)) < p & logFC > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      )
    )
  
  print('Extracting the top 20 genes')
  top_20_genes_up <- labeled_dat %>%
    filter(P.Value <= p & logFC >= lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 20) %>%
    pull(uniprotswissprot)
  
  print('Extracting the bottom 20 genes')
  top_20_genes_dn <- labeled_dat %>%
    filter(P.Value <= p & logFC <= -lfc) %>%
    arrange(eval(as.symbol(sig))) %>%
    slice_head(n = 20) %>%
    pull(uniprotswissprot)
  
  print('Plotting')
  labeled_dat <- labeled_dat %>% mutate(highlight = ifelse(uniprotswissprot %in% c(top_20_genes_up, top_20_genes_dn), uniprotswissprot, NA))
  labeled_dat <- labeled_dat %>% filter(!is.na(logFC), !is.na(P.Value), is.finite(logFC), is.finite(P.Value))
  p <- labeled_dat %>% ggplot(aes(x = logFC, y = -log10(P.Value)),
                              color = color_tag, label = ifelse(highlight == TRUE, uniprotswissprot, NA)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_label_repel(aes(label = highlight), max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = c("firebrick", "steelblue")) +
    geom_hline(yintercept = -log10(p), col = "red", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme(legend.title = element_blank()) +
    labs(x = "Log2 Fold-Change (FC)",
         y = paste0("-log10( P-value )"),
         title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
  return(p)
}

generate_heatmap <- function(results_df, normalized_counts, p = 0.05, lfc = 0.58, 
                             exp_name, ctrl_name, fig_dir) {
  
  # filter the data for certain pvalues and lfc
  filtered_data <- results_df %>% 
    dplyr::filter((adj.P.Val < p & abs(logFC) >= lfc))
  
  # add jitter
  values <- normalized_counts[filtered_data$uniprotswissprot, ] %>%
    as.matrix() %>%
    jitter(factor = 1, amount = 0.00001)
  
  # normalize the values
  zscores <- t(scale(t(values)))
  
  # Replace NaN with row means
  print("WARNING JR: REMOVE THIS REPLACING OF NAN WITH SOMETHING MORE ACCURATE WHEN READY")
  zscores[is.na(zscores) | is.nan(zscores) | is.infinite(zscores)] <- 0
  
  # plot the heatmap
  heatmap.2(zscores,
            col = colorRampPalette(c("blue", "white", "firebrick"))(20),
            density.info = "none",
            dendrogram = "both",
            Colv = TRUE,
            trace = "none",
            margins = c(10, 10),
            labRow = NA)
}

# ==========================
# gsea analysis functions
# ==========================

process_gsea <- function(result, p = 0.05, lfc = 0.58) {
  tryCatch({
    sig_genes <- result %>% arrange(desc(logFC))
    sig_genes <- sig_genes %>% filter(!is.na(logFC))
    
    if(nrow(sig_genes) < 2) {
      message("Not enough significant genes for GSEA analysis")
      return(NULL)
    }
    
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
                        ont = "ALL",
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

run_analysis <- function(comparison, limma_params, normalized_counts, out_dirs) {
  tryCatch({
    print(paste("\nStarting analysis for comparison:", comparison$name))
    
    print('Running limma')
    limma_results <- perform_limma_analysis(limma_params, comparison$exp, comparison$ctrl)
    limma_results <- limma_results %>% rownames_to_column(var = "uniprotswissprot")
    
    print('Annotating the results')
    annotated_results <- limma_results %>% left_join(mapping, by=join_by(uniprotswissprot))
    output_file <- create_file_path(out_dirs$de_data, "limma_", comparison$name)
    write.csv(annotated_results, output_file)
    
    print('Generating the volcano plot')
    volcano_plot <- generate_volcano(annotated_results, comparison$exp, comparison$ctrl)
    save_plot(volcano_plot, create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"))
    
    print('Generating the heatmap')
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 800, height = 1200, res = 150)
    generate_heatmap(limma_results, intensity_matrix,
                     exp_name = comparison$exp, ctrl_name = comparison$ctrl, fig_dir = out_dirs$heatmap)
    dev.off()
    
    print('Running GSEA')
    gse <- NULL # process_gsea(annotated_results)
    gse <- process_gsea(annotated_results)
    
    if(!is.null(gse)) {
      print('GSE has data')
      write.csv(as.data.frame(gse), create_file_path(out_dirs$gsea_data, "GO_Analysis_", comparison$name))
      
      print('create dotplot')
      gsea_plot <- create_dotplot(gse, create_comparison_name(comparison$exp, comparison$ctrl, "GSEA "))
      
      print('save plot')
      save_plot(gsea_plot, create_file_path(out_dirs$gsea, "", comparison$name, "_GSEA.png"))
    }
    
    return(list(limma = annotated_results, gsea = gse))
    
    
  }, error = function(e) {
    message("Error in run_analysis: ", e$message)
    return(NULL)
  })
}















