# ==========================
# limma-voom Differential Expression Pipeline
# ==========================

# Package installation and loading
options(repos = c(CRAN = "https://cran.r-project.org"))

if (!requireNamespace("BiocManager")) install.packages("BiocManager")
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  here, dplyr, data.table, tidyverse, janitor, scales, ggrepel,
  clusterProfiler, enrichplot, readr, DT, edgeR, limma,
  ggplot2, AnnotationDbi, gplots, RColorBrewer, purrr,
  plotly, stats, orca, reticulate, optparse, htmlwidgets
)

# Connect to Ensembl (use the dataset for your species)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ==========================
# Command-line options
# ==========================
option_list = list(
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Required. Path to merged counts.tsv file"),
  make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
              help = "Required. Path to samplesheet.csv file"),
  make_option(c("-o", "--outdir"), type = "character", default = "./output",
              help = "Output directory [default= %default]"),
  make_option(c("-r", "--runid"), type = "character", default = NULL,
              help = "Required. Unique run ID"),
  make_option(c("-a", "--annotation"), type = "character", default = "mouse",
              help = "Genome for annotation: 'mouse' or 'human' [default= %default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

runID <- opt$runid
countData <- opt$counts
samplesheetData <- opt$samplesheet
outDir <- opt$outdir
genome <- opt$annotation

# ==========================
# Debug mode
# ==========================
# debug <- TRUE
# if (debug){
#   runID <- "test_run"
#   countData <- "./example_counts.tsv"
#   samplesheetData <- "./example_samplesheet.csv"
#   outDir <- "./test_output"
#   genome <- "mouse"
# }

debug <- TRUE
if (debug){
  runID <- 'Test'
  countData <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.csv'
  samplesheetData <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.samplesheet.csv'
  outDir <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/analysis/'
  annotation <- 'human'
  genome <- 'human'
}

# ==========================
# Load annotation DB
# ==========================
if (genome == "human") {
  if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  annotation_db <- org.Hs.eg.db
} else if (genome == "mouse") {
  if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
  annotation_db <- org.Mm.eg.db
} else {
  stop("Invalid genome specified. Use 'mouse' or 'human'")
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

# ==========================
# limma-voom Analysis Function
# ==========================
perform_limma_analysis <- function(v, exp, ctrl) {
  tryCatch({
    print(paste("Performing limma-voom analysis for", exp, "vs", ctrl))
    
    # Build contrast
    contrast_matrix <- makeContrasts(contrasts = paste0(exp, "-", ctrl),
                                     levels = colnames(v$design))
    
    # Fit linear model
    fit <- lmFit(v$E, v$design)
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
    # gene_list <- gene_list[!is.na(gene_list)]
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
# Output directories
# ==========================
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
# Read and prepare data
# ==========================
out_dirs <- setup_directories(outDir)
full_prot_levels <- data.frame(read_csv(countData, col_names = TRUE))
full_prot_levels <- full_prot_levels %>% filter(!is.na(PG.Genes)) # remove prots with missing names
full_prot_levels <- full_prot_levels %>% distinct(PG.ProteinGroups, .keep_all = TRUE) # drop duplicate (temporarily)
colnames(full_prot_levels) <- sub("^X\\.\\d+\\.\\.", "", colnames(full_prot_levels)) # Remove the "X.#.." prefix only from columns that start with "X."

# Extract all of the uniprotswiss ids
uniprot_raw <- full_prot_levels$PG.ProteinGroups
uniprot_clean <- unlist(strsplit(uniprot_raw, ";"))
uniprot_clean <- unique(uniprot_clean)

# Query the uniprot ids to obtain the ensembl ids
mapping <- getBM(
  attributes = c("uniprotswissprot", "ensembl_gene_id"), # "external_gene_name"
  filters = "uniprotswissprot",
  values = uniprot_clean,
  mart = ensembl
)
mapping <- mapping %>% distinct(uniprotswissprot, .keep_all = TRUE)

# Merge the Ensemble ID to the full prot df
full_prot_levels <- full_prot_levels %>% left_join(mapping, by=join_by(PG.ProteinGroups == uniprotswissprot))
rownames(full_prot_levels) <- full_prot_levels$PG.ProteinGroups # use gene name as the row name


# Extract a dataframe with just protein levels
counts <- full_prot_levels %>% dplyr::select(where(is.numeric)) 
counts <- counts %>% dplyr::select(-PG.Pvalue, -PG.Qvalue)

# Load the samplesheet
contrasts_raw <- read.csv(samplesheetData)
contrast_cols <- colnames(contrasts_raw)[3:ncol(contrasts_raw)]
comparisons <- list()

for(i in seq_along(contrast_cols)) {
  
  print(i)
  
  contrast_name <- contrast_cols[i]
  contrast_data <- contrasts_raw[[contrast_name]]
  
  exp_group <- unique(contrasts_raw$GroupID[contrast_data == 1 & !is.na(contrast_data)])
  ctrl_group <- unique(contrasts_raw$GroupID[contrast_data == 0 & !is.na(contrast_data)])
  
  if (length(exp_group) > 0 && length(ctrl_group) > 0) {
    comparisons[[length(comparisons) + 1]] <- list(
      name = contrast_name,
      exp = exp_group,
      ctrl = ctrl_group
    )
  }
}

countsdf <- counts %>% mutate(across(where(is.numeric), ~ as.integer(round(.))))

# Optional log2-transform (if values are not log-scaled)
intensity_matrix <- countsdf
if (max(intensity_matrix, na.rm = TRUE) > 100) {  # crude check for raw vs log
  intensity_matrix <- log2(intensity_matrix + 1)
}

# ==========================
# voom transformation
# ==========================
sample_info <- data.frame(sample = contrasts_raw$SampleID,
                          condition = contrasts_raw$GroupID)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(factor(sample_info$condition))
v <- list(E = intensity_matrix, design = design)  # mimic voom object

# ==========================
# Analysis loop
# ==========================
run_analysis <- function(comparison, v, normalized_counts, out_dirs) {
  tryCatch({
    print(paste("\nStarting analysis for comparison:", comparison$name))
    
    print('Running limma')
    limma_results <- perform_limma_analysis(v, comparison$exp, comparison$ctrl)
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
    #gse <- process_gsea(annotated_results)
    
    
    
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

results <- list()
for (i in seq_along(comparisons)) {
  result <- run_analysis(comparisons[[i]], v, intensity_matrix, out_dirs)
  if (!is.null(result)) results[[i]] <- result
}

# ==========================
# PCA (same as before)
# ==========================
print("# PCA (same as before)")
data_reduced_zero <- as.data.frame(intensity_matrix) %>% filter(if_any(where(is.numeric)))
pca_matrix <- as.matrix(data_reduced_zero)

pca_matrix <- apply(pca_matrix, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})


data_pca <- prcomp(pca_matrix)
pc_scores <- data_pca$x
# rownames(pc_scores) <- colnames(data_reduced_zero)
pca_var_per <- round((data_pca$sdev^2)/sum(data_pca$sdev^2)*100, 1)

pca_df <- data.frame(Sample=rownames(pc_scores), # KEEP DEBUGGING HERE
                     X=pc_scores[,1],
                     Y=pc_scores[,2],
                     Group=sample_info$condition)

pca_plot <- ggplot(data=pca_df, aes(x=X, y=Y, label=Sample, color=Group)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
  theme_bw()

ggsave(pca_plot, filename = str_c(out_dirs$pca, "/PCA_plot.png"))

# Interactive PCA 2D/3D
print("# Interactive PCA 2D/3D")
fig <- plot_ly(pca_df, x = ~X, y = ~Y, color = ~Group,
               colors = c('steelblue','firebrick','olivedrab', 'plum'),
               type = 'scatter', mode = 'markers', size = 5.7) %>%
  layout(plot_bgcolor='#e5ecf6')

export_plotly_to_html(fig, paste0(out_dirs$pca, "/allsamples_PCA_plot.html"))

prin_comp <- prcomp(pca_matrix, rank. = 3)
components <- data.frame(prin_comp$x)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components <- cbind(components, Group=pca_df$Group)
tot_explained_variance_ratio <- 100 * sum(summary(prin_comp)[["importance"]]['Proportion of Variance',])
fig3D <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group,
                 colors = c('steelblue','firebrick','olivedrab', 'plum')) %>%
  add_markers(size = 12) %>%
  layout(title = paste('Total Explained Variance = ', tot_explained_variance_ratio),
         scene = list(bgcolor = "#e5ecf6"))
export_plotly_to_html(fig3D, paste0(out_dirs$pca, "/allsamples_PCA_plot3D.html"))

# ==========================
# Save RDS
# ==========================
print('Save RDS')
rds <- list(results, comparisons, out_dirs, pca_plot, fig, fig3D)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_name <- paste0("analysis_results", "_", timestamp, ".rds")
saveRDS(rds, rds_name)
