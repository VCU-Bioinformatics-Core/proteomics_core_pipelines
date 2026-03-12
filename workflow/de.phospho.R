# ==========================
# limma differential expression pipeline for phospopeptides
# ==========================
# In contrast to the protein-level analysis performed in de.regular.R,
# this workflow is specifically designed for phosphopeptide-level analysis
# and is parameterized accordingly.

library(here)
library(dplyr)
library(data.table)
library(tidyverse)
library(janitor)
library(scales)
library(clusterProfiler)
library(enrichplot)
library(readr)
library(DT)
library(limma)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(RColorBrewer)
library(purrr)
library(plotly)
library(stats)
library(orca)
library(reticulate)
library(optparse)
library(htmlwidgets)
library(ComplexHeatmap)
library(circlize)
library(DEP)

# Connect to Ensembl (use the dataset for your species)
library(biomaRt)

# source helper functions
source(here::here('workflow/helpers.R'))
source(here::here("workflow/report_generator.phospho.R"))
#debug(generate_volcano)
#debug(run_analysis)
#debug(generate_heatmap)
#debug(process_gsea)

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
              help = "Genome for annotation: 'mouse' or 'human' [default= %default]"),
  make_option(c("-i", "--imputation"), type = "character", default = "MinProb",
              help = "DEP imputation method: 'MinProb', 'knn', 'bpca', 'QRILC', 'man', or 'none' [default= %default]"),
  make_option(c("-q", "--imputation-q"), type = "double", default = 0.01,
              help = "q parameter for MinProb/QRILC: quantile cutoff for the left-censored distribution [default= %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility of stochastic imputation methods [default= %default]"),
  make_option(c("--heatmap-top-n"), type = "integer", default = 1000,
              help = "Number of top molecules by CV to show in the global heatmap [default= %default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

runID <- opt$runid
countData <- opt$counts
samplesheet <- opt$samplesheet
outDir <- opt$outdir
genome <- opt$annotation
imputation_method <- opt$imputation
imputation_q <- opt$`imputation-q`
imputation_seed <- opt$seed
heatmap_top_n <- opt$`heatmap-top-n`

if (is.null(runID) || is.null(countData) || is.null(samplesheet)) {
  print_help(opt_parser)
  stop("--runid, --counts, and --samplesheet are required.")
}

# ==========================
# Debug mode
# ==========================

#debug <- TRUE
debug <- FALSE
if (debug){
  # runID <- 'Test'
  # countData <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.csv'
  # samplesheet <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/raw/20251001_U266_UT_VS_S+B_LFQ_Report.samplesheet.csv'
  # outDir <- '/global/projects/proteomics_core/analyst_workspace/pipeline_test_data/analyses/'
  # annotation <- 'human'
  # genome <- 'human'
  
  runID <- 'Test'
  countData <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Report.CSV/Phospho_Proteome_DEV/U266_LFQ_phospho_4condition_DEV_Report.Complete.csv'
  samplesheet <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Bioinformatics/results/phospho_proteome/samplesheet.csv'
  outDir <- '/global/projects/proteomics_core/researcher_workspace/steven_grant/Bioinformatics/results/phospho_proteome/'
  annotation <- 'human'
  genome <- 'human'
}

# ==========================
# Load annotation DB
# ==========================
if (genome == "human") {
  if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  annotation_db <- org.Hs.eg.db
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

} else if (genome == "mouse") {
  if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
  annotation_db <- org.Mm.eg.db
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
} else {
  stop("Invalid genome specified. Use 'mouse' or 'human'")
}

# ==========================
# Read and prepare data
# ==========================
out_dirs <- setup_directories(outDir)
full_peptide_levels <- data.frame(read_csv(countData, col_names = TRUE, na = c("", "NA", "Filtered")))

print_peptide_report <- function(df, prefix){
  num_genes <- nrow(df)
  msg <- glue("{prefix}: {num_genes}")
  print(msg)
}
print_peptide_report(full_peptide_levels, 'Peptides - Total Amount')

# remove pepties with missing names (deprecated)
# full_peptide_levels <- full_peptide_levels %>% filter(!is.na(PG.Genes))

# remove peptides that come from the cRAP database
full_peptide_levels <- full_peptide_levels %>% filter(!grepl("cRAP[0-9]+", PG.ProteinAccessions))
print_peptide_report(full_peptide_levels, 'Peptides - without cRAP')

# remove peptides that are not phosphorylated
full_peptide_levels <- full_peptide_levels %>% filter(grepl("Phospho \\(STY\\)", EG.PrecursorId))
print_peptide_report(full_peptide_levels, 'Peptides - only phospho')

# name the rownames with the peptide
rownames(full_peptide_levels) <- paste0(full_peptide_levels$PG.Genes, ' -- ', full_peptide_levels$EG.PrecursorId)

# Extract a dataframe with just peptide levels
#non_numeric_cols <- c('PG.ProteinAccessions', 'PG.Genes', 'PG.UniProtIds', 'PG.ProteinNames', 'PG.FastaHeaders', 'EG.PrecursorId')
non_numeric_cols <- c('PG.ProteinAccessions', 'PG.Genes', 'PG.UniProtIds', 'PG.FASTAName', 'EG.PrecursorId')
counts <- full_peptide_levels %>% dplyr::select(-any_of(non_numeric_cols))

# Metadata table for GSEA ID mapping (peptide_id -> PG.UniProtIds)
peptide_metadata <- full_peptide_levels %>%
  transmute(peptide_id = rownames(full_peptide_levels), PG.UniProtIds)

# Collect the comparisons from the samplesheet
comparisons_raw <- read.csv(samplesheet)
comparisons_cols <- colnames(comparisons_raw)[3:ncol(comparisons_raw)]
comparisons <- list()

for (comparisons_name in comparisons_cols) {

  comparisons_vector <- comparisons_raw[[comparisons_name]]
  
  exp_group <- unique(comparisons_raw$GroupID[comparisons_vector == 1 & !is.na(comparisons_vector)])
  ctrl_group <- unique(comparisons_raw$GroupID[comparisons_vector == 0 & !is.na(comparisons_vector)])
  
  if (length(exp_group) > 0 && length(ctrl_group) > 0) {
    comparisons[[length(comparisons) + 1]] <- list(
      name = comparisons_name,
      exp = exp_group,
      ctrl = ctrl_group
    )
  }
}

# Optional log2-transform (if values are not log-scaled)
intensity_matrix <- log2(counts + 1)
intensity_matrix_raw <- intensity_matrix  # snapshot before imputation

# ==========================
# Imputation with DEP
# ==========================
if (imputation_method != "none") {
  set.seed(imputation_seed)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = as.matrix(intensity_matrix)),
    colData = S4Vectors::DataFrame(
      label     = comparisons_raw$SampleID,
      condition = comparisons_raw$GroupID
    ),
    rowData = S4Vectors::DataFrame(
      name = rownames(intensity_matrix),
      ID   = rownames(intensity_matrix)
    )
  )
  if (imputation_method %in% c("MinProb", "QRILC")) {
    se_imputed <- DEP::impute(se, fun = imputation_method, q = imputation_q)
  } else {
    se_imputed <- DEP::impute(se, fun = imputation_method)
  }
  intensity_matrix <- as.data.frame(SummarizedExperiment::assay(se_imputed))
  print(glue("Imputation applied: {imputation_method} (q={imputation_q})"))
} else {
  print("Imputation skipped (--imputation none)")
}

# ==========================
# Set up the limma design matrix
# ==========================
sample_info <- data.frame(sample = comparisons_raw$SampleID,
                          condition = comparisons_raw$GroupID)
design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(factor(sample_info$condition))
limma_params <- list(E = intensity_matrix, design = design)

# ==========================
# Analysis loop
# ==========================
results <- vector("list", length(comparisons))
for (i in seq_along(comparisons)) {

  # may need to consider what to do about imputations, especially dealing
  # with the limma_param
  
  print(glue('Analysis {i}'))
  
  # run the analysis on the current samples
  curr_result <- run_analysis_phospho(comparisons[[i]], limma_params, intensity_matrix, out_dirs, intensity_matrix_raw, peptide_metadata)
  
  # save the current results if successful
  if (!is.null(curr_result)){
    results[[i]] <- curr_result
  }
}

# ==========================
# Principal Component Analysis
# ==========================
print("# Principal Component Analysis")
input_pca_matrix <- as.matrix(intensity_matrix)

# Values already imputed via DEP above; transpose so samples are rows for prcomp
input_pca_matrix <- t(input_pca_matrix)

# calculate PCA results and scores 
pca_results <- prcomp(input_pca_matrix, rank. = 3)
pc_scores <- pca_results$x

# calculate the variance explained from PCA
pca_var_pct <- round((pca_results$sdev^2)/sum(pca_results$sdev^2)*100, 1)

# create a dataframe to plot PC1 and PC2
two_comps_pca_df <- data.frame(Sample=rownames(pc_scores),
                     X=pc_scores[,1],
                     Y=pc_scores[,2],
                     Group=sample_info$condition)

# plot PC1 versus PC2 values
pca_plot <- ggplot(data=two_comps_pca_df, aes(x=X, y=Y, label=Sample, color=Group)) +
  geom_point(size = 1) +
  # geom_text_repel(aes(label = as.character(Sample)), show.legend = FALSE, size = 2.5) +
  #geom_text(aes(label = as.character(Sample)), show.legend = FALSE, size = 2.5) +
  xlab(paste("PC1 - ", pca_var_pct[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca_var_pct[2], "%", sep=""))
  #theme_bw()
print(pca_plot)
ggsave(pca_plot, filename = str_c(out_dirs$pca, "/PCA_plot.png"))

# # Interactive PCA 2D/3D # removing for proteomic core analyses
# print("# Interactive PCA 2D/3D")
# fig <- plot_ly(two_comps_pca_df, x = ~X, y = ~Y, color = ~Group,
#                colors = c('steelblue','firebrick','olivedrab', 'plum'),
#                type = 'scatter', mode = 'markers', size = 5.7)
# fig <- fig %>% layout(plot_bgcolor='#e5ecf6')
# export_plotly_to_html(fig, paste0(out_dirs$pca, "/allsamples_PCA_plot.html"))
# 
# components <- data.frame(pca_results$x)
# components$PC2 <- -components$PC2
# components$PC3 <- -components$PC3
# components <- cbind(components, Group=pca_df$Group)
# tot_explained_variance_ratio <- 100 * sum(summary(pca_results)[["importance"]]['Proportion of Variance',])
# fig3D <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group,
#                  colors = c('steelblue','firebrick','olivedrab', 'plum')) %>%
#   add_markers(size = 12) %>%
#   layout(title = paste('Total Explained Variance = ', tot_explained_variance_ratio),
#          scene = list(bgcolor = "#e5ecf6"))
# export_plotly_to_html(fig3D, paste0(out_dirs$pca, "/allsamples_PCA_plot3D.html"))

# ==========================
# Generate global heatmap
# ==========================
print('Generating global heatmap')
generate_global_heatmap(intensity_matrix, out_dirs, top_n = heatmap_top_n,
                        molecule_label = "Phosphopeptides")

# ==========================
# Save RDS
# ==========================
print('Save RDS')
#rds <- list(results, comparisons, out_dirs, pca_plot, fig, fig3D)
imputation_params <- list(method = imputation_method, q = imputation_q)
rds <- list(results, comparisons, out_dirs, pca_plot, intensity_matrix_raw, intensity_matrix, imputation_params, sample_info)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#rds_path <- glue("analysis_results_{timestamp}.rds")
rds_path <- file.path(outDir, "data/analysis_results.rds")
saveRDS(rds, rds_path)

# ==========================
# Generate the HTML report
# ==========================
generate_report(rds_path, output_dir = outDir)
